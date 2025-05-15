
use anyhow::{Context, bail};
use log::debug;
use noodles::core::Position;
use noodles::vcf;
use noodles::vcf::header::record::value::{Map, map};
use noodles::vcf::variant::io::Write;
use noodles_util::variant::io::indexed_reader::Builder as VcfBuilder;
use noodles::vcf::variant::record::samples::keys::key as vcf_key;
use noodles::vcf::variant::record_buf;
use std::collections::BTreeMap;
use std::path::Path;

use crate::data_types::compare_benchmark::CompareBenchmark;
use crate::data_types::compare_region::CompareRegion;
use crate::data_types::phase_enums::PhasedZygosity;
use crate::data_types::variant_metrics::{VariantMetrics, VariantSource};
use crate::data_types::variants::Variant;

pub const BENCHMARK_DECISION_KEY: &str = "BD";
pub const EXPECTED_ALLELE_KEY: &str = "EA";
pub const OBSERVED_ALLELE_KEY: &str = "OA";
pub const REGION_ID: &str = "RI";

/// Wrapper struct that will parse the results from the core processes and determine where to write variants
pub struct VariantCategorizer {
    /// Header that goes into the truth VCFs
    truth_header: vcf::Header,
    /// Header that goes into the query VCFs
    query_header: vcf::Header,
    /// Contains links to all of the various VCF writers that we have by their source; e.g. truth
    vcf_writers: BTreeMap<VariantSource, vcf::io::Writer<Box<dyn std::io::Write>>>
}

impl VariantCategorizer {
    /// Contructor that loads the truth and query VCFs and create the files for writing.
    /// # Arguments
    /// * `truth_vcf_fn` - the input Truth VCF
    /// * `truth_sample_name` - the sample name for the truth VCFs
    /// * `query_vcf_fn` - the input Query VCF
    /// * `query_sample_name` - the sample name for the query VCFs
    /// * `debug_folder` - root of the debug folder, this will create files inside it
    /// # Errors
    /// * if there are any problems opening the input files or writing files to the debug folder
    pub fn new(
        truth_vcf_fn: &Path,
        truth_sample_name: String,
        query_vcf_fn: &Path,
        query_sample_name: String,
        debug_folder: &Path
    ) -> anyhow::Result<Self> {
        // Open the VCF files
        let mut truth_vcf_reader = VcfBuilder::default()
            .build_from_path(truth_vcf_fn)
            .with_context(|| format!("Error while opening {truth_vcf_fn:?}:"))?;
        let mut query_vcf_reader = VcfBuilder::default()
            .build_from_path(query_vcf_fn)
            .with_context(|| format!("Error while opening {query_vcf_fn:?}:"))?;

        // get the headers also
        let mut truth_vcf_header = truth_vcf_reader.read_header()
            .with_context(|| format!("Error while reading header of {truth_vcf_fn:?}:"))?;
        let mut query_vcf_header = query_vcf_reader.read_header()
            .with_context(|| format!("Error while reading header of {query_vcf_fn:?}:"))?;

        let ver: &str = crate::cli::core::FULL_VERSION.as_str(); // clippy gets weird about direct access
        let cli_version = format!("\"{}\"", ver);
        let cli_string = format!("\"{}\"", std::env::args().collect::<Vec<String>>().join(" "));
        truth_vcf_header.insert("aardvark_version".parse()?, vcf::header::record::Value::from(cli_version.clone()))?;
        truth_vcf_header.insert("aardvark_command".parse()?, vcf::header::record::Value::from(cli_string.clone()))?;
        query_vcf_header.insert("aardvark_version".parse()?, vcf::header::record::Value::from(cli_version))?;
        query_vcf_header.insert("aardvark_command".parse()?, vcf::header::record::Value::from(cli_string))?;

        // add expected and observed allele counts to our header
        let extra_header = [
            (
                BENCHMARK_DECISION_KEY.to_string(),
                Map::<map::Format>::new(map::format::Number::Count(1), map::format::Type::String, "Benchmark Decision for call (TP/FP/FN)"),
            ),
            (
                EXPECTED_ALLELE_KEY.to_string(),
                Map::<map::Format>::new(map::format::Number::Count(1), map::format::Type::Integer, "Expected Allele count for this genotype"),
            ),
            (
                OBSERVED_ALLELE_KEY.to_string(),
                Map::<map::Format>::new(map::format::Number::Count(1), map::format::Type::Integer, "Observed Allele count for this genotype")
            ),
            (
                REGION_ID.to_string(),
                Map::<map::Format>::new(map::format::Number::Count(1), map::format::Type::Integer, "Region ID for the comparison")
            )
        ];
        
        for (header_key, header_value) in extra_header.iter() {
            truth_vcf_header.formats_mut().insert(header_key.clone(), header_value.clone());
            query_vcf_header.formats_mut().insert(header_key.clone(), header_value.clone());
        }

        // set the truth and query to just have a single sample
        let truth_samples = truth_vcf_header.sample_names_mut();
        truth_samples.clear();
        truth_samples.insert(truth_sample_name);

        let query_samples = query_vcf_header.sample_names_mut();
        query_samples.clear();
        query_samples.insert(query_sample_name);

        // we have a lot of sub-files to write
        let filetypes = [
            (VariantSource::Truth, "truth.vcf.gz"),
            (VariantSource::Query, "query.vcf.gz"),
        ];

        let mut vcf_writers: BTreeMap<VariantSource, vcf::io::Writer<Box<dyn std::io::Write>>> = Default::default();
        for (source, filename) in filetypes.into_iter() {
            // make the VCF and write the header
            let vcf_fn = debug_folder.join(filename);
            debug!("Opening {vcf_fn:?} for writing...");
            let mut vcf_writer = vcf::io::writer::Builder::default()
                .set_compression_method(noodles::vcf::io::CompressionMethod::Bgzf)
                .build_from_path(vcf_fn)?;

            // write the header based on the data source
            match source {
                VariantSource::Truth => vcf_writer.write_header(&truth_vcf_header)?,
                VariantSource::Query => vcf_writer.write_header(&query_vcf_header)?,
            };

            // save it to our writers list
            assert!(vcf_writers.insert(source, vcf_writer).is_none());
        }

        Ok(Self {
            truth_header: truth_vcf_header,
            query_header: query_vcf_header,
            vcf_writers
        })
    }

    /// Core result writer, which currently does not check the input order at all; so user beware.
    /// # Arguments
    /// * `region` - the comparison region, which contains the loaded variants
    /// * `benchmark` - the results from our comparison, which contains the classification information
    pub fn write_results(&mut self, region: &CompareRegion, benchmark: &CompareBenchmark) -> anyhow::Result<()> {
        // sanity checks
        if region.truth_variants().len() != benchmark.truth_variant_data().len() {
            bail!("Region and benchmark truth lengths are different for B#{}", region.region_id());
        }

        if region.query_variants().len() != benchmark.query_variant_data().len() {
            bail!("Region and benchmark query lengths are different for B#{}", region.region_id());
        }

        self.write_variants(VariantSource::Truth, region.region_id(), region.coordinates().chrom(), region.truth_variants(), region.truth_zygosity(), benchmark.truth_variant_data())?;
        self.write_variants(VariantSource::Query, region.region_id(), region.coordinates().chrom(), region.query_variants(), region.query_zygosity(), benchmark.query_variant_data())?;
        Ok(())
    }

    /// Will write out a set of variants to the correct files
    /// # Arguments
    /// * `source` - selects either Truth or Query types
    /// * `region_id` - passed through to FORMAT:RI
    /// * `chrom` - chromosome
    /// * `variants` - the variants that were parsed from the input, not necessarily identical to the input VCF
    /// * `zygosity` - the parsed zygosity
    /// * `bench_data` - the comparison data, which determines whether the variant is a TP, FN, or FP
    fn write_variants(&mut self, source: VariantSource, region_id: u64, chrom: &str, variants: &[Variant], zygosity: &[PhasedZygosity], bench_data: &[VariantMetrics]) -> anyhow::Result<()> {
        let vcf_header = match source {
            VariantSource::Truth => &self.truth_header,
            VariantSource::Query => &self.query_header,
        };
        for ((variant, &zygosity), bd) in variants.iter().zip(zygosity.iter()).zip(bench_data.iter()) {
            // convert REF and ALT back into strings
            let ref_allele = std::str::from_utf8(variant.allele0())?;
            let alternate_bases = record_buf::AlternateBases::from(
                vec![String::from_utf8(variant.allele1().to_vec())?]
            );

            let genotypes = match zygosity {
                PhasedZygosity::Unknown => ".",
                PhasedZygosity::HomozygousReference => "0/0",
                PhasedZygosity::UnphasedHeterozygous => "0/1",
                PhasedZygosity::PhasedHet01 => "0|1",
                PhasedZygosity::PhasedHet10 => "1|0",
                PhasedZygosity::HomozygousAlternate => "1/1",
            };
            let classification = bd.classification();

            // anything going into the final record has a key-value pair for the FORMAT tag
            let format_keys: record_buf::samples::Keys = [
                vcf_key::GENOTYPE.to_string(),
                BENCHMARK_DECISION_KEY.to_string(),
                EXPECTED_ALLELE_KEY.to_string(),
                OBSERVED_ALLELE_KEY.to_string(),
                REGION_ID.to_string()
            ].into_iter().collect();
            let values = vec![
                // nested vec for multi-sample
                vec![
                    Some(record_buf::samples::sample::Value::from(genotypes)),
                    Some(record_buf::samples::sample::Value::from(classification.as_ref())),
                    Some(record_buf::samples::sample::Value::from(bd.expected_allele_count() as i32)),
                    Some(record_buf::samples::sample::Value::from(bd.observed_allele_count() as i32)),
                    Some(record_buf::samples::sample::Value::from(region_id as i32)) // TODO: this might be in danger of overflowing
                ]
            ];

            // save the key values as Samples
            let samples = record_buf::Samples::new(
                format_keys,
                values
            );

            // now fill out the record
            let record = noodles::vcf::variant::RecordBuf::builder()
                .set_reference_sequence_name(chrom)
                .set_variant_start(Position::new(variant.position() as usize + 1).unwrap()) // Position is 1-based in noodles world
                .set_reference_bases(ref_allele)
                .set_alternate_bases(alternate_bases)
                .set_samples(samples)
                .build();

            // save the new records to the correct VCF file
            let vcf_writer = self.vcf_writers.get_mut(&source).unwrap();
            vcf_writer.write_variant_record(vcf_header, &record)?;
        }

        Ok(())
    }
}

/// Helpful utility function to index all outputs from a VariantCategorizer.
/// Critically, it consumes the writer to finalize all outputs prior to indexing.
/// # Arguments
/// * `debug_folder` - the output folder we are indexing
/// * `categorizer` - the writer that gets consumed
/// # Errors
/// * if the noodles indexing throws any errors; this could happen if the BED file provided is not sorted
pub fn index_categorizer(debug_folder: &Path, categorizer: VariantCategorizer) -> anyhow::Result<()> {
    // drop it from memory, this forces all the finalizing to happen
    std::mem::drop(categorizer);

    // now index everything
    let extensions = [
        "truth.vcf.gz",
        "query.vcf.gz",
    ];

    for extension in extensions.iter() {
        let vcf_fn = debug_folder.join(extension);
        debug!("Generating index for {vcf_fn:?}...");
        crate::writers::noodles_idx::index_vcf(&vcf_fn)
            .with_context(|| format!("Error while writing index for {vcf_fn:?}"))?;
    }
    Ok(())
}