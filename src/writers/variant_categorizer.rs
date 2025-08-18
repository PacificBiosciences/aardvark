
use anyhow::{Context, bail};
use log::debug;
use noodles::bgzf;
use noodles::core::Position;
use noodles::vcf;
use noodles::vcf::header::record::value::{Map, map};
use noodles::vcf::variant::io::Write;
use noodles_util::variant::io::indexed_reader::Builder as VcfBuilder;
use noodles::vcf::variant::record::samples::keys::key as vcf_key;
use noodles::vcf::variant::record_buf;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::{Path, PathBuf};

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
    /// Either truth or query
    source: VariantSource,
    /// Track the filename
    filename: PathBuf,
    /// Header that goes into the VCF
    vcf_header: vcf::Header,
    /// The actual writer
    vcf_writer: vcf::io::Writer<bgzf::io::MultithreadedWriter<File>>
}

impl VariantCategorizer {
    /// Contructor that loads the truth and query VCFs and create a pair of writers for the outputs.
    /// # Arguments
    /// * `truth_vcf_fn` - the input Truth VCF
    /// * `truth_sample_name` - the sample name for the truth VCFs
    /// * `query_vcf_fn` - the input Query VCF
    /// * `query_sample_name` - the sample name for the query VCFs
    /// * `debug_folder` - root of the debug folder, this will create files inside it
    /// * `threads` - writer threads per gzip writer
    /// # Errors
    /// * if there are any problems opening the input files or writing files to the debug folder
    pub fn new_writer_pair(
        truth_vcf_fn: &Path,
        truth_sample_name: String,
        query_vcf_fn: &Path,
        query_sample_name: String,
        debug_folder: &Path,
        threads: usize
    ) -> anyhow::Result<(Self, Self)> {
        let mut source_lookup: BTreeMap<VariantSource, Self> = Default::default();

        // constants we add
        let ver: &str = crate::cli::core::FULL_VERSION.as_str(); // clippy gets weird about direct access
        let cli_version = format!("\"{ver}\"");
        let cli_string = format!("\"{}\"", std::env::args().collect::<Vec<String>>().join(" "));
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

        for variant_source in [VariantSource::Truth, VariantSource::Query] {
            let input_fn = match variant_source {
                VariantSource::Truth => truth_vcf_fn,
                VariantSource::Query => query_vcf_fn
            };
            let sample_name = match variant_source {
                VariantSource::Truth => truth_sample_name.clone(),
                VariantSource::Query => query_sample_name.clone()
            };

            // Open the VCF file
            let mut vcf_reader = VcfBuilder::default()
                .build_from_path(input_fn)
                .with_context(|| format!("Error while opening {input_fn:?}:"))?;

            // get the header
            let mut vcf_header = vcf_reader.read_header()
                .with_context(|| format!("Error while reading header of {input_fn:?}:"))?;

            // update the header
            vcf_header.insert("aardvark_version".parse()?, vcf::header::record::Value::from(cli_version.clone()))?;
            vcf_header.insert("aardvark_command".parse()?, vcf::header::record::Value::from(cli_string.clone()))?;
            for (header_key, header_value) in extra_header.iter() {
                vcf_header.formats_mut().insert(header_key.clone(), header_value.clone());
            }

            let samples = vcf_header.sample_names_mut();
            samples.clear();
            samples.insert(sample_name);

            let out_extension = match variant_source {
                VariantSource::Truth => "truth.vcf.gz",
                VariantSource::Query => "query.vcf.gz"
            };

            // make the VCF and write the header
            let vcf_fn = debug_folder.join(out_extension);
            debug!("Opening {vcf_fn:?} for writing...");

            let file = File::create(&vcf_fn)?;
            let w_threads = std::num::NonZeroUsize::new(threads.clamp(1, 4)).unwrap();
            let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(w_threads, file);
            let mut vcf_writer = vcf::io::Writer::new(bgzf_writer);

            // write the header based on the data source
            vcf_writer.write_header(&vcf_header)?;

            // save it to our writers list
            let vc = Self {
                source: variant_source,
                filename: vcf_fn,
                vcf_header,
                vcf_writer
            };
            assert!(source_lookup.insert(variant_source, vc).is_none());
        }

        // now pop and return
        let truth_writer = source_lookup.remove(&VariantSource::Truth).unwrap();
        let query_writer = source_lookup.remove(&VariantSource::Query).unwrap();
        Ok((truth_writer, query_writer))
    }

    /// Core result writer, which currently does not check the input order at all; so user beware.
    /// # Arguments
    /// * `region` - the comparison region, which contains the loaded variants
    /// * `benchmark` - the results from our comparison, which contains the classification information
    pub fn write_results(&mut self, region: &CompareRegion, benchmark: &CompareBenchmark) -> anyhow::Result<()> {
        match self.source {
            VariantSource::Truth => {
                // sanity check then write it
                if region.truth_variants().len() != benchmark.truth_variant_data().len() {
                    bail!("Region and benchmark truth lengths are different for B#{}", region.region_id());
                }
                self.write_variants(VariantSource::Truth, region.region_id(), region.coordinates().chrom(), region.truth_variants(), region.truth_zygosity(), benchmark.truth_variant_data())?;
            },
            VariantSource::Query => {
                // sanity check then write it
                if region.query_variants().len() != benchmark.query_variant_data().len() {
                    bail!("Region and benchmark query lengths are different for B#{}", region.region_id());
                }
                self.write_variants(VariantSource::Query, region.region_id(), region.coordinates().chrom(), region.query_variants(), region.query_zygosity(), benchmark.query_variant_data())?;
            },
        }
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
        // make sure we didn't do something stupid
        assert_eq!(self.source, source);
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
            self.vcf_writer.write_variant_record(&self.vcf_header, &record)?;
        }

        Ok(())
    }

    pub fn filename(&self) -> &Path {
        &self.filename
    }
}

/// Helpful utility function to index all outputs from a VariantCategorizer.
/// Critically, it consumes the writer to finalize all outputs prior to indexing.
/// # Arguments
/// * `categorizer` - the writer that gets consumed
/// # Errors
/// * if the noodles indexing throws any errors; this could happen if the BED file provided is not sorted
pub fn index_categorizer(categorizer: VariantCategorizer) -> anyhow::Result<()> {
    // copy all the filenames, so we can index
    let vcf_fn = categorizer.filename().to_path_buf();

    // drop it from memory, this forces all the finalizing to happen
    std::mem::drop(categorizer);

    // now index everything
    debug!("Generating index for {vcf_fn:?}...");
    crate::writers::noodles_idx::index_vcf(&vcf_fn)
        .with_context(|| format!("Error while writing index for {vcf_fn:?}"))?;
    Ok(())
}