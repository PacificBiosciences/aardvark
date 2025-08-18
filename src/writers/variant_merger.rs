
use anyhow::{Context, bail, ensure};
use log::{debug, info, error};
use noodles::bgzf;
use noodles::core::Position;
use noodles::vcf;
use noodles::vcf::header::record::value::{Map, map};
use noodles::vcf::variant::io::Write;
use noodles_util::variant::io::indexed_reader::Builder as VcfBuilder;
use noodles::vcf::variant::record::samples::keys::key as vcf_key;
use noodles::vcf::variant::record_buf;
use std::fs::File;
use std::path::Path;
use rustc_hash::FxHashMap as HashMap;

use crate::data_types::merge_benchmark::{MergeBenchmark, MergeClassification};
use crate::data_types::multi_region::MultiRegion;
use crate::data_types::phase_enums::PhasedZygosity;
use crate::writers::variant_categorizer::REGION_ID;

/// INFO/SOURCES key label
pub const INFO_KEY_SOURCES: &str = "SOURCES";
/// INFO/MERGE key label
pub const INFO_KEY_MERGE_REASON: &str = "MR";

/// Wrapper struct that will parse the results from the merge processes and determine where to write variants
pub struct VariantMerger {
    /// Header that goes into the merged VCF, copied from the primary input
    merged_header: vcf::Header,
    /// Writer for the primary merged VCF
    merged_writer: vcf::io::Writer<bgzf::io::MultithreadedWriter<File>>,
    /// Pre-computed labels of classification to a joined string
    prelabels: HashMap<MergeClassification, record_buf::info::field::Value>,
    /// Set of all labels by index
    input_labels: Vec<String>,
    /// Writes out the relevant bed regions
    region_writer: noodles::bed::io::Writer<4, std::io::BufWriter<bgzf::io::MultithreadedWriter<std::fs::File>>>,
    /// Writes out the relevant bed regions
    failed_region_writer: noodles::bed::io::Writer<4, std::io::BufWriter<bgzf::io::MultithreadedWriter<std::fs::File>>>
}

// TODO: technically, this is writing 3 big files (well, 2 big, 1 small).
//      If this ever becomes a bottleneck, we can give it the compare_parallel.rs treatment and split the files out.
//      I don't think it's worth the engineering effort ATM.

impl VariantMerger {
    /// Constructor
    /// # Arguments
    /// * `input_vcfs` - used to copy headers
    /// * `input_labels` - used to tag passing variants
    /// * `out_vcf_folder` - defines where all output VCFs go
    /// * `out_sample_name` - sample name for the output
    /// * `threads` - number of threads to use per writer
    pub fn new<P: AsRef<Path> + std::fmt::Debug>(
        input_vcfs: &[P],
        input_labels: Vec<String>,
        out_vcf_folder: &Path,
        out_sample_name: String,
        threads: usize
    ) -> anyhow::Result<Self> {
        // sanity checks
        ensure!(input_vcfs.len() == input_labels.len(), "Expected 1-to-1 relationship between input_vcfs and input_labels");

        // generate the main file first
        info!("Creating output VCF folder at {out_vcf_folder:?}...");
        match std::fs::create_dir_all(out_vcf_folder) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while creating output VCF folder: {e}");
                std::process::exit(exitcode::IOERR);
            }
        }

        // Open the VCF files
        let primary_vcf_fn = &input_vcfs[0];
        let mut vcf_reader = VcfBuilder::default()
            .build_from_path(primary_vcf_fn)
            .with_context(|| format!("Error while opening {primary_vcf_fn:?}:"))?;
        
        // get the headers also
        let mut vcf_header = vcf_reader.read_header()
            .with_context(|| format!("Error while reading header of {primary_vcf_fn:?}:"))?;

        let ver: &str = crate::cli::core::FULL_VERSION.as_str(); // clippy gets weird about direct access
        let cli_version = format!("\"{ver}\"");
        let cli_string = format!("\"{}\"", std::env::args().collect::<Vec<String>>().join(" "));
        vcf_header.insert("aardvark_version".parse()?, vcf::header::record::Value::from(cli_version))?;
        vcf_header.insert("aardvark_command".parse()?, vcf::header::record::Value::from(cli_string))?;

        // add any extra INFO fields
        let info_header = [
            (
                INFO_KEY_SOURCES.to_string(),
                Map::<map::Info>::new(map::info::Number::Unknown, map::info::Type::String, "List of tools or technologies that called the same record")
            ),
            (
                INFO_KEY_MERGE_REASON.to_string(),
                Map::<map::Info>::new(map::info::Number::Count(1), map::info::Type::String, "The reason this record was allowed in the merge")
            )
        ];
        for (header_key, header_value) in info_header.iter() {
            vcf_header.infos_mut().insert(header_key.clone(), header_value.clone());
        }

        // add any extra FORMAT
        let extra_header = [
            (
                REGION_ID.to_string(),
                Map::<map::Format>::new(map::format::Number::Count(1), map::format::Type::Integer, "Region ID for the comparison")
            )
        ];
        for (header_key, header_value) in extra_header.iter() {
            vcf_header.formats_mut().insert(header_key.clone(), header_value.clone());
        }

        // set the sample name for the merge
        let sample_names = vcf_header.sample_names_mut();
        sample_names.clear();
        sample_names.insert(out_sample_name);

        // now open the file up for writing
        let out_vcf_fn = out_vcf_folder.join("passing.vcf.gz");
        debug!("Opening {out_vcf_fn:?} for writing...");
        let file = File::create(&out_vcf_fn)?;
        let w_threads = std::num::NonZeroUsize::new(threads.clamp(1, 4)).unwrap();
        let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(w_threads, file);
        let mut vcf_writer = vcf::io::Writer::new(bgzf_writer);
        vcf_writer.write_header(&vcf_header)?;

        // create any common labels
        let prelabels = [
            (MergeClassification::BasepairIdentical, &input_labels)
        ].into_iter().map(|(c, l)| {
            let v = record_buf::info::field::Value::from(
                l.iter().map(|k| Some(k.clone())).collect::<Vec<Option<String>>>()
            );
            (c, v)
        }).collect();

        // lastly, we need a region writer
        let bed_fn = out_vcf_folder.join("regions.bed.gz");
        let file = File::create(&bed_fn)?;
        let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(w_threads, file);
        #[allow(clippy::default_constructed_unit_structs)]
        let region_writer = noodles::bed::io::writer::Builder::<4>::default()
            .build_from_writer(bgzf_writer);

        let failed_bed_fn = out_vcf_folder.join("failed_regions.bed.gz");
        let file = File::create(&failed_bed_fn)?;
        let bgzf_writer = bgzf::io::MultithreadedWriter::with_worker_count(w_threads, file);
        #[allow(clippy::default_constructed_unit_structs)]
        let failed_region_writer = noodles::bed::io::writer::Builder::<4>::default()
            .build_from_writer(bgzf_writer);

        Ok(Self {
            merged_header: vcf_header,
            merged_writer: vcf_writer,
            prelabels,
            input_labels,
            region_writer,
            failed_region_writer
        })
    }

    /// Core result writer, which currently does not check the input order at all; so user beware.
    /// # Arguments
    /// * `region` - the comparison region, which contains the loaded variants
    /// * `benchmark` - the results from our comparison, which contains the classification information
    pub fn write_results(&mut self, region: &MultiRegion, benchmark: &MergeBenchmark) -> anyhow::Result<()> {
        // sanity checks
        ensure!(region.region_id() == benchmark.region_id(), "region and benchmark have different region_ids");

        // figure out where we are copying from
        let opt_source = match benchmark.merge_classification() {
            MergeClassification::Different => None, // no copy
            MergeClassification::NoConflict { indices } |
            MergeClassification::MajorityAgree { indices } => Some(indices[0]), // the lowest index
            MergeClassification::ConflictSelection { index } => Some(*index),
            MergeClassification::BasepairIdentical => Some(0) // all match, so pick 0
        };

        if let Some(source) = opt_source {
            // these are passing regions, save variants and the region
            self.write_variants(source, region, benchmark.merge_classification())?;
            self.write_region(region, benchmark.merge_classification(), true)?;
        } else {
            // else it's going to the failed regions
            self.write_region(region, benchmark.merge_classification(), false)?;
        }
        Ok(())
    }

    /// Will write out a set of variants to the merged output
    /// # Arguments
    /// * `source` - selects the variant origin, usually 0
    /// * `region` - the region we are copying
    /// * `merge_classification` - determine how things are tagged
    fn write_variants(&mut self, source: usize, region: &MultiRegion, merge_classification: &MergeClassification) -> anyhow::Result<()> {
        let vcf_header = &self.merged_header;
        let region_id = region.region_id();
        let chrom = region.coordinates().chrom();
        let variants = &region.variants()[source];
        let zygosity = &region.zygosity()[source];

        if !self.prelabels.contains_key(merge_classification) {
            // we need to construct a new label
            let new_labels: Vec<Option<String>> = match merge_classification {
                MergeClassification::NoConflict { indices } |
                MergeClassification::MajorityAgree { indices } => {
                    indices.iter()
                        .map(|&i| Some(self.input_labels[i].clone()))
                        .collect()
                },
                MergeClassification::ConflictSelection { index } => {
                    vec![Some(self.input_labels[*index].clone())]
                },
                _ => bail!("No pre-label implementation for {merge_classification:?}")
            };

            // put it in a Value and save it
            let v = record_buf::info::field::Value::from(new_labels);
            self.prelabels.insert(merge_classification.clone(), v);
        }
        let source_tag = self.prelabels.get(merge_classification).unwrap();
        let merge_reason = merge_classification.simplify();
        let merge_reason = record_buf::info::field::Value::from(merge_reason);

        for (variant, &zygosity) in variants.iter().zip(zygosity.iter()) {
            // convert REF and ALT back into strings
            let ref_allele = std::str::from_utf8(variant.allele0())?;
            let alternate_bases = record_buf::AlternateBases::from(
                vec![String::from_utf8(variant.allele1().to_vec())?]
            );

            // construct all info fields
            let info_fields: record_buf::Info = [
                (INFO_KEY_SOURCES.to_string(), Some(source_tag.clone())),
                (INFO_KEY_MERGE_REASON.to_string(), Some(merge_reason.clone()))
            ].into_iter().collect();

            let genotypes = match zygosity {
                PhasedZygosity::Unknown => ".",
                PhasedZygosity::HomozygousReference => "0/0",
                PhasedZygosity::UnphasedHeterozygous => "0/1",
                PhasedZygosity::PhasedHet01 => "0|1",
                PhasedZygosity::PhasedHet10 => "1|0",
                PhasedZygosity::HomozygousAlternate => "1/1",
            };

            // anything going into the final record has a key-value pair for the FORMAT tag
            let format_keys: record_buf::samples::Keys = [
                vcf_key::GENOTYPE.to_string(),
                REGION_ID.to_string()
            ].into_iter().collect();
            let values = vec![
                // nested vec for multi-sample
                vec![
                    Some(record_buf::samples::sample::Value::from(genotypes)),
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
                .set_info(info_fields)
                .set_samples(samples)
                .build();

            // save the new records to the merged VCF file
            self.merged_writer.write_variant_record(vcf_header, &record)?;
        }

        Ok(())
    }

    /// This will save a given region to the output region BED file
    /// # Arguments
    /// * `region` - the region that gets saved
    /// * `merge_classification` - primarily used to build the label
    /// * `is_passing` - if true, this goes into the passing regions, otherwise into failed regions
    fn write_region(&mut self, region: &MultiRegion, merge_classification: &MergeClassification, is_passing: bool) -> anyhow::Result<()> {
        // build the record using the provided inputs; this is most positional
        let record = noodles::bed::feature::record_buf::RecordBuf::<4>::builder()
            .set_reference_sequence_name(region.coordinates().chrom())
            .set_feature_start(Position::try_from(region.coordinates().start() as usize)?)
            .set_feature_end(Position::try_from(region.coordinates().end() as usize)?)
            .set_name(format!("{}_{}", merge_classification.simplify(), region.region_id()))
            .build();

        if is_passing {
            self.region_writer.write_feature_record(&record)?;
        } else {
            self.failed_region_writer.write_feature_record(&record)?;
        }
        Ok(())
    }
}

/// Helpful utility function to index all outputs from a VariantMerger.
/// Critically, it consumes the writer to finalize all outputs prior to indexing.
/// # Arguments
/// * `vcf_folder` - the output folder we are indexing
/// * `merger` - the writer that gets consumed
/// # Errors
/// * if the noodles indexing throws any errors; this could happen if the BED file provided is not sorted
pub fn index_merger(vcf_folder: &Path, merger: VariantMerger) -> anyhow::Result<()> {
    // drop it from memory, this forces all the finalizing to happen
    std::mem::drop(merger);

    // now index everything
    let extensions = [
        "passing.vcf.gz",
    ];

    for extension in extensions.iter() {
        let vcf_fn = vcf_folder.join(extension);
        debug!("Generating index for {vcf_fn:?}...");
        crate::writers::noodles_idx::index_vcf(&vcf_fn)
            .with_context(|| format!("Error while writing index for {vcf_fn:?}"))?;
    }

    // now index all bed files
    let extensions = [
        "regions.bed.gz",
        "failed_regions.bed.gz"
    ];

    for extension in extensions.iter() {
        let bed_fn = vcf_folder.join(extension);
        debug!("Generating index for {bed_fn:?}...");
        crate::writers::noodles_idx::index_bed(&bed_fn)
            .with_context(|| format!("Error while writing index for {bed_fn:?}"))?;
    }

    Ok(())
}
