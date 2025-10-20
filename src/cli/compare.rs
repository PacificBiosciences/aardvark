
use anyhow::bail;
use clap::Args;
use log::info;
use serde::Serialize;
use std::path::PathBuf;

use crate::cli::core::{check_optional_filename, check_required_filename, AFTER_HELP, FULL_VERSION};
use crate::parsing::noodles_helper::get_vcf_sample_name;

#[derive(Args, Clone, Default, Serialize)]
#[clap(author, about, 
    after_help = &**AFTER_HELP
)]
pub struct CompareSettings {
    #[clap(default_value = "")]
    #[clap(hide = true)]
    aardvark_version: String,

    /// Reference FASTA file
    #[clap(required = true)]
    #[clap(short = 'r')]
    #[clap(long = "reference")]
    #[clap(value_name = "FASTA")]
    #[clap(help_heading = Some("Input/Output"))]
    pub reference_fn: PathBuf,

    /// Truth variant call file (VCF)
    #[clap(required = true)]
    #[clap(short = 't')]
    #[clap(long = "truth-vcf")]
    #[clap(value_name = "VCF")]
    #[clap(help_heading = Some("Input/Output"))]
    pub truth_vcf_filename: PathBuf,

    /// Query variant call file (VCF)
    #[clap(required = true)]
    #[clap(short = 'q')]
    #[clap(long = "query-vcf")]
    #[clap(value_name = "VCF")]
    #[clap(help_heading = Some("Input/Output"))]
    pub query_vcf_filename: PathBuf,

    /// Confidence regions (BED)
    #[clap(short = 'b')]
    #[clap(long = "regions")]
    #[clap(value_name = "BED")]
    #[clap(help_heading = Some("Input/Output"))]
    pub regions: Option<PathBuf>,

    /// Stratifications, specifically the root file-of-filenames TSV
    #[clap(short = 's')]
    #[clap(long = "stratification")]
    #[clap(value_name = "TSV")]
    #[clap(help_heading = Some("Input/Output"))]
    pub stratifications: Option<PathBuf>,

    /// Output directory containing summary and VCFs
    #[clap(required = true)]
    #[clap(short = 'o')]
    #[clap(long = "output-dir")]
    #[clap(value_name = "DIR")]
    #[clap(help_heading = Some("Input/Output"))]
    pub output_folder: PathBuf,

    /// Optional output debug folder
    #[clap(long = "output-debug")]
    #[clap(value_name = "DIR")]
    #[clap(help_heading = Some("Input/Output"))]
    pub debug_folder: Option<PathBuf>,

    /// Optional comparison label for the summary output
    #[clap(long = "compare-label")]
    #[clap(value_name = "LABEL")]
    #[clap(help_heading = Some("Input/Output"))]
    #[clap(default_value = "compare")]
    pub compare_label: String,

    /// The sample name to use in the truth VCF [default: first sample]
    #[clap(long = "truth-sample")]
    #[clap(value_name = "SAMPLE")]
    #[clap(help_heading = Some("Input/Output"))]
    #[clap(default_value = "", hide_default_value = true)]
    pub truth_sample: String,

    /// The sample name to use in the query VCF [default: first sample]
    #[clap(long = "query-sample")]
    #[clap(value_name = "SAMPLE")]
    #[clap(help_heading = Some("Input/Output"))]
    #[clap(default_value = "", hide_default_value = true)]
    pub query_sample: String,

    /// The minimum gap (bp) between variants to split into separate sub-regions
    #[clap(long = "min-variant-gap")]
    #[clap(value_name = "BP")]
    #[clap(help_heading = Some("Region generation"))]
    #[clap(default_value = "50")]
    pub min_variant_gap: usize,

    /// Disables variant trimming, which may have a negative impact on accuracy
    #[clap(long = "disable-variant-trimming")]
    #[clap(help_heading = Some("Region generation"))]
    pub disable_variant_trimming: bool,

    /// Maximum edit distance in the WFA comparison before quitting
    #[clap(long = "max-edit-distance")]
    #[clap(value_name = "INT")]
    #[clap(help_heading = Some("Compare parameters"))]
    #[clap(default_value = "5000")]
    #[clap(hide = true)] // if you remove this, make sure you re-enable the CLI outputs
    pub max_edit_distance: usize,

    /// Enables an exact-match compute shortcut at the cost of variant-level assessment accuracy
    #[clap(long = "enable-exact-shortcut")]
    #[clap(help_heading = Some("Compare parameters"))]
    #[clap(hide = true)] // if you remove this, make sure you re-enable the CLI outputs
    pub enable_exact_shortcut: bool,

    /// Enables the haplotype scoring metrics
    #[clap(long = "enable-haplotype-metrics")]
    #[clap(help_heading = Some("Optional metrics"))]
    #[clap(default_value = "false")]
    pub enable_haplotype_scoring: bool,

    /// Enables the weighted haplotype scoring metrics
    #[clap(long = "enable-weighted-haplotype-metrics")]
    #[clap(help_heading = Some("Optional metrics"))]
    #[clap(default_value = "false")]
    pub enable_weighted_haplotype_scoring: bool,

    /// Enables the record-basepair scoring metrics
    #[clap(long = "enable-record-basepair-metrics")]
    #[clap(help_heading = Some("Optional metrics"))]
    #[clap(default_value = "false")]
    pub enable_record_basepair_scoring: bool,

    /// Number of threads to use in the benchmarking step
    #[clap(long = "threads")]
    #[clap(value_name = "THREADS")]
    #[clap(default_value = "1")]
    pub threads: usize,

    /// Enable verbose output.
    #[clap(short = 'v')]
    #[clap(long = "verbose")]
    #[clap(action = clap::ArgAction::Count)]
    pub verbosity: u8,

    // Debug options that are generally hidden and just for quick testing
    /// Skips a number of sub-problems (debug only); non-0 values may create partial values in output
    #[clap(hide = true)]
    #[clap(long = "skip")]
    #[clap(default_value = "0")]
    pub skip_blocks: usize,

    /// Takes a number of sub-problems (debug only); non-0 values may create partial values in output
    #[clap(hide = true)]
    #[clap(long = "take")]
    #[clap(default_value = "0")]
    pub take_blocks: usize,
}

pub fn check_compare_settings(mut settings: CompareSettings) -> anyhow::Result<CompareSettings> {
    // hard code the version in
    settings.aardvark_version = FULL_VERSION.clone();
    info!("Aardvark version: {:?}", &settings.aardvark_version);
    info!("Sub-command: compare");
    info!("Inputs:");

    // check for all the required input files
    check_required_filename(&settings.reference_fn, "Reference FASTA")?;
    check_required_filename(&settings.truth_vcf_filename, "Truth VCF")?;
    check_required_filename(&settings.query_vcf_filename, "Query VCF")?;
    check_optional_filename(settings.regions.as_deref(), "Regions")?;
    check_optional_filename(settings.stratifications.as_deref(), "Stratifications")?;
    
    // dump stuff to the logger
    info!("\tReference: {:?}", &settings.reference_fn);
    info!("\tTruth VCF: {:?}", &settings.truth_vcf_filename);
    if settings.truth_sample.is_empty() {
        settings.truth_sample = get_vcf_sample_name(&settings.truth_vcf_filename, 0)?;
    }
    info!("\tTruth sample: {:?}", &settings.truth_sample);
    info!("\tQuery VCF: {:?}", &settings.query_vcf_filename);
    if settings.query_sample.is_empty() {
        settings.query_sample = get_vcf_sample_name(&settings.query_vcf_filename, 0)?;
    }
    info!("\tQuery sample: {:?}", &settings.query_sample);
    if let Some(hcr_fn) = settings.regions.as_deref() {
        info!("\tRegions: {hcr_fn:?}");
    } else {
        info!("\tRegions: None");
    }
    if let Some(filename) = settings.stratifications.as_deref() {
        info!("\tStratifications: {filename:?}");
    } else {
        info!("\tStratifications: None");
    }

    // check for the haplotype scoring metrics
    if settings.enable_haplotype_scoring || settings.enable_weighted_haplotype_scoring || settings.enable_record_basepair_scoring {
        info!("Optional metrics:");
        info!("\tHaplotype scoring: {}", if settings.enable_haplotype_scoring { "ENABLED" } else { "DISABLED" });
        info!("\tWeighted haplotype scoring: {}", if settings.enable_weighted_haplotype_scoring { "ENABLED" } else { "DISABLED" });
        info!("\tRecord-basepair scoring: {}", if settings.enable_record_basepair_scoring { "ENABLED" } else { "DISABLED" });
    }

    // 0 is just a sentinel for everything
    if settings.take_blocks == 0 {
        settings.take_blocks = usize::MAX;
    }
    
    // outputs
    info!("Outputs:");
    info!("\tCompare label: {:?}", &settings.compare_label);
    info!("\tOutput folder: {:?}", &settings.output_folder);
    if let Some(debug_folder) = settings.debug_folder.as_ref() {
        info!("\tDebug folder: {debug_folder:?}");
    }

    // other misc parameters
    info!("Region generation parameters:");
    if settings.min_variant_gap == 0 {
        bail!("--min-variant-gap must be >0");
    }
    info!("\tMinimum variant gap: {}", settings.min_variant_gap);
    info!("\tVariant trimming: {}", if settings.disable_variant_trimming { "DISABLED "} else { "ENABLED" });

    if settings.enable_exact_shortcut {
        info!("Compare parameters:");
        // info!("\tMax edit distance: {}", settings.max_edit_distance); // we removed this
        info!("\tExact match shortcut: {}", if settings.enable_exact_shortcut { "ENABLED" } else { "DISABLED" });
    }

    if settings.threads == 0 {
        settings.threads = 1;
    }
    info!("Processing threads: {}", settings.threads);

    Ok(settings)
}
