
use anyhow::ensure;
use clap::Args;
use log::info;
use serde::Serialize;
use std::path::PathBuf;
use strum_macros::EnumString;

use crate::cli::core::{check_optional_filename, check_required_filename, AFTER_HELP, FULL_VERSION};
use crate::parsing::noodles_helper::get_vcf_sample_name;

#[derive(Clone, Copy, Default, Debug, strum_macros::Display, EnumString, Serialize, clap::ValueEnum)]
pub enum MergeStrategy {
    /// Must exactly match within the window for all samples
    #[default]
    #[strum(ascii_case_insensitive, serialize = "exact")]
    #[clap(name = "exact")]
    Exact,
    /// Must either exactly match or be missing completely from all samples; i.e., no conflicting variant calls
    #[strum(ascii_case_insensitive, serialize = "no_conflict")]
    #[clap(name = "no_conflict")]
    NoConflict,
    /// Allows for a majority to win even if others disagree
    #[strum(ascii_case_insensitive, serialize = "majority")]
    #[clap(name = "majority")]
    MajorityVote,
    /// Enables all optional merge strategies
    #[strum(ascii_case_insensitive, serialize = "all")]
    #[clap(name = "all")]
    AllOptions,
}

#[derive(Args, Clone, Default, Serialize)]
#[clap(author, about, 
    after_help = &**AFTER_HELP
)]
pub struct MergeSettings {
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

    /// Input variant call file (VCF), provided in priority order
    #[clap(required = true)]
    #[clap(short = 'i')]
    #[clap(long = "input-vcf")]
    #[clap(value_name = "VCF")]
    #[clap(help_heading = Some("Input/Output"))]
    pub vcf_filenames: Vec<PathBuf>,

    /// The sample name to use in the corresponding VCF [default: first sample]
    #[clap(short = 's')]
    #[clap(long = "vcf-sample")]
    #[clap(value_name = "SAMPLE")]
    #[clap(help_heading = Some("Input/Output"))]
    pub vcf_samples: Vec<String>,

    /// The annotation tag to use for the corresponding VCF [default: "vcf_#"]
    #[clap(short = 't')]
    #[clap(long = "vcf-tag")]
    #[clap(value_name = "SAMPLE")]
    #[clap(help_heading = Some("Input/Output"))]
    pub vcf_tags: Vec<String>,

    /// Regions to perform the merge (BED)
    #[clap(short = 'b')]
    #[clap(long = "regions")]
    #[clap(value_name = "BED")]
    #[clap(help_heading = Some("Input/Output"))]
    pub merge_regions: Option<PathBuf>,

    /// Output VCF folder
    #[clap(short = 'o')]
    #[clap(long = "output-vcfs")]
    #[clap(value_name = "DIR")]
    #[clap(help_heading = Some("Input/Output"))]
    pub output_vcf_folder: PathBuf,

    /// Output summary file (CSV/TSV)
    #[clap(long = "output-summary")]
    #[clap(value_name = "TSV")]
    #[clap(help_heading = Some("Input/Output"))]
    pub output_summary_filename: Option<PathBuf>,

    /// Optional output debug folder
    #[clap(long = "output-debug")]
    #[clap(value_name = "DIR")]
    #[clap(help_heading = Some("Input/Output"))]
    pub debug_folder: Option<PathBuf>,

    /*
    // TODO: might add this in later, but it's not needed currently
    /// Maximum edit distance in the WFA comparison before quitting
    #[clap(long = "max-edit-distance")]
    #[clap(value_name = "INT")]
    #[clap(help_heading = Some("Compare parameters"))]
    #[clap(default_value = "5000")]
    pub max_edit_distance: usize,
    */

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

    /// Selects pre-set merge strategy for inclusion of a variant
    #[clap(long = "merge-strategy")]
    #[clap(value_name = "STRAT")]
    #[clap(help_heading = Some("Merge parameters"))]
    pub merge_strategy: Option<MergeStrategy>,

    /// Enables merging if no conflicts are detected
    #[clap(long = "enable-no-conflict")]
    #[clap(help_heading = Some("Merge parameters"))]
    pub enable_no_conflict: bool,

    /// Enables merging if the majority of inputs agree
    #[clap(long = "enable-voting")]
    #[clap(help_heading = Some("Merge parameters"))]
    pub enable_voting: bool,

    /// Sets a VCF index to select to always get selected in the event of conflict
    #[clap(long = "conflict-select")]
    #[clap(value_name = "INDEX")]
    #[clap(help_heading = Some("Merge parameters"))]
    pub conflict_selection: Option<usize>,

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

pub fn check_merge_settings(mut settings: MergeSettings) -> anyhow::Result<MergeSettings> {
    // hard code the version in
    settings.aardvark_version = FULL_VERSION.clone();
    info!("Aardvark version: {:?}", &settings.aardvark_version);
    info!("Sub-command: merge");
    info!("Inputs:");

    // check for all the required input files
    check_required_filename(&settings.reference_fn, "Reference FASTA")?;
    info!("\tReference: {:?}", &settings.reference_fn);
    check_optional_filename(settings.merge_regions.as_deref(), "Merge regions")?;
    if let Some(hcr_fn) = settings.merge_regions.as_deref() {
        info!("\tMerge regions: {hcr_fn:?}");
    } else {
        info!("\tMerge regions: None");
    }
    
    // check the input VCFs and corresponding metadata
    for (i, i_vcf) in settings.vcf_filenames.iter().enumerate() {
        check_required_filename(i_vcf, format!("Input VCF #{i}").as_str())?;
        info!("\tInput VCF #{i}: {i_vcf:?}");

        if settings.vcf_samples.len() <= i {
            settings.vcf_samples.push(get_vcf_sample_name(i_vcf, 0)?);
        }
        info!("\t\tSample name: {:?}", settings.vcf_samples[i]);

        if settings.vcf_tags.len() <= i {
            settings.vcf_tags.push(format!("vcf_{i}"));
        }
        info!("\t\tOutput tag: {:?}", settings.vcf_tags[i]);
    }
    
    // outputs
    info!("Outputs:");
    info!("\tVCF folder: {:?}", &settings.output_vcf_folder);
    info!("\tSummary: {:?}", &settings.output_summary_filename);
    if let Some(debug_folder) = settings.debug_folder.as_ref() {
        info!("\tDebug folder: {debug_folder:?}");
    }

    // other misc parameters
    info!("Region generation parameters:");
    ensure!(settings.min_variant_gap > 0, "--min-variant-gap must be >0");
    info!("\tMinimum variant gap: {}", settings.min_variant_gap);
    info!("\tVariant trimming: {}", if settings.disable_variant_trimming { "DISABLED "} else { "ENABLED" });

    // info!("Compare parameters:");
    // info!("\tMax edit distance: {}", settings.max_edit_distance);

    info!("Merge parameters:");
    if let Some(merge_strat) = settings.merge_strategy {
        info!("\tPre-set merge strategy: {merge_strat}");
        match merge_strat {
            MergeStrategy::Exact => {}, // default
            MergeStrategy::NoConflict => {
                settings.enable_no_conflict = true;
            },
            MergeStrategy::MajorityVote => {
                settings.enable_voting = true;
            },
            MergeStrategy::AllOptions => {
                settings.enable_no_conflict = true;
                settings.enable_voting = true;
            },
        };
    } else {
        info!("\tPre-set merge strategy: None")
    }
    info!("\tNo conflict blocks: {}", if settings.enable_no_conflict { "ENABLED" } else { "DISABLED" });
    info!("\tMajority voting blocks: {}", if settings.enable_voting { "ENABLED" } else { "DISABLED" });

    if let Some(v_index) = settings.conflict_selection {
        ensure!(v_index < settings.vcf_filenames.len(), "--conflict-selection index is greater than number of provided VCFs");
        info!("\tConflict selection: input #{} -> \"{}\"", v_index, settings.vcf_tags[v_index]);
    } else {
        info!("\tConflict selection: None");
    }

    // 0 is just a sentinel for everything
    if settings.take_blocks == 0 {
        settings.take_blocks = usize::MAX;
    }
    if settings.threads == 0 {
        settings.threads = 1;
    }
    info!("Processing threads: {}", settings.threads);

    Ok(settings)
}
