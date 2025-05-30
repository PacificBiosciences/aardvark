
use indicatif::ParallelProgressIterator;
use log::{LevelFilter, debug, error, info, warn};
use rayon::prelude::*;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use std::time::Instant;

use aardvark::cli::compare::{CompareSettings, check_compare_settings};
use aardvark::cli::core::{Commands, get_cli};
use aardvark::cli::merge::{MergeSettings, check_merge_settings};
use aardvark::data_types::compare_benchmark::CompareBenchmark;
use aardvark::data_types::compare_region::CompareRegion;
use aardvark::data_types::merge_benchmark::MergeBenchmark;
use aardvark::data_types::multi_region::MultiRegion;
use aardvark::data_types::summary_metrics::SummaryMetrics;
use aardvark::merge_solver::{MergeConfigBuilder, solve_merge_region};
use aardvark::parsing::region_generation::RegionIterator;
use aardvark::util::json_io::save_json;
use aardvark::util::progress_bar::get_progress_style;
use aardvark::waffle_solver::solve_compare_region;
use aardvark::writers::merge_summary::MergeSummaryWriter;
use aardvark::writers::region_sequence::RegionSequenceWriter;
use aardvark::writers::region_summary::RegionSummaryWriter;
use aardvark::writers::summary::SummaryWriter;
use aardvark::writers::variant_categorizer::{VariantCategorizer, index_categorizer};
use aardvark::writers::variant_merger::{VariantMerger, index_merger};

fn run_compare(settings: CompareSettings) {
    // start the timer
    let start_time = Instant::now();

    // set up logging before we check the other settings
    let filter_level: LevelFilter = match settings.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace
    };
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    let settings = match check_compare_settings(settings) {
        Ok(s) => s,
        Err(e) => {
            error!("Error while verifying settings: {e:#}");
            std::process::exit(exitcode::CONFIG);
        }
    };

    // set up the number of threads for rayon
    match rayon::ThreadPoolBuilder::new().num_threads(settings.threads as usize).build_global() {
        Ok(()) => {},
        Err(e) => {
            error!("Error while building thread pool: {e}");
            std::process::exit(exitcode::OSERR);
        }
    };

    // create the primary output folder
    info!("Creating output folder at {:?}...", settings.output_folder);
    match std::fs::create_dir_all(&settings.output_folder) {
        Ok(()) => {},
        Err(e) => {
            error!("Error while creating output folder: {e}");
            std::process::exit(exitcode::IOERR);
        }
    }

    // create a debug folder if specified, files might get created in sub-routines
    if let Some(debug_folder) = settings.debug_folder.as_ref() {
        info!("Creating debug folder at {debug_folder:?}...");
        match std::fs::create_dir_all(debug_folder) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while creating debug folder: {e}");
                std::process::exit(exitcode::IOERR);
            }
        }

        // save the CLI options
        let cli_json = debug_folder.join("cli_settings.json");
        info!("Saving CLI options to {cli_json:?}...");
        if let Err(e) = save_json(&settings, &cli_json) {
            error!("Error while saving CLI options: {e}");
            std::process::exit(exitcode::IOERR);
        }
    }

    // load the reference genome
    info!("Pre-loading reference genome into memory...");
    let reference_genome = match ReferenceGenome::from_fasta(&settings.reference_fn) {
        Ok(rg) => rg,
        Err(e) => {
            error!("Error while loading reference genome: {e:?}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // build the region iterator
    info!("Generating regions to compare...");
    let mut region_iter = match RegionIterator::new_compare_iterator(
        &settings.truth_vcf_filename,
        &settings.truth_sample,
        &settings.query_vcf_filename,
        &settings.query_sample,
        settings.regions.as_deref(),
        &reference_genome,
        settings.min_variant_gap
    ) {
        Ok(ri) => ri,
        Err(e) => {
            error!("Error while building region iterator: {e:#}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // check if we're in debug mode
    let skip_count = settings.skip_blocks;
    let take_count = settings.take_blocks;
    let debug_run: bool = if skip_count != 0 || take_count != usize::MAX {
        warn!("Debug run detected, disabling file finalizing steps.");
        warn!("Blocks to skip: {}", skip_count);
        warn!("Blocks to process: {}", take_count);
        true
    } else {
        false
    };

    if debug_run {
        warn!("Skip is enabled, output make be truncated.");
    }

    // prep any writer accumulators
    let mut summary_writer = SummaryWriter::new(settings.compare_label.clone());

    info!("Opening output VCF files...");
    let mut vcf_writer = match VariantCategorizer::new(
        &settings.truth_vcf_filename, settings.truth_sample.clone(),
        &settings.query_vcf_filename, settings.query_sample.clone(),
        &settings.output_folder) {
        Ok(vw) => vw,
        Err(e) => {
            error!("Error while building debug VCF writers: {e:#}");
            std::process::exit(exitcode::IOERR);
        }
    };

    let mut region_writer = settings.debug_folder.as_ref().map(|debug_path| {
        info!("Opening debug region writer file...");
        let out_fn = debug_path.join("region_summary.tsv.gz");
        match RegionSummaryWriter::new(&out_fn) {
            Ok(rsw) => rsw,
            Err(e) => {
                error!("Error while building region summary writer: {e:#}");
                std::process::exit(exitcode::IOERR);
            }
        }
    });

    let mut region_seq_writer = settings.debug_folder.as_ref().map(|debug_path| {
        info!("Opening debug sequence writer file...");
        let out_fn = debug_path.join("region_sequences.tsv.gz");
        match RegionSequenceWriter::new(&out_fn) {
            Ok(rsw) => rsw,
            Err(e) => {
                error!("Error while building region summary writer: {e:#}");
                std::process::exit(exitcode::IOERR);
            }
        }
    });

    if settings.threads > 1 && !debug_run {
        // we have parallelization and this is not a debug run, let's pre-load our variant for max efficiency
        info!("Pre-loading all variants...");
        if let Err(e) = region_iter.preload_all_variants() {
            error!("Error while loading all variants into memory: {e:#}");
            std::process::exit(exitcode::IOERR);
        }
    }

    // metrics to track as we iterate
    let mut joint_gt = SummaryMetrics::default();
    let mut joint_hap = SummaryMetrics::default();
    let mut joint_basepair = SummaryMetrics::default();
    let mut error_blocks = 0;
    let mut solved_blocks = 0;

    // first pre-load all the regions we are comparing; this is currently single-threaded
    let all_regions: Vec<CompareRegion> = match region_iter
        .skip(skip_count)
        .take(take_count)
        .map(|r| {
            match r {
                Ok(mr) => mr.try_into(), // convert the MultiRegion into a CompareRegion; this shouldn't ever fail given our CLI setup
                Err(e) => Err(e)
            }
        })
        .collect() {
            Ok(region_vec) => region_vec,
            Err(e) => {
                error!("Error while building regions: {e:#}");
                std::process::exit(exitcode::IOERR);
            }
        };
    info!("Region generation complete.");

    // run the parallel iterator to solve them
    let style = get_progress_style();
    info!("Comparing regions...");
    let mut all_results: Vec<(CompareRegion, Option<CompareBenchmark>)> = all_regions.into_par_iter()
        //.progress_count(num_regions as u64)
        .progress_with_style(style)
        .map(|region| {
            debug!("region = {region:?}");
            let comparison = match solve_compare_region(&region, &reference_genome, region_seq_writer.is_some()) {
                Ok(r) => Some(r),
                Err(e) => {
                    error!("Error while solving compare region #{} ({}): {e:#}", region.region_id(), region.coordinates());
                    None
                }
            };
            debug!("Result = {comparison:?}");
            (region, comparison)
        })
        .collect();

    // sort them by region ID
    all_results.sort_by_key(|(r, _c)| r.region_id());
    info!("Region comparisons complete, saving all outputs...");

    // now save all our result outputs
    let print_delta = 100000;
    for (region, opt_comparison) in all_results.into_iter() {
        if let Some(comparison) = opt_comparison {
            // writer updates
            summary_writer.add_comparison_benchmark(&comparison);
            if let Some(rsw) = region_writer.as_mut() {
                if let Err(e) = rsw.write_region_summary(&region, &comparison) {
                    error!("Error while writing region summary results: {e:#}");
                    std::process::exit(exitcode::IOERR);
                }
            }

            if let Some(rsw) = region_seq_writer.as_mut() {
                if let Err(e) = rsw.write_region_sequences(&region, &comparison) {
                    error!("Error while writing region sequences results: {e:#}");
                    std::process::exit(exitcode::IOERR);
                }
            }

            if let Err(e) = vcf_writer.write_results(&region, &comparison) {
                error!("Error while writing VCF results: {e:#}");
                std::process::exit(exitcode::IOERR);
            }

            // stats this loop tracks
            joint_gt += comparison.bm_gt();
            joint_hap += comparison.bm_hap();
            joint_basepair += comparison.bm_basepair();
            solved_blocks += 1;
        } else {
            error_blocks += 1;
        }

        if region.region_id() % print_delta == 0 {
            info!("Last written block: B#{} {}", region.region_id(), region.coordinates());
        }
    }

    info!("Joint GT: {joint_gt:?}");
    info!("\tRecall: {:?}", joint_gt.recall());
    info!("\tPrecision: {:?}", joint_gt.precision());
    info!("\tF1: {:?}", joint_gt.f1());
    info!("Joint Hap: {joint_hap:?}");
    info!("\tRecall: {:?}", joint_hap.recall());
    info!("\tPrecision: {:?}", joint_hap.precision());
    info!("\tF1: {:?}", joint_hap.f1());
    info!("Joint Basepair: {joint_basepair:?}");
    info!("\tRecall: {:?}", joint_basepair.recall());
    info!("\tPrecision: {:?}", joint_basepair.precision());
    info!("\tF1: {:?}", joint_basepair.f1());
    info!("Solved:error blocks: {solved_blocks} : {error_blocks}");

    // now write things
    let summary_fn = settings.output_folder.join("summary.tsv");
    info!("Saving output summary to {:?}...", summary_fn);
    if let Err(e) = summary_writer.write_summary(&summary_fn) {
        error!("Error while saving summary file: {e:#}");
        std::process::exit(exitcode::IOERR);
    }

    // index the VCFs
    info!("Indexing all output VCF files...");
    if let Err(e) = index_categorizer(&settings.output_folder, vcf_writer) {
        error!("Error while indexing VCF files: {e:#}");
        std::process::exit(exitcode::IOERR);
    }

    info!("Comparisons completed in {} seconds.", start_time.elapsed().as_secs_f64());
}

fn run_merge(settings: MergeSettings) {
    // start the timer
    let start_time = Instant::now();

    // set up logging before we check the other settings
    let filter_level: LevelFilter = match settings.verbosity {
        0 => LevelFilter::Info,
        1 => LevelFilter::Debug,
        _ => LevelFilter::Trace
    };
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    let settings = match check_merge_settings(settings) {
        Ok(s) => s,
        Err(e) => {
            error!("Error while verifying settings: {e:#}");
            std::process::exit(exitcode::CONFIG);
        }
    };

    // set up the number of threads for rayon
    match rayon::ThreadPoolBuilder::new().num_threads(settings.threads as usize).build_global() {
        Ok(()) => {},
        Err(e) => {
            error!("Error while building thread pool: {e}");
            std::process::exit(exitcode::OSERR);
        }
    };

    // create a debug folder if specified, files might get created in sub-routines
    if let Some(debug_folder) = settings.debug_folder.as_ref() {
        info!("Creating debug folder at {debug_folder:?}...");
        match std::fs::create_dir_all(debug_folder) {
            Ok(()) => {},
            Err(e) => {
                error!("Error while creating debug folder: {e}");
                std::process::exit(exitcode::IOERR);
            }
        }

        // save the CLI options
        let cli_json = debug_folder.join("cli_settings.json");
        info!("Saving CLI options to {cli_json:?}...");
        if let Err(e) = save_json(&settings, &cli_json) {
            error!("Error while saving CLI options: {e}");
            std::process::exit(exitcode::IOERR);
        }
    }

    // load the reference genome
    info!("Pre-loading reference genome into memory...");
    let reference_genome = match ReferenceGenome::from_fasta(&settings.reference_fn) {
        Ok(rg) => rg,
        Err(e) => {
            error!("Error while loading reference genome: {e:?}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // build the region iterator
    info!("Generating regions to merge...");
    let mut region_iter = match RegionIterator::new_merge_iterator(
        &settings.vcf_filenames,
        &settings.vcf_samples,
        settings.merge_regions.as_deref(),
        &reference_genome,
        settings.min_variant_gap
    ) {
        Ok(ri) => ri,
        Err(e) => {
            error!("Error while building region iterator: {e:#}");
            std::process::exit(exitcode::IOERR);
        }
    };

    // check if we're in debug mode
    let skip_count = settings.skip_blocks;
    let take_count = settings.take_blocks;
    let debug_run: bool = if skip_count != 0 || take_count != usize::MAX {
        warn!("Debug run detected, disabling file finalizing steps.");
        warn!("Blocks to skip: {}", skip_count);
        warn!("Blocks to process: {}", take_count);
        true
    } else {
        false
    };

    if debug_run {
        warn!("Skip is enabled, output make be truncated.");
    }

    if settings.threads > 1 && !debug_run {
        // we have parallelization and this is not a debug run, let's pre-load our variant for max efficiency
        info!("Pre-loading all variants...");
        if let Err(e) = region_iter.preload_all_variants() {
            error!("Error while loading all variants into memory: {e:#}");
            std::process::exit(exitcode::IOERR);
        }
    }

    // first pre-load all the regions we are comparing; this is currently single-threaded
    let all_regions: Vec<MultiRegion> = match region_iter
        .skip(skip_count)
        .take(take_count)
        .collect() {
            Ok(region_vec) => region_vec,
            Err(e) => {
                error!("Error while building regions: {e:#}");
                std::process::exit(exitcode::IOERR);
            }
        };
    info!("Region generation complete.");

    // build our merge configuration
    let merge_config = match MergeConfigBuilder::default()
        .no_conflict_enabled(settings.enable_no_conflict)
        .majority_voting_enabled(settings.enable_voting)
        .conflict_selection(settings.conflict_selection)
        .build() {
        Ok(mc) => mc,
        Err(e) => {
            error!("Error while building merge config: {e:?}");
            std::process::exit(exitcode::SOFTWARE);
        }
    };

    // run the parallel iterator to solve them
    let style = get_progress_style();
    info!("Merging regions...");
    let mut all_results: Vec<(MultiRegion, Option<MergeBenchmark>)> = all_regions.into_par_iter()
        //.progress_count(num_regions as u64)
        .progress_with_style(style)
        .map(|region| {
            debug!("region = {region:?}");
            let comparison = match solve_merge_region(&region, &reference_genome, merge_config) {
                Ok(r) => Some(r),
                Err(e) => {
                    error!("Error while solving merge region #{} ({}): {e:#}", region.region_id(), region.coordinates());
                    None
                }
            };
            debug!("Result = {comparison:?}");
            (region, comparison)
        })
        .collect();

    // sort them by region ID
    all_results.sort_by_key(|(r, _c)| r.region_id());
    info!("Region merging complete, saving all outputs...");

    // now save all our result outputs
    let mut merged_writer = match VariantMerger::new(
        &settings.vcf_filenames, settings.vcf_tags.clone(),
        &settings.output_vcf_folder, settings.vcf_samples[0].clone()
    ) {
        Ok(mw) => mw,
        Err(e) => {
            error!("Error while creating variant writer: {e:#}");
            std::process::exit(exitcode::IOERR);
        }
    };

    let mut summary_writer = if settings.output_summary_filename.is_some() {
        Some(MergeSummaryWriter::default())
    } else {
        None
    };

    // iterate over each output and save the relevant info
    let print_delta = 100000;
    let mut solved_blocks = 0;
    let mut error_blocks = 0;
    for (region, opt_benchmark) in all_results.into_iter() {
        if let Some(benchmark) = opt_benchmark {
            // writer updates
            if let Err(e) = merged_writer.write_results(&region, &benchmark) {
                error!("Error while writing merged VCF results: {e:#}");
                std::process::exit(exitcode::IOERR);
            }

            // stat update if we're tracking that
            if let Some(writer) = summary_writer.as_mut() {
                writer.add_merge_benchmark(&region, &benchmark);
            }

            // stats this loop tracks
            solved_blocks += 1;
        } else {
            error_blocks += 1;
        }

        if region.region_id() % print_delta == 0 {
            info!("Last written block: B#{} {}", region.region_id(), region.coordinates());
        }
    }
    info!("Solved:error blocks: {solved_blocks} : {error_blocks}");

    // now write things
    if let Some(summary_fn) = settings.output_summary_filename.as_deref() {
        info!("Saving output summary to {:?}...", summary_fn);
        if let Err(e) = summary_writer.unwrap().write_summary(summary_fn, &settings.vcf_tags) {
            error!("Error while saving summary file: {e:#}");
            std::process::exit(exitcode::IOERR);
        }
    }

    // index the VCFs
    info!("Indexing all merged VCF files...");
    if let Err(e) = index_merger(&settings.output_vcf_folder, merged_writer) {
        error!("Error while indexing merged VCF files: {e:#}");
        std::process::exit(exitcode::IOERR);
    }

    info!("Merge completed in {} seconds.", start_time.elapsed().as_secs_f64());
}

fn main() {
    let cli = get_cli();
    match cli.command {
        Commands::Compare(settings) => {
            run_compare(*settings);
        },
        Commands::Merge(settings) => {
            run_merge(*settings);
        }
    }

    info!("Process finished successfully.");
}
