
use anyhow::Context;
use indicatif::{MultiProgress, ProgressBar, ProgressIterator};
use log::error;
use std::sync::{mpsc, Arc};

use crate::data_types::compare_benchmark::CompareBenchmark;
use crate::data_types::compare_region::CompareRegion;
use crate::util::progress_bar::get_progress_style;
use crate::writers::region_sequence::RegionSequenceWriter;
use crate::writers::region_summary::RegionSummaryWriter;
use crate::writers::summary::SummaryWriter;
use crate::writers::variant_categorizer::{index_categorizer, VariantCategorizer};

/// Wrapper function that will basically generate a bunch of writer threads that each independently iterate over the results and perform file I/O.
/// This cuts file I/O time roughly in half if we generate all outputs.
/// If any Errors occur, unaffected files will continue to write, but this function will ultimately generate an Error also.
/// The function also generates a multi-progress bar for monitoring everything.
/// # Arguments
/// * `all_results` - All of the regions and results from the core analysis
/// * `summary_writer` - Collector of summary data, this is a reference since main will need it for printing afterwards
/// * `vcf_writer` - Saves all the output VCFs, this function also calls the indexing of the VCFs
/// * `region_writer` - Writes summary statistics per-region
/// * `region_seq_writer` - Writes the haplotype sequences generated per-region, this tends to use the most time
pub fn write_compare_outputs(
    all_results: Vec<(CompareRegion, Option<CompareBenchmark>)>,
    summary_writer: &mut SummaryWriter,
    mut truth_vcf_writer: VariantCategorizer,
    mut query_vcf_writer: VariantCategorizer,
    mut region_writer: Option<RegionSummaryWriter>,
    mut region_seq_writer: Option<RegionSequenceWriter>
) -> anyhow::Result<()> {
    // wrap the inputs in an Arc reference for reading in each thread
    let arc_all_results = Arc::new(all_results);

    // create a MultiProgress bar
    let multi_bar = MultiProgress::new();

    // create a receiver in case any of the sub-threads have Errors
    let (arc_tx, rx) = mpsc::channel();
    rayon::scope(|s| {
        // This loop just gathers the summary stats into the summary writer
        let all_results = arc_all_results.clone();
        let pb = ProgressBar::new(all_results.len() as u64)
            .with_style(get_progress_style())
            .with_message("Gathering summary stats...");
        let pb = multi_bar.add(pb);
        pb.tick();
        s.spawn(move|_| {
            for (_region, opt_comparison) in all_results.iter() {
                if let Some(comparison) = opt_comparison {
                    // summary update
                    summary_writer.add_comparison_benchmark(comparison);
                } else {
                    // one block failed, mark the error
                    summary_writer.inc_error_blocks();
                }
                pb.inc(1);
            }
            pb.finish_with_message("Summary stats complete.");
        });

        // this loop handles truth VCF writing
        let all_results = arc_all_results.clone();
        let tx = arc_tx.clone();
        let pb = ProgressBar::new(all_results.len() as u64)
            .with_style(get_progress_style())
            .with_message("Truth VCF writing...");
        let pb = multi_bar.add(pb);
        pb.tick();
        s.spawn(move|_| {
            let mut error = false;
            for (region, opt_comparison) in all_results.iter() {
                if let Some(comparison) = opt_comparison {
                    // Write VCF results
                    if let Err(e) = truth_vcf_writer.write_results(region, comparison)
                        .with_context(|| "Error while writing truth VCF results:") {
                        tx.send(e).unwrap();
                        pb.abandon_with_message("Truth VCF writing - ERROR");
                        error = true;
                        break;
                    }
                }
                pb.inc(1);
            }

            if !error {
                // index the VCFs
                pb.set_message("Truth VCF indexing...");
                if let Err(e) = index_categorizer(truth_vcf_writer)
                    .with_context(|| "Error while indexing truth VCF file:") {
                    tx.send(e).unwrap();
                    pb.abandon_with_message("Truth VCF indexing - ERROR");
                    error = true;
                }
            }

            if !error {
                pb.finish_with_message("Truth VCF writing complete.");
            }
        });

        // this loop handles query VCF writing
        let all_results = arc_all_results.clone();
        let tx = arc_tx.clone();
        let pb = ProgressBar::new(all_results.len() as u64)
            .with_style(get_progress_style())
            .with_message("Query VCF writing...");
        let pb = multi_bar.add(pb);
        pb.tick();
        s.spawn(move|_| {
            let mut error = false;
            for (region, opt_comparison) in all_results.iter() {
                if let Some(comparison) = opt_comparison {
                    // Write VCF results
                    if let Err(e) = query_vcf_writer.write_results(region, comparison)
                        .with_context(|| "Error while writing query VCF results:") {
                        tx.send(e).unwrap();
                        pb.abandon_with_message("Query VCF writing - ERROR");
                        error = true;
                        break;
                    }
                }
                pb.inc(1);
            }

            if !error {
                // index the VCFs
                pb.set_message("Query VCF indexing...");
                if let Err(e) = index_categorizer(query_vcf_writer)
                    .with_context(|| "Error while indexing query VCF file:") {
                    tx.send(e).unwrap();
                    pb.abandon_with_message("Query VCF indexing - ERROR");
                    error = true;
                }
            }

            if !error {
                pb.finish_with_message("Query VCF writing complete.");
            }
        });

        // this loop handles region writing
        if let Some(rsw) = region_writer.as_mut() {
            let all_results = arc_all_results.clone();
            let tx = arc_tx.clone();
            let pb = ProgressBar::new(all_results.len() as u64)
                .with_style(get_progress_style())
                .with_message("Region stats writing...");
            let pb = multi_bar.add(pb);
            pb.tick();
            s.spawn(move|_| {
                let mut error = false;
                for (region, opt_comparison) in all_results.iter() {
                    if let Some(comparison) = opt_comparison {
                        // Write region results
                        if let Err(e) = rsw.write_region_summary(region, comparison)
                            .with_context(|| "Error while writing region summary results:") {
                            tx.send(e).unwrap();
                            pb.abandon_with_message("Region stats writing - ERROR");
                            error = true;
                            break;
                        }
                    }
                    pb.inc(1);
                }
                if !error {
                    pb.finish_with_message("Region stats writing complete.");
                }
            });
        }

        // this loop handles the sequence writer
        if let Some(rsw) = region_seq_writer.as_mut() {
            let all_results = arc_all_results.clone();
            let tx = arc_tx.clone();
            let pb = ProgressBar::new(all_results.len() as u64)
                .with_style(get_progress_style())
                .with_message("Region sequence writing...");
            let pb = multi_bar.add(pb);
            pb.tick();
            s.spawn(move|_| {
                let mut error = false;
                for (region, opt_comparison) in all_results.iter() {
                    if let Some(comparison) = opt_comparison {
                        // Write region results
                        if let Err(e) = rsw.write_region_sequences(region, comparison)
                            .with_context(|| "Error while writing region sequences results:") {
                            tx.send(e).unwrap();
                            pb.abandon_with_message("Region sequence writing - ERROR");
                            error = true;
                            break;
                        }
                    }
                    pb.inc(1);
                }
                if !error {
                    pb.finish_with_message("Region sequence writing complete.");
                }
            });
        }
    });

    // drop this so rx closes
    std::mem::drop(arc_tx);

    // check if anything had an error
    let mut result = Ok(());
    for e in rx {
        error!("Error while writing files: {e:#}");
        result = Err(e);
    }
    result
}

/// Single threaded version of the above function. Clearly much simpler, but also slower if threads are available.
/// This is usually ~30 seconds or so, with the above version closer to ~10.
/// # Arguments
/// * `all_results` - All of the regions and results from the core analysis
/// * `summary_writer` - Collector of summary data, this is a reference since main will need it for printing afterwards
/// * `vcf_writer` - Saves all the output VCFs, this function also calls the indexing of the VCFs
/// * `region_writer` - Writes summary statistics per-region
/// * `region_seq_writer` - Writes the haplotype sequences generated per-region, this tends to use the most time
pub fn single_write_compare_outputs(
    all_results: Vec<(CompareRegion, Option<CompareBenchmark>)>,
    summary_writer: &mut SummaryWriter,
    mut truth_vcf_writer: VariantCategorizer,
    mut query_vcf_writer: VariantCategorizer,
    mut region_writer: Option<RegionSummaryWriter>,
    mut region_seq_writer: Option<RegionSequenceWriter>
) -> anyhow::Result<()> {
    for (region, opt_comparison) in all_results.into_iter()
        .progress_with_style(get_progress_style()) {
        if let Some(comparison) = opt_comparison {
            // summary update
            summary_writer.add_comparison_benchmark(&comparison);

            // Write VCF results
            truth_vcf_writer.write_results(&region, &comparison)
                .with_context(|| "Error while writing truth VCF results:")?;

            // Write VCF results
            query_vcf_writer.write_results(&region, &comparison)
                .with_context(|| "Error while writing query VCF results:")?;

            // Write region summary results
            if let Some(rsw) = region_writer.as_mut() {
                rsw.write_region_summary(&region, &comparison)
                    .with_context(|| "Error while writing region summary results:")?;
            }

            // Write region sequence results
            if let Some(rsw) = region_seq_writer.as_mut() {
                rsw.write_region_sequences(&region, &comparison)
                    .with_context(|| "Error while writing region sequences results:")?;
            }
        } else {
            // one block failed, mark the error
            summary_writer.inc_error_blocks();
        }
    }

    // finalize by indexing
    index_categorizer(truth_vcf_writer)?;
    index_categorizer(query_vcf_writer)?;

    Ok(())
}

