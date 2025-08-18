
use anyhow::ensure;
use indexmap::IndexMap;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

use crate::data_types::compare_benchmark::CompareBenchmark;
use crate::data_types::grouped_metrics::GroupMetrics;
use crate::data_types::summary_metrics::{SummaryGtMetrics, SummaryMetrics};
use crate::data_types::variants::VariantType;
use crate::parsing::stratifications::Stratifications;

pub const COMPARE_GT: &str = "GT";
pub const COMPARE_HAP: &str = "HAP";
pub const COMPARE_WEIGHTED_HAP: &str = "WEIGHTED_HAP";
pub const COMPARE_BASEPAIR: &str = "BASEPAIR";

/// This is a wrapper for writing out summary stats to a file
#[derive(Default)]
pub struct SummaryWriter {
    /// Comparison label to go on each row
    compare_label: String,
    /// Metrics for the "ALL" category
    all_metrics: GroupMetrics,
    /// Metrics for the stratifications
    strat_metrics: IndexMap<String, GroupMetrics>,
    /// Tracks the number of passing blocks added
    solved_blocks: u64,
    /// Tracks the number of blocks that failed
    error_blocks: u64
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct SummaryRow {
    /// User provided label
    compare_label: String,
    /// Comparison type, which will be one of GT, HAP, or BASEPAIR
    comparison: String,
    /// Region label, which is ALL or a label from stratification
    region_label: String,
    /// Any applied filters
    filter: String,
    /// The type of variant represented by this row
    variant_type: String,
    /// Total number of variants in the truth set
    // #[serde(rename="TRUTH.TOTAL")]
    truth_total: u64,
    /// Total number of true positives in truth
    truth_tp: u64,
    /// Total number of false negatives
    truth_fn: u64,
    /// Total number of variants in the query
    query_total: u64,
    /// Total number of true positives in query
    query_tp: u64,
    /// Total number of false positives
    query_fp: u64,
    /// Recall = truth.TP / (truth.TP+truth.FN)
    metric_recall: Option<f64>,
    /// Precision = query.TP / (query.TP + query.FP)
    metric_precision: Option<f64>,
    /// F1 = combination score of recall and precision
    metric_f1: Option<f64>,
    /// FN.GT = false negatives that are GT errors, so 1/1 -> 0/1; only present for GT type
    truth_fn_gt: Option<u64>,
    /// FP.GT = false positives that are GT errors, so 0/1 -> 1/1; only present for GT type
    query_fp_gt: Option<u64>
}

impl SummaryRow {
    /// Creates a new row from labels and summary metrics
    pub fn new(
        compare_label: String, comparison: String, region_label: String, filter: String, variant_type: String,
        metrics: &SummaryMetrics
    ) -> Self {
        Self {
            compare_label,
            comparison, region_label, filter, variant_type, 
            truth_total: metrics.truth_tp + metrics.truth_fn,
            truth_tp: metrics.truth_tp,
            truth_fn: metrics.truth_fn,
            query_total: metrics.query_tp + metrics.query_fp,
            query_tp: metrics.query_tp,
            query_fp: metrics.query_fp,
            metric_recall: metrics.recall(),
            metric_precision: metrics.precision(),
            metric_f1: metrics.f1(),
            truth_fn_gt: None,
            query_fp_gt: None
        }
    }

    /// Creates a new row for a GT entry
    pub fn new_gt(
        compare_label: String, comparison: String, region_label: String, filter: String, variant_type: String,
        metrics: &SummaryGtMetrics
    ) -> Self {
        Self {
            compare_label,
            comparison, region_label, filter, variant_type, 
            truth_total: metrics.summary_metrics.truth_tp + metrics.summary_metrics.truth_fn,
            truth_tp: metrics.summary_metrics.truth_tp,
            truth_fn: metrics.summary_metrics.truth_fn,
            query_total: metrics.summary_metrics.query_tp + metrics.summary_metrics.query_fp,
            query_tp: metrics.summary_metrics.query_tp,
            query_fp: metrics.summary_metrics.query_fp,
            metric_recall: metrics.summary_metrics.recall(),
            metric_precision: metrics.summary_metrics.precision(),
            metric_f1: metrics.summary_metrics.f1(),
            truth_fn_gt: Some(metrics.truth_fn_gt),
            query_fp_gt: Some(metrics.query_fp_gt)
        }
    }
}

impl SummaryWriter {
    /// Creates a new writer to accumulate stats
    /// # Arguments
    /// * `compare_label` - User-provided string that gets filled into a column
    /// * `stratifications` - Optional, if present the labels (and order) get copied
    pub fn new(compare_label: String, stratifications: Option<&Stratifications>) -> Self {
        let strat_metrics = if let Some(strat) = stratifications {
            let labels = strat.labels();
            labels.into_iter()
                .map(|l| (l, Default::default()))
                .collect()
        } else {
            // no stratifications
            Default::default()
        };

        Self {
            compare_label,
            all_metrics: Default::default(),
            strat_metrics,
            solved_blocks: 0,
            error_blocks: 0
        }
    }

    /// Adds a set of metrics to our collection
    /// # Arguments
    /// * `comparison` - the results from our benchmarking
    pub fn add_comparison_benchmark(&mut self, comparison: &CompareBenchmark) {
        // always add to ALL metrics
        let group_metrics = comparison.group_metrics();
        self.all_metrics += group_metrics;
        self.solved_blocks += 1;

        // add to any containment indices
        if let Some(contained_indices) = comparison.containment_regions() {
            for &ci in contained_indices.iter() {
                self.strat_metrics[ci] += group_metrics;
            }
        }
    }

    /// Increments the number of errors blocks we found
    pub fn inc_error_blocks(&mut self) {
        self.error_blocks += 1;
    }

    /// Will write the summary out to the given file path
    /// # Arguments
    /// * `filename` - the filename for the output (tsv/csv)
    pub fn write_summary(&mut self, filename: &Path) -> anyhow::Result<()> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let mut csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;

        // joint indel sub-categories
        let joint_indel_label = "JointIndel".to_string();
        let joint_indel_types = [VariantType::Insertion, VariantType::Deletion, VariantType::Indel];

        // joint tandem repeat sub-categories
        let joint_tr_label = "JointTandemRepeat".to_string();
        let joint_tr_types = [VariantType::TrExpansion, VariantType::TrContraction];

        // joint structural variant sub-categories
        let joint_sv_label = "JointStructuralVariant".to_string();
        let joint_sv_types = [
            VariantType::SvInsertion,
            VariantType::SvDeletion,
            VariantType::SvDuplication,
            VariantType::SvInversion,
            VariantType::SvBreakend,
        ];

        // construct the joint categories
        let joint_categories = [
            (joint_indel_label.clone(), joint_indel_types.as_slice()),
            (joint_sv_label.clone(), joint_sv_types.as_slice()),
            (joint_tr_label.clone(), joint_tr_types.as_slice()),
        ];

        // first, write the ALL group
        write_group(
            &mut csv_writer, self.compare_label.clone(), "ALL".to_string(), "ALL".to_string(),
            &self.all_metrics,
            &joint_categories
        )?;

        for (strat_label, strat_result) in self.strat_metrics.iter() {
            write_group(
                &mut csv_writer, self.compare_label.clone(), "ALL".to_string(), strat_label.clone(),
                strat_result,
                &joint_categories
            )?;
        }

        // save everything
        csv_writer.flush()?;
        Ok(())
    }

    // getters
    pub fn all_metrics(&self) -> &GroupMetrics {
        &self.all_metrics
    }

    pub fn solved_blocks(&self) -> u64 {
        self.solved_blocks
    }

    pub fn error_blocks(&self) -> u64 {
        self.error_blocks
    }
}

/// Wrapper function for write out everything for a particular group of statistics
/// # Arguments
/// * `csv_writer` - the writer handle
/// * `compare_label` - user provided comparison label, fixed
/// * `filter` - pass through to filter field of row, this is either "ALL" or one of the stratification labels
/// * `group_metrics` - the group metrics
/// * `joint_categories` - vector of (label, types) tuples for joint categories
fn write_group(
    csv_writer: &mut csv::Writer<File>,
    compare_label: String, filter: String, region_label: String,
    group_metrics: &GroupMetrics,
    joint_categories: &[(String, &[VariantType])],
) -> anyhow::Result<()> {
    // GT level analysis
    write_gt_category(
        csv_writer, compare_label.clone(), filter.clone(), COMPARE_GT.to_string(), region_label.clone(),
        group_metrics.gt(), group_metrics.variant_gt(),
        joint_categories
    )?;

    // Haplotype level analysis
    write_category(
        csv_writer, compare_label.clone(), filter.clone(), COMPARE_HAP.to_string(), region_label.clone(),
        group_metrics.hap(), group_metrics.variant_hap(),
        joint_categories
    )?;

    // Weighted Haplotype level analysis
    write_category(
        csv_writer, compare_label.clone(), filter.clone(), COMPARE_WEIGHTED_HAP.to_string(), region_label.clone(),
        group_metrics.weighted_hap(), group_metrics.variant_weighted_hap(),
        joint_categories
    )?;

    // Basepair level analysis
    write_category(
        csv_writer, compare_label, filter, COMPARE_BASEPAIR.to_string(), region_label,
        group_metrics.basepair(), group_metrics.variant_basepair(),
        joint_categories
    )?;
    Ok(())
}

/// Wrapper function for write out everything for a particular category
/// # Arguments
/// * `csv_writer` - the writer handle
/// * `compare_label` - user provided comparison label, fixed
/// * `filter` - pass through to filter field of row
/// * `comparison_type` - pass through to comparison in row
/// * `full_metrics` - the summary for all metric types
/// * `type_metrics` - variant-specific metrics
/// * `joint_categories` - vector of (label, types) tuples for joint categories
#[allow(clippy::too_many_arguments)]
fn write_category(
    csv_writer: &mut csv::Writer<File>,
    compare_label: String, filter: String, comparison_type: String, region_label: String,
    full_metrics: &SummaryMetrics,
    type_metrics: &BTreeMap<VariantType, SummaryMetrics>,
    joint_categories: &[(String, &[VariantType])],
) -> csv::Result<()> {
    // write the row for all variants
    let all_row = SummaryRow::new(
        compare_label.clone(), comparison_type.clone(), region_label.clone(), filter.clone(), "ALL".to_string(),
        full_metrics
    );
    csv_writer.serialize(&all_row)?;

    // variant-specific metrics
    for (variant_type, metrics) in type_metrics.iter() {
        if metrics.is_empty() {
            // do not write a row if there's nothing to write
            continue;
        }

        let v_row = SummaryRow::new(
            compare_label.clone(), comparison_type.clone(), region_label.clone(), filter.clone(), format!("{variant_type:?}"),
            metrics
        );
        csv_writer.serialize(&v_row)?;
    }

    // write joint category rows
    for (joint_label, joint_types) in joint_categories {
        let mut joint_metrics = SummaryMetrics::default();
        for variant_type in joint_types.iter() {
            if let Some(metrics) = type_metrics.get(variant_type) {
                joint_metrics += *metrics;
            }
        }

        // only write the row if there's something to write
        if !joint_metrics.is_empty() {
            let joint_row = SummaryRow::new(
                compare_label.clone(), comparison_type.clone(), region_label.clone(), filter.clone(), joint_label.clone(),
                &joint_metrics
            );
            csv_writer.serialize(&joint_row)?;
        }
    }

    Ok(())
}

/// Wrapper function for write out everything for a GT category
/// # Arguments
/// * `csv_writer` - the writer handle
/// * `compare_label` - user provided comparison label, fixed
/// * `filter` - pass through to filter field of row
/// * `comparison_type` - pass through to comparison in row
/// * `full_metrics` - the summary for all metric types
/// * `type_metrics` - variant-specific metrics
/// * `joint_categories` - vector of (label, types) tuples for joint categories
#[allow(clippy::too_many_arguments)]
fn write_gt_category(
    csv_writer: &mut csv::Writer<File>,
    compare_label: String, filter: String, comparison_type: String, region_label: String,
    full_metrics: &SummaryGtMetrics,
    type_metrics: &BTreeMap<VariantType, SummaryGtMetrics>,
    joint_categories: &[(String, &[VariantType])],
) -> anyhow::Result<()> {
    ensure!(comparison_type.as_str() == "GT", "write_gt_category requires a GT input");

    // write the row for all variants
    let all_row = SummaryRow::new_gt(
        compare_label.clone(), comparison_type.clone(), region_label.clone(), filter.clone(), "ALL".to_string(),
        full_metrics
    );
    csv_writer.serialize(&all_row)?;

    // variant-specific metrics
    for (variant_type, metrics) in type_metrics.iter() {
        let v_row = SummaryRow::new_gt(
            compare_label.clone(), comparison_type.clone(), region_label.clone(), filter.clone(), format!("{variant_type:?}"),
            metrics
        );
        csv_writer.serialize(&v_row)?;
    }

    // write joint category rows
    for (joint_label, joint_types) in joint_categories {
        let mut joint_metrics = SummaryGtMetrics::default();
        for variant_type in joint_types.iter() {
            if let Some(metrics) = type_metrics.get(variant_type) {
                joint_metrics += *metrics;
            }
        }

        // only write the row if there's something to write
        if !joint_metrics.summary_metrics.is_empty() {
            let joint_row = SummaryRow::new_gt(
                compare_label.clone(), comparison_type.clone(), region_label.clone(), filter.clone(), joint_label.clone(),
                &joint_metrics
            );
            csv_writer.serialize(&joint_row)?;
        }
    }

    Ok(())
}
