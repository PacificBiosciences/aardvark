

use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

use crate::data_types::compare_benchmark::CompareBenchmark;
use crate::data_types::summary_metrics::SummaryMetrics;
use crate::data_types::variants::VariantType;

pub const COMPARE_GT: &str = "GT";
pub const COMPARE_HAP: &str = "HAP";
pub const COMPARE_BASEPAIR: &str = "BASEPAIR";

/// This is a wrapper for writing out summary stats to a file
#[derive(Default)]
pub struct SummaryWriter {
    /// Comparison label to go on each row
    compare_label: String,
    /// Stores GT-level summary metrics
    gt_summary_metrics: SummaryMetrics,
    /// Stores haplotype-level summary metrics; i.e., a homozygous TP counts for 2
    hap_summary_metrics: SummaryMetrics,
    /// Stores alignment-level summary metrics (basepairs)
    basepair_summary_metrics: SummaryMetrics,
    /// Stores the variant-level stats for GT comparison
    variant_gt_summary_metrics: BTreeMap<VariantType, SummaryMetrics>,
    /// Stores the variant-level stats for hap comparison
    variant_hap_summary_metrics: BTreeMap<VariantType, SummaryMetrics>,
    /// Stores the variant-level stats for basepair comparison
    variant_basepair_summary_metrics: BTreeMap<VariantType, SummaryMetrics>
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct SummaryRow {
    /// User provided label
    compare_label: String,
    /// Comparison type
    comparison: String,
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
    /// Total number of true positives in query
    query_tp: u64,
    /// Total number of false positives
    query_fp: u64,
    /// Recall = truth.TP / (truth.TP+truth.FN)
    metric_recall: Option<f64>,
    /// Precision = query.TP / (query.TP + query.FP)
    metric_precision: Option<f64>,
    /// F1 = combination score of recall and precision
    metric_f1: Option<f64>
}

impl SummaryRow {
    /// Creates a new row from labels and summary metrics
    pub fn new(compare_label: String, variant_type: String, filter: String, comparison: String, metrics: &SummaryMetrics) -> Self {
        Self {
            compare_label,
            variant_type, filter, comparison,
            truth_total: metrics.truth_tp + metrics.truth_fn,
            truth_tp: metrics.truth_tp,
            truth_fn: metrics.truth_fn,
            query_tp: metrics.query_tp,
            query_fp: metrics.query_fp,
            metric_recall: metrics.recall(),
            metric_precision: metrics.precision(),
            metric_f1: metrics.f1(),
        }
    }
}

impl SummaryWriter {
    /// Creates a new writer to accumulate stats
    pub fn new(compare_label: String) -> Self {
        Self {
            compare_label,
            ..Default::default()
        }
    }

    /// Adds a set of metrics to our collection
    /// # Arguments
    /// * `comparison` - the results from our benchmarking
    pub fn add_comparison_benchmark(&mut self, comparison: &CompareBenchmark) {
        self.gt_summary_metrics += comparison.bm_gt();
        self.hap_summary_metrics += comparison.bm_hap();
        self.basepair_summary_metrics += comparison.bm_basepair();

        for (&variant_type, &variant_metrics) in comparison.variant_gt() {
            let entry = self.variant_gt_summary_metrics.entry(variant_type).or_default();
            *entry += variant_metrics;
        }

        for (&variant_type, &variant_metrics) in comparison.variant_hap() {
            let entry = self.variant_hap_summary_metrics.entry(variant_type).or_default();
            *entry += variant_metrics;
        }

        for (&variant_type, &variant_metrics) in comparison.variant_basepair() {
            let entry = self.variant_basepair_summary_metrics.entry(variant_type).or_default();
            *entry += variant_metrics;
        }
    }

    /// Will write the summary out to the given file path
    /// # Arguments
    /// * `filename` - the filename for the output (tsv/csv)
    pub fn write_summary(&mut self, filename: &Path) -> csv::Result<()> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let mut csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;

        // joint indel sub-categories
        let joint_label = "JointIndel".to_string();
        let joint_types = [VariantType::Insertion, VariantType::Deletion, VariantType::Indel];

        // GT level analysis
        write_category(
            &mut csv_writer, self.compare_label.clone(), "ALL".to_string(), "GT".to_string(),
            &self.gt_summary_metrics, &self.variant_gt_summary_metrics,
            joint_label.clone(), &joint_types
        )?;

        // Haplotype level analysis
        write_category(
            &mut csv_writer, self.compare_label.clone(), "ALL".to_string(), "HAP".to_string(),
            &self.hap_summary_metrics, &self.variant_hap_summary_metrics,
            joint_label.clone(), &joint_types
        )?;

        // Basepair level analysis
        write_category(
            &mut csv_writer, self.compare_label.clone(), "ALL".to_string(), "BASEPAIR".to_string(),
            &self.basepair_summary_metrics, &self.variant_basepair_summary_metrics,
            joint_label.clone(), &joint_types
        )?;

        // save everything
        csv_writer.flush()?;
        Ok(())
    }
}

/// Wrapper function for write out everything for a particular category
/// # Arguments
/// * `csv_writer` - the writer handle
/// * `compare_label` - user provided comparison label, fixed
/// * `filter` - pass through to filter field of row
/// * `comparison_type` - pass through to comparison in row
/// * `full_metrics` - the summary for all metric types
/// * `type_metrics` - variant-specific metrics
/// * `joint_label` - a joint label for a special row
/// * `joint_types` - the variant types that get added together for the joint row
#[allow(clippy::too_many_arguments)]
fn write_category(
    csv_writer: &mut csv::Writer<File>,
    compare_label: String, filter: String, comparison_type: String,
    full_metrics: &SummaryMetrics,
    type_metrics: &BTreeMap<VariantType, SummaryMetrics>,
    joint_label: String, joint_types: &[VariantType],
) -> csv::Result<()> {
    // write the row for all variants
    let all_row = SummaryRow::new(
        compare_label.clone(), "ALL".to_string(), filter.clone(), comparison_type.clone(), full_metrics
    );
    csv_writer.serialize(&all_row)?;

    // variant-specific metrics
    let mut joint_metrics = SummaryMetrics::default();
    for (variant_type, metrics) in type_metrics.iter() {
        let v_row = SummaryRow::new(
            compare_label.clone(), format!("{variant_type:?}"), filter.clone(), comparison_type.clone(), metrics
        );
        csv_writer.serialize(&v_row)?;

        if joint_types.contains(variant_type) {
            joint_metrics += *metrics;
        }
    }

    // we also have the joint indel row
    let joint_row = SummaryRow::new(
        compare_label.clone(), joint_label, filter.clone(), comparison_type.clone(), &joint_metrics
    );
    csv_writer.serialize(&joint_row)?;

    Ok(())
}
