

use serde::Serialize;
use std::fs::File;
use std::path::Path;

use crate::data_types::compare_benchmark::CompareBenchmark;
use crate::data_types::compare_region::CompareRegion;
use crate::data_types::summary_metrics::SummaryMetrics;
use crate::writers::summary::{COMPARE_BASEPAIR, COMPARE_GT, COMPARE_HAP};

/// This is a wrapper for writing out summary stats to a file
pub struct RegionSummaryWriter {
    /// Handle on the writer
    csv_writer: csv::Writer<File>,
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct RegionSummaryRow {
    /// Unique region identifier
    region_id: u64,
    /// Coordinates for simple tracking
    coordinates: String,
    /// Comparison type
    comparison: String,
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

impl RegionSummaryRow {
    /// Creates a new row from labels and summary metrics
    pub fn new(region: &CompareRegion, comparison: String, metrics: &SummaryMetrics) -> Self {
        let region_id = region.region_id();
        let coordinates = format!("{}", region.coordinates());
        Self {
            region_id, coordinates, comparison,
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

impl RegionSummaryWriter {
    /// Creates a new writer to accumulate stats
    /// # Arguments
    /// * `filename` - path to the filename that will get opened, must be .csv/.tsv
    pub fn new(filename: &Path) -> csv::Result<Self> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;
        Ok(Self {
            csv_writer
        })
    }

    /// Will write the summary out to the given file path
    /// # Arguments
    /// * `region` - the region of interest
    /// * `comparison` - results for this region
    pub fn write_region_summary(&mut self, region: &CompareRegion, comparison: &CompareBenchmark) -> csv::Result<()> {
        // GT level analysis
        let gt_row = RegionSummaryRow::new(
            region, COMPARE_GT.to_string(), &comparison.bm_gt()
        );        
        self.csv_writer.serialize(&gt_row)?;

        // HAP level analysis
        let gt_row = RegionSummaryRow::new(
            region, COMPARE_HAP.to_string(), &comparison.bm_hap()
        );        
        self.csv_writer.serialize(&gt_row)?;

        // basepair level analysis
        let gt_row = RegionSummaryRow::new(
            region, COMPARE_BASEPAIR.to_string(), &comparison.bm_basepair()
        );        
        self.csv_writer.serialize(&gt_row)?;

        Ok(())
    }
}
