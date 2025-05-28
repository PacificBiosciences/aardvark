
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

use crate::data_types::merge_benchmark::{MergeBenchmark, MergeClassification};
use crate::data_types::multi_region::MultiRegion;
use crate::data_types::variants::VariantType;

/// we want a key that is (merge reason, variant type, vcf index)
type LookupKey = (MergeClassification, VariantType, usize);

/// This is a wrapper for writing out summary stats to a file
#[derive(Default)]
pub struct MergeSummaryWriter {
    /// Passing and failing counts for each lookup type
    pass_fail_counts: BTreeMap<LookupKey, (u64, u64)>,
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct MergeSummaryRow {
    /// The merge type
    merge_reason: String,
    /// The type of variant represented by this row
    variant_type: String,
    /// The vcf index
    vcf_index: usize,
    /// The vcf label
    vcf_label: String,
    /// Number of variants that passed
    pass_variants: u64,
    /// Number of variants that did not pass
    fail_variants: u64,
}

impl MergeSummaryRow {
    /// Creates a new row from labels and summary metrics
    pub fn new(
        merge_classification: &MergeClassification, variant_type: &VariantType, vcf_index: usize, 
        vcf_label: String, pass_variants: u64, fail_variants: u64
    ) -> Self {
        Self {
            merge_reason: merge_classification.to_string(),
            variant_type: format!("{variant_type:?}"),
            vcf_index,
            vcf_label,
            pass_variants,
            fail_variants
        }
    }
}

impl MergeSummaryWriter {
    /// Adds a set of metrics to our collection
    /// # Arguments
    /// * `comparison` - the results from our benchmarking
    pub fn add_merge_benchmark(&mut self, region: &MultiRegion, comparison: &MergeBenchmark) {
        // get the class, and figure out which indices are passing
        let merge_classification = comparison.merge_classification();
        let passing_indices = match merge_classification {
            MergeClassification::Different => vec![],
            MergeClassification::NoConflict { indices } |
            MergeClassification::MajorityAgree { indices } => indices.clone(),
            MergeClassification::ConflictSelection { index } => vec![*index],
            MergeClassification::BasepairIdentical => (0..region.variants().len()).collect(),
        };

        // now accumulate the pass/fail based on the input variants for the region
        for (i, variants) in region.variants().iter().enumerate() {
            let is_passing = passing_indices.contains(&i);
            for v in variants.iter() {
                let k = (merge_classification.clone(), v.variant_type(), i);
                let entry = self.pass_fail_counts.entry(k).or_default();
                if is_passing {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
    }

    /// Will write the summary out to the given file path
    /// # Arguments
    /// * `filename` - the filename for the output (tsv/csv)
    /// * `tags` - the set of tags for the input VCFs
    pub fn write_summary<T: std::fmt::Display>(&mut self, filename: &Path, tags: &[T]) -> csv::Result<()> {
        // modify the delimiter to "," if it ends with .csv
        let is_csv: bool = filename.extension().unwrap_or_default() == "csv";
        let delimiter: u8 = if is_csv { b',' } else { b'\t' };
        let mut csv_writer: csv::Writer<File> = csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_path(filename)?;

        // go through each entry in order and output the results
        for ((merge_class, variant_type, vcf_index), (pass_count, fail_count)) in self.pass_fail_counts.iter() {
            let row = MergeSummaryRow::new(
                merge_class, variant_type, *vcf_index,
                tags[*vcf_index].to_string(), *pass_count, *fail_count
            );
            csv_writer.serialize(&row)?;
        }

        // TODO: if we ever do JointIndel stats, we'll likely need to restructure to do the grouping properly

        // save everything
        csv_writer.flush()?;
        Ok(())
    }
}
