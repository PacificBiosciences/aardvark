
use noodles::bgzf;
use serde::Serialize;
use std::fs::File;
use std::path::Path;

use crate::data_types::compare_benchmark::CompareBenchmark;
use crate::data_types::compare_region::CompareRegion;

/// This is a wrapper for writing out summary stats to a file
pub struct RegionSequenceWriter {
    /// Handle on the writer
    csv_writer: csv::Writer<bgzf::io::MultithreadedWriter<File>>,
}

/// Contains all the data written to each row of our stats file.
/// Sequence fields borrow the bundle, or a `'static` `=column` alias.
#[derive(Serialize)]
struct RegionSequenceRow<'a> {
    /// Unique region identifier
    region_id: u64,
    /// Coordinates for simple tracking
    coordinates: &'a str,

    // each of the sequences that were generated
    ref_seq: &'a str,
    truth_seq1: &'a str,
    truth_seq2: &'a str,
    query_seq1: &'a str,
    query_seq2: &'a str
}

// these are labels for replicate sequences
const SAME_AS_REF: &str = "=ref_seq";
const SAME_AS_TRUTH1: &str = "=truth_seq1";
const SAME_AS_TRUTH2: &str = "=truth_seq2";
const SAME_AS_QUERY1: &str = "=query_seq1";

impl<'a> RegionSequenceRow<'a> {
    /// Creates a new row from labels and sequences.
    /// A column after `ref_seq` that repeats an earlier column is written as `=column`, using the earliest match.
    /// # Arguments
    /// * `region_id` - unique region identifier
    /// * `coordinates` - already formatted region coordinates
    /// * `benchmark` - comparison result holding the haplotype sequences
    /// # Errors
    /// * if `benchmark` has no sequence bundle
    pub fn new(region_id: u64, coordinates: &'a str, benchmark: &'a CompareBenchmark) -> anyhow::Result<Self> {
        let bundle = benchmark.sequence_bundle()
            .ok_or_else(|| anyhow::anyhow!("no sequence bundle for region {region_id} ({coordinates})"))?;
        let ref_seq = bundle.ref_seq.as_str();
        let truth1 = bundle.truth_seq1.as_str();
        let truth2 = bundle.truth_seq2.as_str();
        let query1 = bundle.query_seq1.as_str();
        let query2 = bundle.query_seq2.as_str();

        // point a repeat at the earliest earlier column
        let truth_seq1 = repeated_sequence(truth1, &[(SAME_AS_REF, ref_seq)]);
        let truth_seq2 = repeated_sequence(truth2, &[(SAME_AS_REF, ref_seq), (SAME_AS_TRUTH1, truth1)]);
        let query_seq1 = repeated_sequence(query1, &[(SAME_AS_REF, ref_seq), (SAME_AS_TRUTH1, truth1), (SAME_AS_TRUTH2, truth2)]);
        let query_seq2 = repeated_sequence(query2, &[(SAME_AS_REF, ref_seq), (SAME_AS_TRUTH1, truth1), (SAME_AS_TRUTH2, truth2), (SAME_AS_QUERY1, query1)]);

        Ok(Self {
            region_id, coordinates, ref_seq,
            truth_seq1, truth_seq2, query_seq1, query_seq2,
        })
    }
}

/// Returns `=column` for the earliest earlier column that `sequence` repeats.
/// # Arguments
/// * `sequence` - the haplotype being written
/// * `earlier` - earlier columns, as `("=column", sequence)`, in file order
fn repeated_sequence<'a>(sequence: &'a str, earlier: &[(&'static str, &str)]) -> &'a str {
    for &(label, value) in earlier {
        if sequence == value {
            return label;
        }
    }
    sequence
}

impl RegionSequenceWriter {
    /// Creates a new writer to save sequences. The output will be tab-delimited and gzipped.
    /// # Arguments
    /// * `filename` - path to the filename that will get opened; expected to be .tsv.gz
    /// * `threads` - worker threads for the gzip writing
    pub fn new(filename: &Path, threads: usize) -> anyhow::Result<Self> {
        let delimiter: u8 = b'\t';
        let w_threads = std::num::NonZeroUsize::new(threads.clamp(1, 32)).unwrap();
        let gzip_writer = bgzf::io::MultithreadedWriter::with_worker_count(w_threads, File::create(filename)?);
        let csv_writer= csv::WriterBuilder::new()
            .delimiter(delimiter)
            .from_writer(gzip_writer);
        Ok(Self {
            csv_writer
        })
    }

    /// Will write the summary out to the given file path
    /// # Arguments
    /// * `region` - the region of interest
    /// * `comparison` - results for this region
    /// # Errors
    /// * if `comparison` has no sequence bundle
    /// * if the row cannot be written
    pub fn write_region_sequences(&mut self, region: &CompareRegion, comparison: &CompareBenchmark) -> anyhow::Result<()> {
        // GT level analysis
        let coordinates = format!("{}", region.coordinates());
        let row = RegionSequenceRow::new(region.region_id(), &coordinates, comparison)?;
        self.csv_writer.serialize(&row)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data_types::compare_benchmark::{CompareBenchmark, SequenceBundle};
    use crate::data_types::coordinates::Coordinates;

    /// Builds one TSV data row for a region whose haplotypes are `sequences`.
    /// # Arguments
    /// * `sequences` - reference, truth 1, truth 2, query 1, and query 2
    fn written_row(sequences: [&str; 5]) -> String {
        let region = CompareRegion::new(
            7,
            Coordinates::new("chr1".to_string(), 10, 20),
            vec![], vec![],
            vec![], vec![],
        ).unwrap();
        let mut benchmark = CompareBenchmark::new(7, 0, 0);
        benchmark.add_sequence_bundle(SequenceBundle::new(
            sequences[0].to_string(),
            sequences[1].to_string(),
            sequences[2].to_string(),
            sequences[3].to_string(),
            sequences[4].to_string(),
        ));
        let coordinates = format!("{}", region.coordinates());
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_writer(vec![]);
        let row = RegionSequenceRow::new(region.region_id(), &coordinates, &benchmark).unwrap();
        writer.serialize(row).unwrap();
        String::from_utf8(writer.into_inner().unwrap()).unwrap()
    }

    /// A repeated haplotype is written as "=column" for the earliest earlier column it copies.
    #[test]
    fn test_repeated_haplotypes_collapse() {
        // truth 2 and query 1 copy truth 1; query 2 is different
        let homozygous_match = written_row(["ACGT", "AGGT", "AGGT", "AGGT", "ATTT"]);
        assert_eq!(homozygous_match, "7\tchr1:11-20\tACGT\tAGGT\t=truth_seq1\t=truth_seq1\tATTT\n");

        // query 2 copies truth 2, which itself differs from truth 1
        let query_copies_truth2 = written_row(["ACGT", "AGGT", "ACCT", "ATTT", "ACCT"]);
        assert_eq!(query_copies_truth2, "7\tchr1:11-20\tACGT\tAGGT\tACCT\tATTT\t=truth_seq2\n");

        // anything equal to the reference points at ref_seq, even when a later column also matches
        let copies_reference = written_row(["ACGT", "ACGT", "ACGT", "ATTT", "ATTT"]);
        assert_eq!(copies_reference, "7\tchr1:11-20\tACGT\t=ref_seq\t=ref_seq\tATTT\t=query_seq1\n");
    }

    /// A comparison without a sequence bundle is an error.
    #[test]
    fn test_missing_sequence_bundle_is_an_error() {
        let benchmark = CompareBenchmark::new(7, 0, 0);
        let result = RegionSequenceRow::new(7, "chr1:11-20", &benchmark);
        let Err(error) = result else {
            panic!("expected a missing sequence bundle");
        };
        assert_eq!(error.to_string(), "no sequence bundle for region 7 (chr1:11-20)");
    }
}
