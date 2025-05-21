
use flate2::write::GzEncoder;
use serde::Serialize;
use std::fs::File;
use std::path::Path;

use crate::data_types::compare_benchmark::CompareBenchmark;
use crate::data_types::compare_region::CompareRegion;

/// This is a wrapper for writing out summary stats to a file
pub struct RegionSequenceWriter {
    /// Handle on the writer
    csv_writer: csv::Writer<GzEncoder<File>>,
}

/// Contains all the data written to each row of our stats file
#[derive(Serialize)]
struct RegionSequenceRow {
    /// Unique region identifier
    region_id: u64,
    /// Coordinates for simple tracking
    coordinates: String,
    
    ref_seq: String,
    truth_seq1: String,
    truth_seq2: String,
    query_seq1: String,
    query_seq2: String
}

impl RegionSequenceRow {
    /// Creates a new row from labels and sequences
    pub fn new(region: &CompareRegion, benchmark: &CompareBenchmark) -> Self {
        let region_id = region.region_id();
        let coordinates = format!("{}", region.coordinates());
        let sequence_bundle = benchmark.sequence_bundle().unwrap();
        Self {
            region_id, coordinates,
            ref_seq: sequence_bundle.ref_seq.clone(),
            truth_seq1: sequence_bundle.truth_seq1.clone(),
            truth_seq2: sequence_bundle.truth_seq2.clone(),
            query_seq1: sequence_bundle.query_seq1.clone(),
            query_seq2: sequence_bundle.query_seq2.clone(),
        }
    }
}

impl RegionSequenceWriter {
    /// Creates a new writer to save sequences. The output will be tab-delimited and gzipped.
    /// # Arguments
    /// * `filename` - path to the filename that will get opened; expected to be .tsv.gz
    pub fn new(filename: &Path) -> anyhow::Result<Self> {
        let delimiter: u8 = b'\t';
        let gzip_writer = GzEncoder::new(
            File::create(filename)?,
            // default compression = 6; the "best" mode was 2x slow with very little gains
            flate2::Compression::default()
        );

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
    pub fn write_region_sequences(&mut self, region: &CompareRegion, comparison: &CompareBenchmark) -> csv::Result<()> {
        // GT level analysis
        let row = RegionSequenceRow::new(
            region, comparison
        );
        self.csv_writer.serialize(&row)?;
        Ok(())
    }
}
