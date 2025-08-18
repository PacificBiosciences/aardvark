
use anyhow::{Context, anyhow};
use indexmap::IndexMap;
use log::debug;
use noodles::bed::io::reader::Builder as BedBuilder;
use noodles::bed::{io::Reader as BedReader, Record as BedRecord};
use noodles::core::region::Interval;
use std::io::BufReader;
use std::path::Path;

/// Wrapper function that handles both gzip compressed and uncompressed BED files
/// # Arguments
/// * `filename` - path to the .bed(.gz) file to open
pub fn open_bed_file(filename: &Path) -> anyhow::Result<BedReader<3, BufReader<Box<dyn std::io::Read>>>> {
    let is_compressed = match filename.extension() {
        Some(extension) => {
            extension == "gz"
        },
        None => false
    };

    let buf_reader: Box<dyn std::io::Read> = if is_compressed {
        #[allow(clippy::default_constructed_unit_structs)]
        let bgzf_reader = noodles::bgzf::io::reader::Builder::default()
            .build_from_path(filename)
            .with_context(|| format!("Error while loading {filename:?}:"))?;
        Box::new(bgzf_reader)
    } else {
        Box::new(std::fs::File::open(filename)?)
    };

    #[allow(clippy::default_constructed_unit_structs)]
    let bed_reader = BedBuilder::<3>::default()
        .build_from_reader(buf_reader);
    Ok(bed_reader)
}

/// A pre-loaded BED file where chromosome order is supported and the intervals are sorted.
pub struct LoadedBed {
    /// Map from chromosome to the sorted intervals
    chrom_lookup: IndexMap<String, Vec<Interval>>
}

impl LoadedBed {
    /// This will load an entire BED file into memory, preserving chromosome order but also sorting any intervals if they are not sorted already.
    /// # Arguments
    /// * `filename` - path to the .bed(.gz) file to open
    pub fn preload_bed_file(filename: &Path) -> anyhow::Result<Self> {
        // open the file
        debug!("Pre-loading {filename:?}...");
        let mut bed_handle = open_bed_file(filename)?;

        // we need to build more regions
        let mut record = BedRecord::<3>::default();
        let mut chrom_lookup: IndexMap<String, Vec<Interval>> = Default::default();
        while bed_handle.read_record(&mut record)? > 0 {
            let chrom = record.reference_sequence_name().to_string();
            let start = record.feature_start()
                .with_context(|| format!("Error while parsing start for record: {record:?}"))?;
            let end = record.feature_end()
                .unwrap_or(Err(std::io::Error::other("Missing end")))
                .with_context(|| format!("Error while parsing end for record: {record:?}"))?;
            let interval = Interval::from(start..=end);

            let entry = chrom_lookup.entry(chrom.clone()).or_default();
            entry.push(interval);
        }

        // sort all the intervals
        for (chrom, interval_set) in chrom_lookup.iter_mut() {
            let num_entries = interval_set.len();
            if !interval_set.is_sorted_by_key(|i| (i.start().unwrap(), i.end().unwrap())) {
                debug!("Sorting {num_entries} BED entries for {chrom}...");
                interval_set.sort_by_key(|i| (i.start().unwrap(), i.end().unwrap()));
            } else {
                debug!("Found {num_entries} sorted BED entries for {chrom}.");
            }
        }

        Ok(Self {
            chrom_lookup
        })
    }

    /// Get the entry as the specified index, returning both the chromosome (key) and the sorted interval set
    pub fn get_index(&self, index: usize) -> Option<(&String, &Vec<Interval>)> {
        self.chrom_lookup.get_index(index)
    }

    // getters
    pub fn chrom_lookup(&self) -> &IndexMap<String, Vec<Interval>> {
        &self.chrom_lookup
    }
}

/// This will open a VCF file and retrieve the sample name at the given index
/// # Arguments
/// * `vcf_fn` - the VCF filename to open
/// * `index` - the index of the sample to return; 0 = first sample
pub fn get_vcf_sample_name(vcf_fn: &Path, index: usize) -> anyhow::Result<String> {
    use noodles_util::variant::io::indexed_reader::Builder as VcfBuilder;

    // Open the VCF files
    let mut vcf_reader = VcfBuilder::default()
        .build_from_path(vcf_fn)
        .with_context(|| format!("Error while opening {vcf_fn:?} (or associated index):"))?;

    // get the headers also
    let vcf_header = vcf_reader.read_header()
        .with_context(|| format!("Error while reading header of {vcf_fn:?}:"))?;

    let sample_name = vcf_header.sample_names().get_index(0)
        .ok_or(anyhow!("Sample index {index} does not exist."))?
        .clone();

    Ok(sample_name)
}

// TODO: we might want some end-to-end tests here eventually, but we are getting identical results atm so I'm less concerned
