
use anyhow::{anyhow, bail, Context};
use coitrees::{COITree, Interval, IntervalTree};
use indexmap::IndexMap;
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use crate::parsing::noodles_helper::LoadedBed;
use crate::util::progress_bar::get_progress_style;

/// Wrapper for a collection of stratification BED files.
#[derive(Clone, Debug)]
pub struct Stratifications {
    /// Lookup from a label to a loaded stratification BED file
    datasets: IndexMap<String, StratData>
}

impl Stratifications {
    /// This will open a TSV file that is expected to have two columns and no header.
    /// The first column is a label, and the second column is a file path which may be relative.
    /// This will parse the provided file, open each relative file path, and construct the stratification search trees.
    pub fn from_tsv_batch(core_fn: &Path) -> anyhow::Result<Self> {
        let mut csv_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false) // no headers in the file, disable so we do not skip first row
            .from_path(core_fn)
            .with_context(|| format!("Error while opening {core_fn:?}:"))?;

        let core_folder = match core_fn.parent() {
            Some(parent) => parent.to_path_buf(),
            None => PathBuf::default()
        };

        // parse the TSV and fill in the loaded data
        let mut filenames: BTreeMap<String, PathBuf> = Default::default();
        for result in csv_reader.records() {
            let row = result.with_context(|| format!("Error while reading {core_fn:?}"))?;

            // make sure this is not a duplicate
            let label = row.get(0).ok_or(anyhow!("Missing label on row: {row:?}"))?;
            if filenames.contains_key(label) {
                bail!("Duplicate label found: {label}");
            }

            // now find the file
            let filename = row.get(1).ok_or(anyhow!("Missing filename on row: {row:?}"))?;
            let raw_path = PathBuf::from(filename);
            let full_path = if raw_path.has_root() {
                // this will not really happen often
                raw_path
            } else {
                // this is normal approach, with relative paths
                core_folder.join(raw_path)
            };

            // insert it
            assert!(filenames.insert(label.to_string(), full_path).is_none());
        }

        // par_iter with progress does not like the BTreeMap, so convert to Vec first
        let filenames: Vec<(String, PathBuf)> = filenames.into_iter().collect();

        // now lets actually load these files
        let style = get_progress_style();
        let datasets: Vec<(String, StratData)> = filenames.into_par_iter()
            .map(|(label, filepath)| {
                // load the BED file data
                let strat_data = StratData::from_bed(&filepath)
                    .with_context(|| format!("Error while loading {filepath:?}:"))?;
                Ok((label, strat_data))
            })
            .progress_with_style(style)
            .collect::<anyhow::Result<_>>()?;

        // reformat to an IndexMap
        let datasets = datasets.into_iter().collect();
        Ok(Self {
            datasets
        })
    }

    /// Returns the set of label indices that overlap the provided region.
    /// These are based on 0-based inclusive lookups.
    /// # Arguments
    /// * `chrom` - the chromosome interval
    /// * `first` - the first included base, 0-based
    /// * `last` - the last included base, 0-based
    pub fn overlaps(&self, chrom: &str, first: i32, last: i32) -> Vec<usize> {
        self.datasets.values().enumerate()
            .filter_map(|(l_index, strat_data)| {
                if strat_data.is_overlapping(chrom, first, last) {
                    Some(l_index)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Returns the set of label indices that fully contain the provided region.
    /// These are based on 0-based inclusive lookups.
    /// # Arguments
    /// * `chrom` - the chromosome interval
    /// * `first` - the first included base, 0-based
    /// * `last` - the last included base, 0-based
    pub fn containments(&self, chrom: &str, first: i32, last: i32) -> Vec<usize> {
        self.datasets.values().enumerate()
            .filter_map(|(l_index, strat_data)| {
                if strat_data.is_contained(chrom, first, last) {
                    Some(l_index)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Gets the labels for our dataset, which will be in the IndexMap order.
    pub fn labels(&self) -> Vec<String> {
        self.datasets.keys()
            .cloned().collect()
    }
}

#[derive(Clone)]
struct StratData {
    /// Lookup from a chromosome to a COITree, which has 0-based inclusive ranges
    lookup_trees: BTreeMap<String, COITree<(), usize>>
}

impl std::fmt::Debug for StratData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // COITree does not have Debug, so lets just convert it to a length for simplicity
        let lookup_counts: BTreeMap<String, usize> = self.lookup_trees.iter()
            .map(|(s, c)| {
                (s.clone(), c.len())
            })
            .collect();
        f.debug_struct("StratData").field("lookup_trees_len", &lookup_counts).finish()
    }
}

impl StratData {
    /// Loads a BED file and converts all the entries to the COI trees for lookup
    fn from_bed(bed_fn: &Path) -> anyhow::Result<Self> {
        // first, load the bed file into memory
        let loaded_bed = LoadedBed::preload_bed_file(bed_fn)?;
        
        // now we need to iterate over each chromosome and build a COITree from each
        let mut lookup_trees: BTreeMap<String, COITree<(), usize>> = Default::default();
        for (chrom, intervals) in loaded_bed.chrom_lookup().iter() {
            // convert the BED into COI intervals
            let coi_intervals: Vec<Interval<()>> = intervals.iter()
                .map(|i| {
                    // the positions are 1-based inclusive
                    // convert to 0-based inclusive
                    let start = i.start().ok_or(anyhow!("Missing start"))?.get() as i32 - 1;
                    let end = i.end().ok_or(anyhow!("Missing end"))?.get() as i32 - 1;
                    Ok(Interval::new(start, end, ()))
                })
                .collect::<anyhow::Result<_>>()?;
            
            // now put them into the tree and save it
            let coi_tree = COITree::new(&coi_intervals);
            assert!(lookup_trees.insert(chrom.clone(), coi_tree).is_none());
        }

        // we succeeded
        Ok(Self {
            lookup_trees
        })
    }

    /// Returns true if the provided interval overlaps at least one interval in this dataset.
    /// These are based on 0-based inclusive lookups.
    /// # Arguments
    /// * `chrom` - the chromosome interval
    /// * `first` - the first included base, 0-based
    /// * `last` - the last included base, 0-based
    fn is_overlapping(&self, chrom: &str, first: i32, last: i32) -> bool {
        match self.lookup_trees.get(chrom) {
            Some(coi_tree) => {
                coi_tree.query_count(first, last) > 0
            },
            None => false
        }
    }

    /// Returns true if the provided interval is fully contained in at least one interval in this dataset.
    /// These are based on 0-based inclusive lookups.
    /// # Arguments
    /// * `chrom` - the chromosome interval
    /// * `first` - the first included base, 0-based
    /// * `last` - the last included base, 0-based
    fn is_contained(&self, chrom: &str, first: i32, last: i32) -> bool {
        match self.lookup_trees.get(chrom) {
            Some(coi_tree) => {
                let mut included = false;
                coi_tree.query(first, last, |i| {
                    if i.first <= first && i.last >= last {
                        included = true;
                    }
                });
                included
            },
            None => false
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_example_stratification() {
        let strat_fn = PathBuf::from("test_data/example_stratification/strat.tsv");
        let strat = Stratifications::from_tsv_batch(&strat_fn).unwrap();

        // first some simple checks based on counts
        assert_eq!(strat.datasets.len(), 2);
        assert_eq!(strat.datasets[0].lookup_trees.get("mock").unwrap().len(), 2);
        assert_eq!(strat.datasets[0].lookup_trees.get("mock2").unwrap().len(), 1);
        assert_eq!(strat.datasets[1].lookup_trees.get("mock").unwrap().len(), 2);
        assert_eq!(strat.datasets[1].lookup_trees.get("mock2").unwrap().len(), 2);

        // now some functional checks - single bases
        assert_eq!(strat.labels(), vec!["example1".to_string(), "example2".to_string()]);
        assert!(strat.containments("mock", 9, 9).is_empty());
        assert_eq!(strat.containments("mock", 10, 10), vec![0]);
        assert_eq!(strat.containments("mock", 14, 14), vec![0]);
        assert_eq!(strat.containments("mock", 15, 15), vec![0, 1]);
        assert_eq!(strat.containments("mock", 20, 20), vec![1]);
        assert_eq!(strat.containments("mock", 25, 25), vec![0]);

        // now try some ranges
        assert_eq!(strat.containments("mock", 10, 19), vec![0]);
        assert_eq!(strat.overlaps("mock", 10, 19), vec![0, 1]);
        assert_eq!(strat.containments("mock", 15, 24), vec![1]);
        assert_eq!(strat.overlaps("mock", 15, 24), vec![0, 1]);
    }
}