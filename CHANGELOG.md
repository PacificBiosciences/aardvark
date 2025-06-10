# v0.7.2
## Fixed
* Added changes to build script to enable bioconda building from source

# v0.7.1
## Changes
- Added two new `GT`-specific statistics to the summary files:
  - `truth_fn_gt` - The number of `truth_fn` where the allelic sequence was matched in both inputs, but with the wrong genotype (e.g., 1/1 in truth, 0/1 in query). This value is only populated if the `comparison` is `GT`.
  - `query_fp_gt` - The number of `query_fp` where the allelic sequence was matched in both inputs, but with the wrong genotype (e.g., 0/1 in truth, 1/1 in query). This value is only populated if the `comparison` is `GT`.
- Added parallelization to the writers for both `compare` and `merge`.

# v0.7.0
## Changes
- Added a new option to `compare` mode: `--stratifications`. If provided, this will post-annotate all regions by the provided stratifications and add additional rows to the output summary TSV, see documentation for more details.
- Added `query_total` to the main summary and region summary output files. This metric is `query_tp + query_fp`, and is analogous to `truth_total`.
- Replaced the info statements during file writing with a progress bar.

## Fixed
- Fixed an issue where "\*" ALT alleles were treated as alternate sequence. They are now ignored.

# v0.6.0
## Changes
- In `compare` mode, the `--confidence-regions` parameter has been replaced with `--regions` for consistency.
- The progress bar has been added to the parallelized variant loading step.

## Fixed
- Fixed an issue where variants that were hemizygous would be ignored entirely. They are now treated as homozygous variants for comparison.

# v0.5.4
## Changes
- Added a short-circuit in the `merge` routine that checks variant lengths prior to trying a full comparison. This reduces run-time significantly when larger events are unique to an input.
- Adds parallelization to the file loading for both `compare` and `merge`. The parallelization is across both file and chromosome for the typical process, significantly reducing initially variant loading times.

# v0.5.3
## Changes
- Added a new optional output file for merging (`--output-summary`) which contains statistics on the merge. See documentation for details.

# v0.5.2
## Changes
- Added a new file (`region_sequences.tsv.gz`) to the debug output for `aardvark compare`. This file contains the constructed haplotype sequences for each region.
- Replaced the region summary file in the debug output with a gzip-compressed version (`region_summary.tsv` -> `region_summary.tsv.gz`) for `aardvark compare`. In internal tests, this reduced the disk footprint by ~90% for this file.

## Fixed
- Updated `crossbeam-channel` for security update

# v0.5.1
Initial release.
