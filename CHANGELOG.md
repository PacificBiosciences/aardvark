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
