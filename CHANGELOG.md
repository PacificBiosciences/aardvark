# v0.5.2
## Changes
- Added a new file (`region_sequences.tsv.gz`) to the debug output for `aardvark compare`. This file contains the constructed haplotype sequences for each region.
- Replaced the region summary file in the debug output with a gzip-compressed version (`region_summary.tsv` -> `region_summary.tsv.gz`) for `aardvark compare`. In internal tests, this reduced the disk footprint by ~90% for this file.

## Fixed
- Updated `crossbeam-channel` for security update

# v0.5.1
Initial release.
