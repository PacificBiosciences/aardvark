# v0.10.1
## Fixed
- Fixed a versioning issue with dependencies

# v0.10.0
## Changes
- GitHub release workflow updated to auto-publish to crates.io
- The library has been renamed to `aardvark-bio` for publishing to crates.io, binary name is unchanged

# v0.9.0
## Changes
- Adds the optional `RECORD_BP` comparison type to all summary files. Review our [methods](./docs/methods.md#record-basepair) for details on this metric. This metric is enabled with `--enable-record-basepair-metrics`.
- The `GT` and `BASEPAIR` metrics are recommend for most users. Thus, the `HAP` and `WEIGHTED_HAP` metrics are now disabled by default. Two new options have been added to enable these secondary metrics: `--enable-hapotype-metrics` and `--enabled-weighted-haplotype-metrics`.
- Refactored the way group metrics are internally represented for easier long-term maintenance
- Added new GitHub automations for developer ease-of-use

# v0.8.1
## Changes
- Adds automated workflow for static binary build and release. Aardvark functionality identical to v0.8.0.

# v0.8.0
## Changes
- Adds beta support for medium-sized variant types (<= 10 kbp):
  - `SvDeletion` and `SvInsertion` - Structural variants identified by the `SVTYPE` tag in the provided VCF files
  - `TrContraction` and `TrExpansion` - Tandem repeat variants identified by a `TRID` tag in the provided VCF files
  - Add `JointStructuralVariant` and `JointTandemRepeat` types, which are similar to `JointIndel`
  - Including these larger variant types increase the compute cost of Aardvark
- Adds additional heuristic limitation to reduce compute costs from larger variant types:
  - Adds a filter on variants larger than 10 kbp
  - Adds a short-circuit for GT candidate scoring such that any edits immediately ends the candidate exploration
  - Adds a limit on the number of candidates that can get generated before a sync point is identified in GT scoring, preventing exponentially blow-up of GT scoring when large variants mismatch
  - Adds filter on variants that are not sequence resolved (e.g., "<DEL>" entry is not supported)

## Fixed
- Updates `noodles` crate to resolve some parsing issues:
  - Fixed issue with parsing `SVLEN` field
  - Fixed a VCF parsing issue where absent chromosomes would error instead of gracefully returning no variants
- Updated build to `rust:1.88.0` to resolve some compilation issues

# v0.7.4
## Changes
- Adds a new scoring mode, `WEIGHTED_HAP`. This scoring mode is similar to `HAP` scoring, but variants are weighted by the number of changes between the REF and ALT sequences. For SNPs, the `HAP` and `WEIGHTED_HAP` scores should be identical since all SNPs have the same weight. For indels, each variant is effectively weighted by its length, so longer variants have an increased weight.

## Fixed
- Removed an exact-match shortcut in `compare` mode that would occasionally under-estimate the variant-level errors

# v0.7.3
## Changes
- Added REF/ALT trimming of identical tail basepair sequences, which resolves some edge-case variant conflicts in the T2T truth set. This minor change tends to slightly improve the overall accuracy by reducing variant conflicts, allowing for overlapping changes when the trimming removes previously conflicting bases. Some Indel variants may be classified differently compared to previous versions. See the output VCF files for trimmed representations. Examples:
  - AC->CC : this is trimmed to A->C and classified as a SNV now.
  - ACC->AC : this is trimmed to AC->A and classified as a Deletion now.
  - AC->ACC : this is trimmed to A->AC and classified as an Insertion now.
- Adds a new option `--disable-variant-trimming` to both `compare` and `merge` modes, which disable the above trimming behavior.
- Updated the install documentation to reflect bioconda support

# v0.7.2
## Fixed
- Added changes to build script to enable bioconda building from source

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
