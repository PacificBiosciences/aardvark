/*!
# Writers module
Contains the logic for writing the output files for the compare and merge commands.
*/
/// Wrapper that handles the parallel writing of all outputs at the compare step
pub mod compare_parallel;
/// Generates the summary file for merging
pub mod merge_summary;
/// Helper functions for indexing file
pub mod noodles_idx;
/// Generates a large sequence file; each line corresponds to a sub-region
pub mod region_sequence;
/// Generates the much larger region summary file; each line correspond to a sub-region
pub mod region_summary;
/// Generates the summary file
pub mod summary;
/// Generates the categorized variant files
pub mod variant_categorizer;
/// Generates the merged VCF from the merge process
pub mod variant_merger;
