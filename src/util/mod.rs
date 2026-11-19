/*!
# Utility module
Contains generic utility functions that are not particularly specific to aardvark-bio.
*/

/// Exact edit-distance grid boundaries that can be extended as both sequences grow
pub mod boundary_alignment;
/// Lazily filled value compared only when both sides are set
pub mod cache;
/// Helper functions for read/writing JSON via serde
pub mod json_io;
/// Binary-heap entry ordered only by priority
pub mod priority_entry;
/// Helper functions for generating the progress bars
pub mod progress_bar;
/// Helper functions for sequence manipulation
pub mod sequence;
/// Functions for basic sequence alignment scoring
pub mod sequence_alignment;
