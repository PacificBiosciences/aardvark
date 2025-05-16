
/// Command line interface functionality
pub mod cli;
/// Contains various shared data types
pub mod data_types;
/// Contains various Dynamic WFA implementations we leverage
pub mod dwfa;
/// Tooling for optimizing the exact GT category
pub mod exact_gt_optimizer;
/// Contains the core logic for identifying regions to merge
pub mod merge_solver;
/// Tooling for parsing input files into meaningful structs / data
pub mod parsing;
/// Optimizes query sequence relative to truth
pub mod query_optimizer;
/// Various utility functions that tend to be very generic
pub mod util;
/// Core logic for solving a problem using waffle_graph
pub mod waffle_solver;
/// All output writers
pub mod writers;
