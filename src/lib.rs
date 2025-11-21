
/*!
# Aardvark-bio
Aardvark-bio is the underlying library that supports the Aardvark command line tool.
The library provides the core functionality for comparing and merging variant calls.
The waffle solver contains the main entry point for the compare command, with example usage below:

## Example compare usage
```rust
use aardvark_bio::data_types::compare_region::CompareRegion;
use aardvark_bio::data_types::coordinates::Coordinates;
use aardvark_bio::data_types::phase_enums::PhasedZygosity;
use aardvark_bio::data_types::variants::Variant;
use aardvark_bio::data_types::summary_metrics::{SummaryGtMetrics, SummaryMetrics};
use aardvark_bio::data_types::variant_metrics::{VariantMetrics, VariantSource};
use aardvark_bio::waffle_solver::{CompareConfig, solve_compare_region};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;

// create a basic reference genome in memory
let mut reference_genome = ReferenceGenome::empty_reference();
reference_genome.add_contig(
    "mock_chr1".to_string(), "ACCGTTACCAGGACTTGACAAACCG"
).unwrap();

// create a problem to solve; this is a simple SNV in the middle of the sequence
let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
let truth_variants = vec![
    Variant::new_snv(0, 2, b"C".to_vec(), b"G".to_vec()).unwrap()
];
let truth_zygosity = vec![
    PhasedZygosity::PhasedHet10
];

// create a query set that is the same as the truth set
let query_variants = truth_variants.clone();
let query_zygosity = truth_zygosity.clone();

// put it all into a CompareRegion object for the solver
let problem = CompareRegion::new(
    0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
).unwrap();

// solve the problem
let compare_config = CompareConfig::default();
let result = solve_compare_region(&problem, &reference_genome, compare_config, None).unwrap();

// check the results
assert_eq!(result.total_ed(), 0); // should be exact paths
let group_metrics = result.group_metrics();
assert_eq!(*group_metrics.joint_metrics().gt(), SummaryGtMetrics::new(1, 0, 1, 0, 0, 0));
assert_eq!(*group_metrics.joint_metrics().hap(), SummaryMetrics::new(1, 0, 1, 0));
assert_eq!(*group_metrics.joint_metrics().basepair(), SummaryMetrics::new(2*1, 0, 2*1, 0));
assert_eq!(result.truth_variant_data(), &[VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap()]);
assert_eq!(result.query_variant_data(), &[VariantMetrics::new(VariantSource::Query, 1, 1).unwrap()]);

// check the sequences also
let sequence_bundle = result.sequence_bundle().unwrap();
assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
assert_eq!(&sequence_bundle.truth_seq1, "ACGGTTACCA");
assert_eq!(&sequence_bundle.query_seq1, "ACGGTTACCA");
assert_eq!(&sequence_bundle.truth_seq2, "ACCGTTACCA");
assert_eq!(&sequence_bundle.query_seq2, "ACCGTTACCA");
```
*/

/// Command line interface functionality that is specific to Aardvark
pub mod cli;
/// Contains various shared data types
pub mod data_types;
/// Contains Dynamic WFA implementations that enable fast quantification of edit distance between two dynamic sequences
pub mod dwfa;
/// Tooling for optimizing the exact GT category
pub mod exact_gt_optimizer;
/// Contains the core logic for identifying regions to merge
pub mod merge_solver;
/// Tooling for parsing input files into meaningful structs / data
pub mod parsing;
/// Optimizes query sequence relative to truth
pub mod query_optimizer;
/// Contains generic utility functions
pub mod util;
/// Contains the entry point for evaluating truth and query variants relative to each other
pub mod waffle_solver;
/// All output writers
pub mod writers;
