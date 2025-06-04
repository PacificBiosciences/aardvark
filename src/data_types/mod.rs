
/// Contains comparison benchmark results for a given sub-unit of the full task
pub mod compare_benchmark;
/// Contains comparison regions which serve as a sub-unit for solving
pub mod compare_region;
/// Wrapper for coordinates with some additional functionalities
pub mod coordinates;
/// Wrapper containing SummaryMetrics for GT, HAP, BASEPAIR, and the variant-level metrics
pub mod grouped_metrics;
/// Contains merging benchmark results for VCF merging
pub mod merge_benchmark;
/// Contains comparison regions for multiple inputs, which serve as a sub-unit for solving
pub mod multi_region;
/// Different phasing-based enumerations
pub mod phase_enums;
/// Contains tracker for TP, FP, FN and derived metrics
pub mod summary_metrics;
/// Contains variant-specific metrics that will end up in a VCF debug file
pub mod variant_metrics;
/// Contains variant definition functionality and checks
pub mod variants;
