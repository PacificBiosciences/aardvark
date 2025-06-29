
use itertools::Itertools;

/// Classification of output for merging.
/// In general, anything that is not "Different" will go into the outputs.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum MergeClassification {
    /// Indicates they are different
    Different,
    /// Indicates there is no conflict amongst the list indices; all others had no variants
    NoConflict { indices: Vec<usize> },
    /// Indicates that the listed indices agree, but others do not
    MajorityAgree { indices: Vec<usize> },
    /// Indicates a user has selected a single index to report for conflict
    ConflictSelection { index: usize },
    /// Indicates they are basepair-level identical
    BasepairIdentical
    // TODO: phasing identical if we add phasing
    // TODO: non-conflicting, with more info
}

impl std::fmt::Display for MergeClassification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let simple_form = self.simplify();
        let index_seq = match self {
            MergeClassification::Different |
            MergeClassification::BasepairIdentical => String::default(),
            MergeClassification::NoConflict { indices } |
            MergeClassification::MajorityAgree { indices } => {
                indices.iter()
                    .map(|i| format!("_{i}"))
                    .join("")
            },
            MergeClassification::ConflictSelection { index } => format!("_{index}"),
        };

        f.write_str(simple_form)?;
        f.write_str(&index_seq)
    }
}

impl MergeClassification {
    /// Helper function to convert to a simplified string representation
    pub fn simplify(&self) -> &str {
        match self {
            MergeClassification::Different => "different",
            MergeClassification::NoConflict { .. } => "no_conflict",
            MergeClassification::MajorityAgree { .. } => "majority",
            MergeClassification::ConflictSelection { .. } => "conflict_select",
            MergeClassification::BasepairIdentical => "identical",
        }
    }
}

/// Intended to capture all of the results from a comparison
#[derive(Debug)]
pub struct MergeBenchmark {
    /// Unique identifier for the comparison region
    region_id: u64,
    /// The classification of this region
    merge_classification: MergeClassification
}

impl MergeBenchmark {
    /// Constructor
    pub fn new(region_id: u64, merge_classification: MergeClassification) -> Self {
        Self {
            region_id,
            merge_classification
        }
    }

    // getters
    pub fn region_id(&self) -> u64 {
        self.region_id
    }

    pub fn merge_classification(&self) -> &MergeClassification {
        &self.merge_classification
    }
}
