
use std::ops::AddAssign;

/// High-level summary metrics we expect to use frequently
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SummaryMetrics {
    /// Number of truth entries found in the query
    pub truth_tp: u64,
    /// Number of truth entries missing in the query
    pub truth_fn: u64,
    /// Number of query entries that match truth
    pub query_tp: u64,
    /// Number of query entries that are not in truth
    pub query_fp: u64,
}

impl AddAssign for SummaryMetrics {
    // Enables += with stats
    fn add_assign(&mut self, rhs: Self) {
        self.truth_tp += rhs.truth_tp;
        self.truth_fn += rhs.truth_fn;
        self.query_tp += rhs.query_tp;
        self.query_fp += rhs.query_fp;
    }
}

impl SummaryMetrics {
    /// Constructor
    pub fn new(truth_tp: u64, truth_fn: u64, query_tp: u64, query_fp: u64) -> Self {
        Self {
            truth_tp,
            truth_fn,
            query_tp,
            query_fp
        }
    }

    /// Copies the truth values from `other` to our query values
    pub fn set_query_from_truth(&mut self, other: &Self) {
        self.query_tp = other.truth_tp;
        self.query_fp = other.truth_fn;
    }

    /// Copies the query values from `other` to our truth values
    pub fn set_truth_from_query(&mut self, other: &Self) {
        self.truth_tp = other.query_tp;
        self.truth_fn = other.query_fp;
    }

    /// Calculates recall if it can, which is relative to truth
    pub fn recall(&self) -> Option<f64> {
        let denom = self.truth_tp + self.truth_fn;
        if denom > 0 {
            Some(self.truth_tp as f64 / denom as f64)
        } else {
            None
        }
    }

    /// Calculates precision if it can, which is relative to query
    pub fn precision(&self) -> Option<f64> {
        let denom = self.query_tp + self.query_fp;
        if denom > 0 {
            Some(self.query_tp as f64 / denom as f64)
        } else {
            None
        }
    }

    /// Calculates F1 score if possible
    pub fn f1(&self) -> Option<f64> {
        if let (Some(recall), Some(precision)) = (self.recall(), self.precision()) {
            Some(2.0 * recall * precision / (recall + precision))
        } else {
            None
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use approx_eq::assert_approx_eq;

    #[test]
    fn test_scores() {
        let summary = SummaryMetrics { truth_tp: 10, truth_fn: 2, query_tp: 7, query_fp: 5};
        assert_approx_eq!(summary.recall().unwrap(), 10.0 / 12.0);
        assert_approx_eq!(summary.precision().unwrap(), 7.0 / 12.0);
        assert_approx_eq!(summary.f1().unwrap(), 2.0 * (10.0 / 12.0) * (7.0 / 12.0) / (17.0 / 12.0));
    }

    #[test]
    fn test_add_assign() {
        let mut summary = SummaryMetrics { truth_tp: 10, truth_fn: 2, query_tp: 3, query_fp: 4};
        let summary2 = SummaryMetrics { truth_tp: 3, truth_fn: 1, query_tp: 10, query_fp: 2};
        summary += summary2;
        assert_eq!(summary, SummaryMetrics { truth_tp: 13, truth_fn: 3, query_tp: 13, query_fp: 6})
    }
}
