
use anyhow::{bail, ensure};
use std::collections::BTreeMap;
use std::ops::AddAssign;

use crate::data_types::summary_metrics::{SummaryGtMetrics, SummaryMetrics};
use crate::data_types::variants::{Variant, VariantType};

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct GroupMetrics {
    /// Stores GT-level summary metrics
    gt: SummaryGtMetrics,
    /// Stores haplotype-level summary metrics; i.e., a homozygous TP counts for 2
    hap: SummaryMetrics,
    /// Stores weighted haplotype-level summary metrics; weights are based on ED between REF/ALT
    weighted_hap: SummaryMetrics,
    /// Stores alignment-level summary metrics (basepairs)
    basepair: SummaryMetrics,
    /// Stores the variant-level stats for GT comparison
    variant_gt: BTreeMap<VariantType, SummaryGtMetrics>,
    /// Stores the variant-level stats for hap comparison
    variant_hap: BTreeMap<VariantType, SummaryMetrics>,
    /// Stores the variant-level stats for weighted hap comparison
    variant_weighted_hap: BTreeMap<VariantType, SummaryMetrics>,
    /// Stores the variant-level stats for basepair comparison
    variant_basepair: BTreeMap<VariantType, SummaryMetrics>
}

impl GroupMetrics {
    /// Constructor
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        gt: SummaryGtMetrics, hap: SummaryMetrics, weighted_hap: SummaryMetrics, basepair: SummaryMetrics,
        variant_gt: BTreeMap<VariantType, SummaryGtMetrics>,
        variant_hap: BTreeMap<VariantType, SummaryMetrics>,
        variant_weighted_hap: BTreeMap<VariantType, SummaryMetrics>,
        variant_basepair: BTreeMap<VariantType, SummaryMetrics>
    ) -> Self {
        Self {
            gt, hap, weighted_hap, basepair,
            variant_gt, variant_hap, variant_weighted_hap, variant_basepair
        }
    }

    /// Adds a truth variant comparison to the tracking
    /// # Arguments
    /// * `variant` - the truth variant to add
    /// * `expected_zygosity_count` - number of times this allele was expected in the truth set; e.g. 0/1 => 1, 1/1 => 2
    /// * `observed_zygosity_count` - number of times this allele was identified in the best-matching query haplotype sequences
    pub fn add_truth_zygosity(&mut self, variant: &Variant, expected_zygosity_count: u8, observed_zygosity_count: u8) -> anyhow::Result<()> {
        // we are assuming that we're only looking at non-reference genotypes
        ensure!(expected_zygosity_count > 0, "No implementation for expecting all reference");

        // get the variant type for classification, and the variant metrics for mutation
        let variant_type = variant.variant_type();
        let variant_weight = variant.alt_ed()? as u64;
        let v_gt = self.variant_gt.entry(variant_type).or_default();
        let v_hap = self.variant_hap.entry(variant_type).or_default();
        let v_w_hap = self.variant_weighted_hap.entry(variant_type).or_default();

        /*
        In our modified exact setup, variants are either found or they are not; meaning that a false positive does not exist from a "truth" perspective.
        Instead, a false positive is just a query variant that is not found in truth (i.e., the inverse calculations).
        This means we should bail! if we find more observed alleles than expected.
         */

        // compare expected to observed
        match expected_zygosity_count.cmp(&observed_zygosity_count) {
            std::cmp::Ordering::Less => {
                // we found too many observed relative to expected, so add the extra to FP
                bail!("No implementation for truth false positives");
            },
            std::cmp::Ordering::Equal => {
                // they match, so add one per expected
                self.hap.truth_tp += expected_zygosity_count as u64;
                v_hap.truth_tp += expected_zygosity_count as u64;

                // weighted is same as hap but multiplied by the variant weight
                self.weighted_hap.truth_tp += expected_zygosity_count as u64 * variant_weight;
                v_w_hap.truth_tp += expected_zygosity_count as u64 * variant_weight;

                self.gt.summary_metrics.truth_tp += 1;
                v_gt.summary_metrics.truth_tp += 1;
            },
            std::cmp::Ordering::Greater => {
                // we found too few relative to expected, so add the missing to FN
                self.hap.truth_tp += observed_zygosity_count as u64;
                self.hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64;
                v_hap.truth_tp += observed_zygosity_count as u64;
                v_hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64;

                // weighted is same as hap but multiplied by the variant weight
                self.weighted_hap.truth_tp += observed_zygosity_count as u64 * variant_weight;
                self.weighted_hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64 * variant_weight;
                v_w_hap.truth_tp += observed_zygosity_count as u64 * variant_weight;
                v_w_hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64 * variant_weight;

                // expected > observed, this is a false negative
                self.gt.summary_metrics.truth_fn += 1;
                v_gt.summary_metrics.truth_fn += 1;
                if observed_zygosity_count > 0 {
                    // observed was >0, this is purely a GT difference like 1/1 -> 0/1
                    self.gt.truth_fn_gt += 1;
                    v_gt.truth_fn_gt += 1;
                }
            },
        };

        Ok(())
    }

    /// Adds a query variant comparison to the tracking
    /// # Arguments
    /// * `variant` - the query variant to add
    /// * `expected_zygosity_count` - number of times this allele was expected in the query set; e.g. 0/1 => 1, 1/1 => 2
    /// * `observed_zygosity_count` - number of times this allele was identified in the best-matching truth haplotype sequences
    pub fn add_query_zygosity(&mut self, variant: &Variant, expected_zygosity_count: u8, observed_zygosity_count: u8) -> anyhow::Result<()> {
        // we are assuming that we're only looking at non-reference genotypes
        ensure!(expected_zygosity_count > 0, "No implementation for expecting all reference");

        // get the variant type for classification, and the variant metrics for mutation
        let variant_type = variant.variant_type();
        let variant_weight = variant.alt_ed()? as u64;
        let v_gt = self.variant_gt.entry(variant_type).or_default();
        let v_hap = self.variant_hap.entry(variant_type).or_default();
        let v_w_hap = self.variant_weighted_hap.entry(variant_type).or_default();

        // normally, we would compare expected to observed, but this function is only used for exact matches currently
        ensure!(expected_zygosity_count == observed_zygosity_count, "No implementation for non-equal query zygosities");

        // they match, so add one per expected
        self.hap.query_tp += expected_zygosity_count as u64;
        v_hap.query_tp += expected_zygosity_count as u64;

        self.weighted_hap.query_tp += expected_zygosity_count as u64 * variant_weight;
        v_w_hap.query_tp += expected_zygosity_count as u64 * variant_weight;

        self.gt.summary_metrics.query_tp += 1;
        v_gt.summary_metrics.query_tp += 1;

        Ok(())
    }

    /// Adds basepair metrics to the current tracking
    /// # Arguments
    /// * `basepair_metrics` - the metrics to get added into the current values
    /// * `variant_type` - optional variant type these stats are added to; if none, it goes into the "All" grouping
    pub fn add_basepair_metrics(&mut self, basepair_metrics: SummaryMetrics, variant_type: Option<VariantType>) {
        if let Some(vt) = variant_type {
            let v_basepair = self.variant_basepair.entry(vt).or_default();
            *v_basepair += basepair_metrics;
        } else {
            // all variant category
            self.basepair += basepair_metrics;
        }
    }

    /// Sets the values for a reverse benchmark. This is primarily to set query_tp and query_fp values from a reverse comparison.
    /// # Arguments
    /// * `other` - Results from a benchmark where truth and query have been swapped.
    pub fn add_swap_benchmark(&mut self, other: &Self) -> anyhow::Result<()> {
        // set the best match values
        self.gt.set_query_from_truth(&other.gt);
        self.hap.set_query_from_truth(&other.hap);
        self.weighted_hap.set_query_from_truth(&other.weighted_hap);

        // bm_basepair and variant_basepair do not compute the inverse separately

        // set the per-variant values as well for both
        for (vt, metrics) in other.variant_gt.iter() {
            let entry = self.variant_gt.entry(*vt).or_default();
            entry.set_query_from_truth(metrics);
        }

        for (vt, metrics) in other.variant_hap.iter() {
            let entry = self.variant_hap.entry(*vt).or_default();
            entry.set_query_from_truth(metrics);
        }

        for (vt, metrics) in other.variant_weighted_hap.iter() {
            let entry = self.variant_weighted_hap.entry(*vt).or_default();
            entry.set_query_from_truth(metrics);
        }

        Ok(())
    }

    // getters
    pub fn gt(&self) -> &SummaryGtMetrics {
        &self.gt
    }

    pub fn hap(&self) -> &SummaryMetrics {
        &self.hap
    }

    pub fn weighted_hap(&self) -> &SummaryMetrics {
        &self.weighted_hap
    }

    pub fn basepair(&self) -> &SummaryMetrics {
        &self.basepair
    }

    pub fn variant_gt(&self) -> &BTreeMap<VariantType, SummaryGtMetrics> {
        &self.variant_gt
    }

    pub fn variant_hap(&self) -> &BTreeMap<VariantType, SummaryMetrics> {
        &self.variant_hap
    }

    pub fn variant_weighted_hap(&self) -> &BTreeMap<VariantType, SummaryMetrics> {
        &self.variant_weighted_hap
    }

    pub fn variant_basepair(&self) -> &BTreeMap<VariantType, SummaryMetrics> {
        &self.variant_basepair
    }
}

impl AddAssign<Self> for GroupMetrics {
    fn add_assign(&mut self, rhs: Self) {
        self.add_assign(&rhs);
    }
}

impl AddAssign<&Self> for GroupMetrics {
    // Enables += for the grouped statistics, this is basically calling AddAssign on SummaryMetrics a ton
    fn add_assign(&mut self, rhs: &Self) {
        // overall stats
        self.gt += rhs.gt;
        self.hap += rhs.hap;
        self.weighted_hap += rhs.weighted_hap;
        self.basepair += rhs.basepair;

        // variant type stats
        for (vt, metrics) in rhs.variant_gt.iter() {
            let entry = self.variant_gt.entry(*vt).or_default();
            *entry += *metrics;
        }

        for (vt, metrics) in rhs.variant_hap.iter() {
            let entry = self.variant_hap.entry(*vt).or_default();
            *entry += *metrics;
        }

        for (vt, metrics) in rhs.variant_weighted_hap.iter() {
            let entry = self.variant_weighted_hap.entry(*vt).or_default();
            *entry += *metrics;
        }

        for (vt, metrics) in rhs.variant_basepair.iter() {
            let entry = self.variant_basepair.entry(*vt).or_default();
            *entry += *metrics;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_add_assign() {
        let mut g1 = GroupMetrics::new(
            SummaryGtMetrics::new(1, 0, 1, 0, 0, 2),
            SummaryMetrics::new(1, 2, 1, 2),
            SummaryMetrics::new(3, 2, 1, 2),
            SummaryMetrics::new(3, 0, 1, 4),
            [(VariantType::Snv, SummaryGtMetrics::new(1, 2, 3, 4, 0, 0))].into_iter().collect(),
            Default::default(),
            Default::default(),
            [(VariantType::Insertion, SummaryMetrics::new(1, 2, 3, 4))].into_iter().collect(),
        );
        let g2 = GroupMetrics::new(
            SummaryGtMetrics::new(1, 2, 3, 4, 1, 0),
            SummaryMetrics::new(4, 3, 2, 1),
            SummaryMetrics::new(0, 1, 0, 2),
            SummaryMetrics::new(0, 1, 0, 2),
            Default::default(),
            Default::default(),
            Default::default(),
            [(VariantType::Insertion, SummaryMetrics::new(2, 2, 3, 4))].into_iter().collect(),
        );

        let expected_sum = GroupMetrics::new(
            SummaryGtMetrics::new(2, 2, 4, 4, 1, 2),
            SummaryMetrics::new(5, 5, 3, 3),
            SummaryMetrics::new(3, 3, 1, 4),
            SummaryMetrics::new(3, 1, 1, 6),
            [(VariantType::Snv, SummaryGtMetrics::new(1, 2, 3, 4, 0, 0))].into_iter().collect(),
            Default::default(),
            Default::default(),
            [(VariantType::Insertion, SummaryMetrics::new(3, 4, 6, 8))].into_iter().collect(),
        );

        // now add and test
        g1 += g2;
        assert_eq!(g1, expected_sum);
    }

    // the following are some very simple tests, higher tests are in the `waffle_solver.rs` file

    #[test]
    fn test_add_truth_zygosity() {
        let mut group_metrics = GroupMetrics::default();
        group_metrics.add_truth_zygosity(
            &Variant::new_insertion(0, 10, b"A".to_vec(), b"ACCC".to_vec()).unwrap(),
            2, 2
        ).unwrap();
        assert_eq!(group_metrics.gt, SummaryGtMetrics::new(1, 0, 0, 0, 0, 0));
        assert_eq!(group_metrics.hap, SummaryMetrics::new(2, 0, 0, 0));
        assert_eq!(group_metrics.weighted_hap, SummaryMetrics::new(6, 0, 0, 0));
        assert_eq!(group_metrics.basepair, SummaryMetrics::new(0, 0, 0, 0)); // basepair is unaffected
    }

    #[test]
    fn test_add_query_zygosity() {
        let mut group_metrics = GroupMetrics::default();
        group_metrics.add_query_zygosity(
            &Variant::new_insertion(0, 10, b"A".to_vec(), b"ACC".to_vec()).unwrap(),
            1, 1
        ).unwrap();
        assert_eq!(group_metrics.gt, SummaryGtMetrics::new(0, 0, 1, 0, 0, 0));
        assert_eq!(group_metrics.hap, SummaryMetrics::new(0, 0, 1, 0));
        assert_eq!(group_metrics.weighted_hap, SummaryMetrics::new(0, 0, 2, 0));
        assert_eq!(group_metrics.basepair, SummaryMetrics::new(0, 0, 0, 0)); // basepair is unaffected
    }

    #[test]
    fn test_add_basepair_metrics() {
        let mut group_metrics = GroupMetrics::default();
        group_metrics.add_basepair_metrics(SummaryMetrics::new(1, 2, 3, 4), None);
        assert_eq!(group_metrics.basepair, SummaryMetrics::new(1, 2, 3, 4));
    }

    #[test]
    fn test_add_swap_benchmark() {
        let mut group_metrics = GroupMetrics::default();
        let mut other_metrics = GroupMetrics::default();
        other_metrics.add_truth_zygosity(
            &Variant::new_insertion(0, 10, b"A".to_vec(), b"ACC".to_vec()).unwrap(),
            1, 1
        ).unwrap();
        group_metrics.add_swap_benchmark(&other_metrics).unwrap();

        // one truth TP becomes one query TP
        assert_eq!(group_metrics.gt, SummaryGtMetrics::new(0, 0, 1, 0, 0, 0));
        assert_eq!(group_metrics.hap, SummaryMetrics::new(0, 0, 1, 0));
        assert_eq!(group_metrics.weighted_hap, SummaryMetrics::new(0, 0, 2, 0));
        assert_eq!(group_metrics.basepair, SummaryMetrics::new(0, 0, 0, 0)); // basepair is unaffected
    }
}
