
use std::collections::BTreeMap;
use std::ops::AddAssign;

use crate::data_types::summary_metrics::SummaryMetrics;
use crate::data_types::variants::VariantType;

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct GroupMetrics {
    /// Stores GT-level summary metrics
    gt: SummaryMetrics,
    /// Stores haplotype-level summary metrics; i.e., a homozygous TP counts for 2
    hap: SummaryMetrics,
    /// Stores alignment-level summary metrics (basepairs)
    basepair: SummaryMetrics,
    /// Stores the variant-level stats for GT comparison
    variant_gt: BTreeMap<VariantType, SummaryMetrics>,
    /// Stores the variant-level stats for hap comparison
    variant_hap: BTreeMap<VariantType, SummaryMetrics>,
    /// Stores the variant-level stats for basepair comparison
    variant_basepair: BTreeMap<VariantType, SummaryMetrics>
}

impl GroupMetrics {
    /// Constructor
    pub fn new(
        gt: SummaryMetrics, hap: SummaryMetrics, basepair: SummaryMetrics,
        variant_gt: BTreeMap<VariantType, SummaryMetrics>,
        variant_hap: BTreeMap<VariantType, SummaryMetrics>,
        variant_basepair: BTreeMap<VariantType, SummaryMetrics>
    ) -> Self {
        Self {
            gt, hap, basepair,
            variant_gt, variant_hap, variant_basepair
        }
    }

    // getters
    pub fn gt(&self) -> &SummaryMetrics {
        &self.gt
    }

    pub fn hap(&self) -> &SummaryMetrics {
        &self.hap
    }

    pub fn basepair(&self) -> &SummaryMetrics {
        &self.basepair
    }

    pub fn variant_gt(&self) -> &BTreeMap<VariantType, SummaryMetrics> {
        &self.variant_gt
    }

    pub fn variant_hap(&self) -> &BTreeMap<VariantType, SummaryMetrics> {
        &self.variant_hap
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
            SummaryMetrics::new(1, 0, 1, 0),
            SummaryMetrics::new(1, 2, 1, 2),
            SummaryMetrics::new(3, 0, 1, 4),
            [(VariantType::Snv, SummaryMetrics::new(1, 2, 3, 4))].into_iter().collect(),
            Default::default(),
            [(VariantType::Insertion, SummaryMetrics::new(1, 2, 3, 4))].into_iter().collect(),
        );
        let g2 = GroupMetrics::new(
            SummaryMetrics::new(1, 2, 3, 4),
            SummaryMetrics::new(4, 3, 2, 1),
            SummaryMetrics::new(0, 1, 0, 2),
            Default::default(),
            Default::default(),
            [(VariantType::Insertion, SummaryMetrics::new(2, 2, 3, 4))].into_iter().collect(),
        );

        let expected_sum = GroupMetrics::new(
            SummaryMetrics::new(2, 2, 4, 4),
            SummaryMetrics::new(5, 5, 3, 3),
            SummaryMetrics::new(3, 1, 1, 6),
            [(VariantType::Snv, SummaryMetrics::new(1, 2, 3, 4))].into_iter().collect(),
            Default::default(),
            [(VariantType::Insertion, SummaryMetrics::new(3, 4, 6, 8))].into_iter().collect(),
        );

        // now add and test
        g1 += g2;
        assert_eq!(g1, expected_sum);
    }
}
