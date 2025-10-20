
use anyhow::{bail, ensure};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::ops::AddAssign;

use crate::data_types::summary_metrics::{SummaryGtMetrics, SummaryMetrics};
use crate::data_types::variants::{Variant, VariantType};

#[derive(Clone, Copy, Debug, Eq, PartialEq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum MetricsType {
    /// Genotype-level metrics
    #[serde(rename="GT")]
    Genotype,
    /// Haplotype-level metrics
    #[serde(rename="HAP")]
    Haplotype,
    /// Weighted haplotype-level metrics
    #[serde(rename="WEIGHTED_HAP")]
    WeightedHaplotype,
    /// Basepair-level metrics
    #[serde(rename="BASEPAIR")]
    Basepair,
    /// Record-basepair-level metrics
    #[serde(rename="RECORD_BP")]
    RecordBasepair
}

/// Stores metrics for a given region.
/// This includes joint metrics for all variants, and variant-specific metrics for each variant type.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct GroupTypeMetrics {
    /// Stores joint metrics for the entire group
    joint_metrics: GroupMetrics,
    /// Stores variant-level metrics for the group
    variant_metrics: BTreeMap<VariantType, GroupMetrics>
}

impl GroupTypeMetrics {
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

        // add to the joint metrics
        self.joint_metrics.add_truth_zygosity(variant_weight, expected_zygosity_count, observed_zygosity_count)?;

        // add to the variant metrics
        self.variant_metrics.entry(variant_type).or_default()
            .add_truth_zygosity(variant_weight, expected_zygosity_count, observed_zygosity_count)?;

        Ok(())
    }

    /// Adds a query variant comparison to the tracking
    /// # Arguments
    /// * `variant` - the query variant to add
    /// * `expected_zygosity_count` - number of times this allele was expected in the query set; e.g. 0/1 => 1, 1/1 => 2
    /// * `observed_zygosity_count` - number of times this allele was identified in the best-matching truth haplotype sequences
    pub fn add_query_zygosity(&mut self, variant: &Variant, expected_zygosity_count: u8, observed_zygosity_count: u8) -> anyhow::Result<()> {
        // get the variant type for classification, and the variant metrics for mutation
        let variant_type = variant.variant_type();
        let variant_weight = variant.alt_ed()? as u64;

        // add to the joint metrics
        self.joint_metrics.add_query_zygosity(variant_weight, expected_zygosity_count, observed_zygosity_count)?;

        // add to the variant metrics
        self.variant_metrics.entry(variant_type).or_default()
            .add_query_zygosity(variant_weight, expected_zygosity_count, observed_zygosity_count)?;

        Ok(())
    }

    /// Adds basepair metrics to the current tracking
    /// # Arguments
    /// * `basepair_metrics` - the metrics to get added into the current values
    /// * `variant_type` - optional variant type these stats are added to; if none, it goes into the "All" grouping
    pub fn add_basepair_metrics(&mut self, basepair_metrics: SummaryMetrics, variant_type: Option<VariantType>) {
        if let Some(vt) = variant_type {
            let v_metrics = self.variant_metrics.entry(vt).or_default();
            v_metrics.add_basepair_metrics(basepair_metrics);
        } else {
            // all variant category
            self.joint_metrics.add_basepair_metrics(basepair_metrics);
        }
    }

    /// Adds record-basepair metrics to the current tracking
    /// # Arguments
    /// * `record_bp_metrics` - the metrics to get added into the current values
    /// * `variant_type` - optional variant type these stats are added to; if none, it goes into the "All" grouping
    pub fn add_record_bp_metrics(&mut self, record_bp_metrics: SummaryMetrics, variant_type: Option<VariantType>) {
        if let Some(vt) = variant_type {
            let v_metrics = self.variant_metrics.entry(vt).or_default();
            v_metrics.add_record_bp_metrics(record_bp_metrics);
        } else {
            self.joint_metrics.add_record_bp_metrics(record_bp_metrics);
        }
    }

    /// Sets the values for a reverse benchmark. This is primarily to set query_tp and query_fp values from a reverse comparison.
    /// # Arguments
    /// * `other` - Results from a benchmark where truth and query have been swapped.
    pub fn add_swap_benchmark(&mut self, other: &Self) -> anyhow::Result<()> {
        self.joint_metrics.add_swap_benchmark(&other.joint_metrics)?;
        for (vt, metrics) in other.variant_metrics.iter() {
            let entry = self.variant_metrics.entry(*vt).or_default();
            entry.add_swap_benchmark(metrics)?;
        }
        Ok(())
    }

    // getters
    pub fn joint_metrics(&self) -> &GroupMetrics {
        &self.joint_metrics
    }

    pub fn variant_metrics(&self) -> &BTreeMap<VariantType, GroupMetrics> {
        &self.variant_metrics
    }
}

impl AddAssign<Self> for GroupTypeMetrics {
    fn add_assign(&mut self, rhs: Self) {
        self.add_assign(&rhs);
    }
}

impl AddAssign<&Self> for GroupTypeMetrics {
    fn add_assign(&mut self, rhs: &Self) {
        self.joint_metrics += &rhs.joint_metrics;
        for (vt, metrics) in rhs.variant_metrics.iter() {
            let entry = self.variant_metrics.entry(*vt).or_default();
            *entry += metrics;
        }
    }
}


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
    /// Stores record-basepair summary metrics ("raw" basepairs)
    record_bp: SummaryMetrics
}

impl GroupMetrics {
    /// Constructor
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        gt: SummaryGtMetrics,
        hap: SummaryMetrics,
        weighted_hap: SummaryMetrics,
        basepair: SummaryMetrics,
        record_bp: SummaryMetrics
    ) -> Self {
        Self {
            gt, hap, weighted_hap, basepair, record_bp
        }
    }

    /// Adds a truth variant comparison to the tracking
    /// # Arguments
    /// * `variant_weight` - the weight of the variant, this is used to calculate weighted metrics
    /// * `expected_zygosity_count` - number of times this allele was expected in the truth set; e.g. 0/1 => 1, 1/1 => 2
    /// * `observed_zygosity_count` - number of times this allele was identified in the best-matching query haplotype sequences
    pub fn add_truth_zygosity(&mut self, variant_weight: u64, expected_zygosity_count: u8, observed_zygosity_count: u8) -> anyhow::Result<()> {
        // we are assuming that we're only looking at non-reference genotypes
        ensure!(expected_zygosity_count > 0, "No implementation for expecting all reference");

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

                // weighted is same as hap but multiplied by the variant weight
                self.weighted_hap.truth_tp += expected_zygosity_count as u64 * variant_weight;

                self.gt.summary_metrics.truth_tp += 1;
            },
            std::cmp::Ordering::Greater => {
                // we found too few relative to expected, so add the missing to FN
                self.hap.truth_tp += observed_zygosity_count as u64;
                self.hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64;

                // weighted is same as hap but multiplied by the variant weight
                self.weighted_hap.truth_tp += observed_zygosity_count as u64 * variant_weight;
                self.weighted_hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64 * variant_weight;

                // expected > observed, this is a false negative
                self.gt.summary_metrics.truth_fn += 1;
                if observed_zygosity_count > 0 {
                    // observed was >0, this is purely a GT difference like 1/1 -> 0/1
                    self.gt.truth_fn_gt += 1;
                }
            },
        };

        Ok(())
    }

    /// Adds a query variant comparison to the tracking
    /// # Arguments
    /// * `variant_weight` - the weight of the variant, this is used to calculate weighted metrics
    /// * `expected_zygosity_count` - number of times this allele was expected in the query set; e.g. 0/1 => 1, 1/1 => 2
    /// * `observed_zygosity_count` - number of times this allele was identified in the best-matching truth haplotype sequences
    pub fn add_query_zygosity(&mut self, variant_weight: u64, expected_zygosity_count: u8, observed_zygosity_count: u8) -> anyhow::Result<()> {
        // we are assuming that we're only looking at non-reference genotypes
        ensure!(expected_zygosity_count > 0, "No implementation for expecting all reference");

        // normally, we would compare expected to observed, but this function is only used for exact matches currently
        ensure!(expected_zygosity_count == observed_zygosity_count, "No implementation for non-equal query zygosities");

        // they match, so add one per expected
        self.hap.query_tp += expected_zygosity_count as u64;

        self.weighted_hap.query_tp += expected_zygosity_count as u64 * variant_weight;

        self.gt.summary_metrics.query_tp += 1;

        Ok(())
    }

    /// Adds basepair metrics to the current tracking
    /// # Arguments
    /// * `basepair_metrics` - the metrics to get added into the current values
    pub fn add_basepair_metrics(&mut self, basepair_metrics: SummaryMetrics) {
        self.basepair += basepair_metrics;
    }

    /// Adds alt-basepair metrics to the current tracking
    /// # Arguments
    /// * `record_bp_metrics` - the metrics to get added into the current values
    pub fn add_record_bp_metrics(&mut self, record_bp_metrics: SummaryMetrics) {
        self.record_bp += record_bp_metrics;
    }

    /// Sets the values for a reverse benchmark. This is primarily to set query_tp and query_fp values from a reverse comparison.
    /// # Arguments
    /// * `other` - Results from a benchmark where truth and query have been swapped.
    pub fn add_swap_benchmark(&mut self, other: &Self) -> anyhow::Result<()> {
        // set the best match values
        self.gt.set_query_from_truth(&other.gt);
        self.hap.set_query_from_truth(&other.hap);
        self.weighted_hap.set_query_from_truth(&other.weighted_hap);

        // basepair and record_bp are computed separately, so we don't need to do anything here

        Ok(())
    }

    /// Gets the metrics for a given type
    /// # Arguments
    /// * `metrics_type` - the type of metrics to get
    pub fn get_metrics(&self, metrics_type: MetricsType) -> &SummaryMetrics {
        match metrics_type {
            MetricsType::Genotype => &self.gt.summary_metrics,
            MetricsType::Haplotype => &self.hap,
            MetricsType::WeightedHaplotype => &self.weighted_hap,
            MetricsType::Basepair => &self.basepair,
            MetricsType::RecordBasepair => &self.record_bp,
        }
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

    pub fn record_bp(&self) -> &SummaryMetrics {
        &self.record_bp
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
        self.record_bp += rhs.record_bp;
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
            SummaryMetrics::new(4, 0, 2, 4),
        );
        let g2 = GroupMetrics::new(
            SummaryGtMetrics::new(1, 2, 3, 4, 1, 0),
            SummaryMetrics::new(4, 3, 2, 1),
            SummaryMetrics::new(0, 1, 0, 2),
            SummaryMetrics::new(0, 1, 0, 2),
            SummaryMetrics::new(1, 1, 1, 2),
        );

        let expected_sum = GroupMetrics::new(
            SummaryGtMetrics::new(2, 2, 4, 4, 1, 2),
            SummaryMetrics::new(5, 5, 3, 3),
            SummaryMetrics::new(3, 3, 1, 4),
            SummaryMetrics::new(3, 1, 1, 6),
            SummaryMetrics::new(5, 1, 3, 6),
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
            3,
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
            2,
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
        group_metrics.add_basepair_metrics(SummaryMetrics::new(1, 2, 3, 4));
        assert_eq!(group_metrics.basepair, SummaryMetrics::new(1, 2, 3, 4));
    }

    #[test]
    fn test_add_record_bp_metrics() {
        let mut group_metrics = GroupMetrics::default();
        group_metrics.add_record_bp_metrics(SummaryMetrics::new(1, 2, 3, 4));
        assert_eq!(group_metrics.record_bp, SummaryMetrics::new(1, 2, 3, 4));
    }

    #[test]
    fn test_add_swap_benchmark() {
        let mut group_metrics = GroupMetrics::default();
        let mut other_metrics = GroupMetrics::default();
        other_metrics.add_truth_zygosity(
            2,
            1, 1
        ).unwrap();
        group_metrics.add_swap_benchmark(&other_metrics).unwrap();

        // one truth TP becomes one query TP
        assert_eq!(group_metrics.gt, SummaryGtMetrics::new(0, 0, 1, 0, 0, 0));
        assert_eq!(group_metrics.hap, SummaryMetrics::new(0, 0, 1, 0));
        assert_eq!(group_metrics.weighted_hap, SummaryMetrics::new(0, 0, 2, 0));
        assert_eq!(group_metrics.basepair, SummaryMetrics::new(0, 0, 0, 0)); // basepair is unaffected
        assert_eq!(group_metrics.record_bp, SummaryMetrics::new(0, 0, 0, 0)); // record_bp is unaffected
    }

    #[test]
    fn test_group_type_metrics_default() {
        let metrics = GroupTypeMetrics::default();
        assert_eq!(metrics.joint_metrics(), &GroupMetrics::default());
        assert!(metrics.variant_metrics().is_empty());
    }

    #[test]
    fn test_group_type_add_truth_zygosity() {
        let mut group_type_metrics = GroupTypeMetrics::default();

        // Create a mock variant using the correct constructor
        let variant = Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec()).unwrap();

        group_type_metrics.add_truth_zygosity(&variant, 2, 2).unwrap();

        // Check that joint metrics were updated
        let joint = group_type_metrics.joint_metrics();
        assert_eq!(joint.gt().summary_metrics.truth_tp, 1);
        assert_eq!(joint.hap().truth_tp, 2);
        assert_eq!(joint.weighted_hap().truth_tp, 2);

        // Check that variant-specific metrics were created and updated
        let variant_metrics = group_type_metrics.variant_metrics();
        assert_eq!(variant_metrics.len(), 1);
        assert!(variant_metrics.contains_key(&variant.variant_type()));

        let variant_specific = variant_metrics.get(&variant.variant_type()).unwrap();
        assert_eq!(variant_specific.gt().summary_metrics.truth_tp, 1);
        assert_eq!(variant_specific.hap().truth_tp, 2);
        assert_eq!(variant_specific.weighted_hap().truth_tp, 2);
    }

    #[test]
    fn test_group_type_add_query_zygosity() {
        let mut group_type_metrics = GroupTypeMetrics::default();

        let variant = Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec()).unwrap();

        group_type_metrics.add_query_zygosity(&variant, 1, 1).unwrap();

        // Check that joint metrics were updated
        let joint = group_type_metrics.joint_metrics();
        assert_eq!(joint.gt().summary_metrics.query_tp, 1);
        assert_eq!(joint.hap().query_tp, 1);
        assert_eq!(joint.weighted_hap().query_tp, 1);

        // Check that variant-specific metrics were created and updated
        let variant_metrics = group_type_metrics.variant_metrics();
        assert!(variant_metrics.contains_key(&variant.variant_type()));

        let variant_specific = variant_metrics.get(&variant.variant_type()).unwrap();
        assert_eq!(variant_specific.gt().summary_metrics.query_tp, 1);
        assert_eq!(variant_specific.hap().query_tp, 1);
        assert_eq!(variant_specific.weighted_hap().query_tp, 1);
    }

    #[test]
    fn test_add_basepair_metrics_with_variant_type() {
        let mut group_type_metrics = GroupTypeMetrics::default();
        let variant_type = VariantType::Snv;
        let basepair_metrics = SummaryMetrics::new(10, 5, 3, 2);

        group_type_metrics.add_basepair_metrics(basepair_metrics.clone(), Some(variant_type));

        // Check that variant-specific metrics were created and updated
        let variant_metrics = group_type_metrics.variant_metrics();
        assert_eq!(variant_metrics.len(), 1);
        assert!(variant_metrics.contains_key(&variant_type));

        let variant_specific = variant_metrics.get(&variant_type).unwrap();
        assert_eq!(variant_specific.basepair(), &basepair_metrics);

        // Joint metrics should be unchanged
        assert_eq!(group_type_metrics.joint_metrics().basepair(), &SummaryMetrics::default());
    }

    #[test]
    fn test_add_basepair_metrics_without_variant_type() {
        let mut group_type_metrics = GroupTypeMetrics::default();
        let basepair_metrics = SummaryMetrics::new(10, 5, 3, 2);

        group_type_metrics.add_basepair_metrics(basepair_metrics.clone(), None);

        // Check that joint metrics were updated
        let joint = group_type_metrics.joint_metrics();
        assert_eq!(joint.basepair(), &basepair_metrics);

        // Variant metrics should be unchanged
        assert!(group_type_metrics.variant_metrics().is_empty());
    }

    #[test]
    fn test_add_record_bp_metrics_with_variant_type() {
        let mut group_type_metrics = GroupTypeMetrics::default();
        let variant_type = VariantType::Insertion;
        let record_bp_metrics = SummaryMetrics::new(15, 7, 4, 1);

        group_type_metrics.add_record_bp_metrics(record_bp_metrics.clone(), Some(variant_type));

        // Check that variant-specific metrics were created and updated
        let variant_metrics = group_type_metrics.variant_metrics();
        assert_eq!(variant_metrics.len(), 1);
        assert!(variant_metrics.contains_key(&variant_type));

        let variant_specific = variant_metrics.get(&variant_type).unwrap();
        assert_eq!(variant_specific.record_bp(), &record_bp_metrics);

        // Joint metrics should be unchanged
        assert_eq!(group_type_metrics.joint_metrics().record_bp(), &SummaryMetrics::default());
    }

    #[test]
    fn test_add_record_bp_metrics_without_variant_type() {
        let mut group_type_metrics = GroupTypeMetrics::default();
        let record_bp_metrics = SummaryMetrics::new(15, 7, 4, 1);

        group_type_metrics.add_record_bp_metrics(record_bp_metrics.clone(), None);

        // Check that joint metrics were updated
        let joint = group_type_metrics.joint_metrics();
        assert_eq!(joint.record_bp(), &record_bp_metrics);

        // Variant metrics should be unchanged
        assert!(group_type_metrics.variant_metrics().is_empty());
    }

    #[test]
    fn test_group_type_add_swap_benchmark() {
        let mut group_type_metrics = GroupTypeMetrics::default();
        let mut other_metrics = GroupTypeMetrics::default();

        // Add some data to other_metrics
        let variant = Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec()).unwrap();
        other_metrics.add_truth_zygosity(&variant, 1, 1).unwrap();
        other_metrics.add_basepair_metrics(SummaryMetrics::new(5, 2, 1, 0), Some(variant.variant_type()));

        group_type_metrics.add_swap_benchmark(&other_metrics).unwrap();
        
        // Check that joint metrics were updated with swap benchmark
        let joint = group_type_metrics.joint_metrics();
        assert_eq!(joint.gt(), &SummaryGtMetrics::new(0, 0, 1, 0, 0, 0));
        assert_eq!(joint.hap(), &SummaryMetrics::new(0, 0, 1, 0));
        assert_eq!(joint.weighted_hap(), &SummaryMetrics::new(0, 0, 1, 0));

        // basepair and record_bp are computed separately, so these should still be empty
        assert!(joint.basepair().is_empty());
        assert!(joint.record_bp().is_empty());

        // Check that variant-specific metrics were updated
        let variant_metrics = group_type_metrics.variant_metrics();
        assert!(variant_metrics.contains_key(&variant.variant_type()));

        let variant_specific = variant_metrics.get(&variant.variant_type()).unwrap();
        assert_eq!(variant_specific.gt().summary_metrics.query_tp, 1);
        assert_eq!(variant_specific.hap().query_tp, 1);
    }

    #[test]
    fn test_group_type_metrics_add_assign() {
        let mut g1 = GroupTypeMetrics::default();
        let mut g2 = GroupTypeMetrics::default();

        // Add some data to both
        let variant1 = Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec()).unwrap();
        let variant2 = Variant::new_insertion(0, 200, b"C".to_vec(), b"CG".to_vec()).unwrap();
        let variant3 = Variant::new_snv(0, 300, b"T".to_vec(), b"A".to_vec()).unwrap();

        g1.add_truth_zygosity(&variant1, 1, 1).unwrap();
        g1.add_basepair_metrics(SummaryMetrics::new(10, 5, 3, 2), None);
        g1.add_basepair_metrics(SummaryMetrics::new(10, 5, 3, 2), Some(variant1.variant_type()));

        g2.add_truth_zygosity(&variant2, 2, 2).unwrap();
        g2.add_basepair_metrics(SummaryMetrics::new(15, 7, 4, 1), None);
        g2.add_basepair_metrics(SummaryMetrics::new(15, 7, 4, 1), Some(variant2.variant_type()));
        g2.add_truth_zygosity(&variant3, 1, 1).unwrap();
        g2.add_basepair_metrics(SummaryMetrics::new(15, 7, 4, 1), None);
        g2.add_basepair_metrics(SummaryMetrics::new(15, 7, 4, 1), Some(variant3.variant_type()));

        // Test AddAssign
        g1 += g2;

        // Check that joint metrics were combined
        let joint = g1.joint_metrics();
        assert_eq!(joint.gt(), &SummaryGtMetrics::new(3, 0, 0, 0, 0, 0));
        assert_eq!(joint.hap(), &SummaryMetrics::new(4, 0, 0, 0));
        assert_eq!(joint.basepair(), &SummaryMetrics::new(40, 19, 11, 4));
        assert_eq!(joint.record_bp(), &SummaryMetrics::new(0, 0, 0, 0));

        // Check that variant metrics were combined
        let variant_metrics = g1.variant_metrics();
        assert_eq!(variant_metrics.len(), 2);

        // Both variants should be present
        assert_eq!(variant_metrics.get(&VariantType::Snv).unwrap().basepair(), &SummaryMetrics::new(25, 12, 7, 3));
        assert_eq!(variant_metrics.get(&VariantType::Insertion).unwrap().basepair(), &SummaryMetrics::new(15, 7, 4, 1));
    }

    #[test]
    fn test_error_cases() {
        let mut group_type_metrics = GroupTypeMetrics::default();
        let variant = Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec()).unwrap();

        // Test that adding truth zygosity with 0 expected count fails
        let result = group_type_metrics.add_truth_zygosity(&variant, 0, 1);
        assert!(result.is_err());

        // Test that adding query zygosity with mismatched counts fails
        let result = group_type_metrics.add_query_zygosity(&variant, 2, 1);
        assert!(result.is_err());
    }

    #[test]
    fn test_multiple_variant_types() {
        let mut group_type_metrics = GroupTypeMetrics::default();

        let snp = Variant::new_snv(0, 100, b"A".to_vec(), b"T".to_vec()).unwrap();
        let insertion = Variant::new_insertion(0, 200, b"C".to_vec(), b"CGG".to_vec()).unwrap();
        let deletion = Variant::new_deletion(0, 300, b"AT".to_vec(), b"A".to_vec()).unwrap();

        // Add metrics for different variant types
        group_type_metrics.add_truth_zygosity(&snp, 1, 1).unwrap();
        group_type_metrics.add_truth_zygosity(&insertion, 2, 2).unwrap();
        group_type_metrics.add_truth_zygosity(&deletion, 1, 1).unwrap();

        // Check that all variant types are present
        let variant_metrics = group_type_metrics.variant_metrics();
        assert_eq!(variant_metrics.len(), 3);
        assert!(variant_metrics.contains_key(&snp.variant_type()));
        assert!(variant_metrics.contains_key(&insertion.variant_type()));
        assert!(variant_metrics.contains_key(&deletion.variant_type()));

        // Check that joint metrics accumulated all
        let joint = group_type_metrics.joint_metrics();
        assert_eq!(joint.gt(), &SummaryGtMetrics::new(3, 0, 0, 0, 0, 0));
        assert_eq!(joint.hap(), &SummaryMetrics::new(4, 0, 0, 0));
        assert_eq!(joint.weighted_hap(), &SummaryMetrics::new(6, 0, 0, 0));
        assert_eq!(joint.basepair(), &SummaryMetrics::new(0, 0, 0, 0));
        assert_eq!(joint.record_bp(), &SummaryMetrics::new(0, 0, 0, 0));
    }
}
