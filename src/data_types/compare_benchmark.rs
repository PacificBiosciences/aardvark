
use anyhow::{bail, ensure};
use std::collections::BTreeMap;

use crate::data_types::grouped_metrics::GroupMetrics;
use crate::data_types::summary_metrics::{SummaryGtMetrics, SummaryMetrics};
use crate::data_types::variant_metrics::{VariantMetrics, VariantSource};
use crate::data_types::variants::{Variant, VariantType};

/// Intended to capture all of the results from a comparison
#[derive(Debug)]
pub struct CompareBenchmark {
    /// Unique identifier for the comparison region
    region_id: u64,
    /// Best match, edit distance for hap1
    bm_edit_distance_h1: usize,
    /// Best match, edit distance for hap2
    bm_edit_distance_h2: usize,
    /// Best match genotype metrics
    bm_gt: SummaryGtMetrics,
    /// Best match haplotype metrics
    bm_hap: SummaryMetrics,
    /// Best match basepair metrics
    bm_basepair: SummaryMetrics,

    /// Variant genotype statistics
    variant_gt: BTreeMap<VariantType, SummaryGtMetrics>,
    /// Variant haplotype statistics
    variant_hap: BTreeMap<VariantType, SummaryMetrics>,
    /// Variant basepair statistics
    variant_basepair: BTreeMap<VariantType, SummaryMetrics>,

    // variant level statistics
    /// Variant metrics for truth variants
    truth_variant_data: Vec<VariantMetrics>,
    /// Variant metrics for query variant
    query_variant_data: Vec<VariantMetrics>,

    /// optional sequence bundle, this can use a lot of memory
    sequence_bundle: Option<SequenceBundle>,

    /// optional containment regions which share an index with the Stratifications, this often has 20+ entries in the GIAB test data
    containment_regions: Option<Vec<usize>>

    // TODO: phasing stats, if we care to report it
}

impl CompareBenchmark {
    /// Constructor
    /// # Arguments
    /// * `region_id` - unique region identifier
    /// * `bm_edit_distance_h1` - edit distance between the query haplotype 1 and the best path through the truth graph
    /// * `bm_edit_distance_h2` - edit distance between the query haplotype 2 and the best path through the truth graph
    pub fn new(
        region_id: u64, bm_edit_distance_h1: usize, bm_edit_distance_h2: usize,
    ) -> Self {
        Self {
            region_id,
            bm_edit_distance_h1,
            bm_edit_distance_h2,
            bm_gt: Default::default(),
            bm_hap: Default::default(),
            bm_basepair: Default::default(),
            variant_gt: Default::default(),
            variant_hap: Default::default(),
            variant_basepair: Default::default(),
            truth_variant_data: Default::default(),
            query_variant_data: Default::default(),
            sequence_bundle: None,
            containment_regions: None
        }
    }

    /// Adds a truth variant comparison to the tracking
    /// # Arguments
    /// * `variant` - the truth variant to add
    /// * `expected_zygosity_count` - number of times this allele was expected in the truth set; e.g. 0/1 => 1, 1/1 => 2
    /// * `observed_zygosity_count` - number of times this allele was identified in the best-matching query haplotype sequences
    pub fn add_truth_zygosity(&mut self, variant: &Variant, expected_zygosity_count: u8, observed_zygosity_count: u8) -> anyhow::Result<()> {
        // we are assuming that we're only looking at non-reference genotypes
        assert!(expected_zygosity_count > 0);

        // get the variant type for classification, and the variant metrics for mutation
        let variant_type = variant.variant_type();
        let v_gt = self.variant_gt.entry(variant_type).or_default();
        let v_hap = self.variant_hap.entry(variant_type).or_default();

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
                self.bm_hap.truth_tp += expected_zygosity_count as u64;
                v_hap.truth_tp += expected_zygosity_count as u64;

                self.bm_gt.summary_metrics.truth_tp += 1;
                v_gt.summary_metrics.truth_tp += 1;
            },
            std::cmp::Ordering::Greater => {
                // we found too few relative to expected, so add the missing to FN
                self.bm_hap.truth_tp += observed_zygosity_count as u64;
                self.bm_hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64;
                v_hap.truth_tp += observed_zygosity_count as u64;
                v_hap.truth_fn += (expected_zygosity_count - observed_zygosity_count) as u64;

                // expected > observed, this is a false negative
                self.bm_gt.summary_metrics.truth_fn += 1;
                v_gt.summary_metrics.truth_fn += 1;
                if observed_zygosity_count > 0 {
                    // observed was >0, this is purely a GT difference like 1/1 -> 0/1
                    self.bm_gt.truth_fn_gt += 1;
                    v_gt.truth_fn_gt += 1;
                }
            },
        };

        // add the variant metrics
        let variant_metrics = VariantMetrics::new(VariantSource::Truth, expected_zygosity_count, observed_zygosity_count)?;
        self.truth_variant_data.push(variant_metrics);

        Ok(())
    }

    /// Adds a query variant comparison to the tracking
    /// # Arguments
    /// * `variant` - the query variant to add
    /// * `expected_zygosity_count` - number of times this allele was expected in the query set; e.g. 0/1 => 1, 1/1 => 2
    /// * `observed_zygosity_count` - number of times this allele was identified in the best-matching truth haplotype sequences
    pub fn add_query_zygosity(&mut self, variant: &Variant, expected_zygosity_count: u8, observed_zygosity_count: u8) -> anyhow::Result<()> {
        // we are assuming that we're only looking at non-reference genotypes
        assert!(expected_zygosity_count > 0);

        // get the variant type for classification, and the variant metrics for mutation
        let variant_type = variant.variant_type();
        let v_gt = self.variant_gt.entry(variant_type).or_default();
        let v_hap = self.variant_hap.entry(variant_type).or_default();

        // normally, we would compare expected to observed, but this function is only used for exact matches currently
        ensure!(expected_zygosity_count == observed_zygosity_count, "No implementation for non-equal query zygosities");

        // they match, so add one per expected
        self.bm_hap.query_tp += expected_zygosity_count as u64;
        v_hap.query_tp += expected_zygosity_count as u64;

        self.bm_gt.summary_metrics.query_tp += 1;
        v_gt.summary_metrics.query_tp += 1;

        // add the variant metrics, since this is a query variant, toggle the source info
        let variant_metrics = VariantMetrics::toggle_source(
            &VariantMetrics::new(VariantSource::Truth, expected_zygosity_count, observed_zygosity_count)?
        );
        self.query_variant_data.push(variant_metrics);
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
            self.bm_basepair += basepair_metrics;
        }
    }

    /// Sets the values for a reverse benchmark. This is primarily to set query_tp and query_fp values from a reverse comparison.
    /// # Arguments
    /// * `other` - Results from a benchmark where truth and query have been swapped.
    pub fn add_swap_benchmark(&mut self, other: &Self) -> anyhow::Result<()> {
        // set the best match values
        self.bm_gt.set_query_from_truth(&other.bm_gt);
        self.bm_hap.set_query_from_truth(&other.bm_hap);

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

        // we need to add any truth variants from other to our query variants, and also invert them as we go
        self.query_variant_data.extend(
            other.truth_variant_data().iter()
                .map(|original| {
                    assert_eq!(original.source(), VariantSource::Truth);
                    VariantMetrics::toggle_source(original)
                })
        );

        Ok(())
    }

    /// Adds a sequence bundle to this return value
    pub fn add_sequence_bundle(&mut self, sequence_bundle: SequenceBundle) {
        self.sequence_bundle = Some(sequence_bundle);
    }

    /// Adds the annotated containment regions for this set
    pub fn add_containment_regions(&mut self, containment_regions: Vec<usize>) {
        self.containment_regions = Some(containment_regions);
    }

    // wrappers that are useful for summaries
    /// Total edit distance: if 0, then it means both haplotypes exactly match *a* path in the graph.
    /// Note that this does not mean the genotypes are perfect.
    pub fn total_ed(&self) -> usize {
        self.bm_edit_distance_h1 + self.bm_edit_distance_h2
    }

    /// Utility function that creates the grouped metrics for us by copying out of this benchmark.
    pub fn group_metrics(&self) -> GroupMetrics {
        GroupMetrics::new(
            self.bm_gt, self.bm_hap, self.bm_basepair,
            self.variant_gt.clone(), self.variant_hap.clone(), self.variant_basepair.clone()
        )
    }

    // getters
    pub fn region_id(&self) -> u64 {
        self.region_id
    }

    pub fn bm_hap(&self) -> SummaryMetrics {
        self.bm_hap
    }

    pub fn bm_gt(&self) -> SummaryGtMetrics {
        self.bm_gt
    }

    pub fn bm_basepair(&self) -> SummaryMetrics {
        self.bm_basepair
    }

    pub fn truth_variant_data(&self) -> &[VariantMetrics] {
        &self.truth_variant_data
    }

    pub fn query_variant_data(&self) -> &[VariantMetrics] {
        &self.query_variant_data
    }

    pub fn variant_gt(&self) -> &BTreeMap<VariantType, SummaryGtMetrics> {
        &self.variant_gt
    }

    pub fn variant_hap(&self) -> &BTreeMap<VariantType, SummaryMetrics> {
        &self.variant_hap
    }

    pub fn variant_basepair(&self) -> &BTreeMap<VariantType, SummaryMetrics> {
        &self.variant_basepair
    }

    pub fn sequence_bundle(&self) -> Option<&SequenceBundle> {
        self.sequence_bundle.as_ref()
    }

    pub fn containment_regions(&self) -> Option<&[usize]> {
        self.containment_regions.as_deref()
    }
}

/// Wrapper that just contains a ton of sequences
#[derive(Clone, Debug)]
pub struct SequenceBundle {
    pub ref_seq: String,
    pub truth_seq1: String,
    pub truth_seq2: String,
    pub query_seq1: String,
    pub query_seq2: String
}

impl SequenceBundle {
    /// Constructor
    pub fn new(ref_seq: String, truth_seq1: String, truth_seq2: String, query_seq1: String, query_seq2: String) -> Self {
        Self {
            ref_seq, truth_seq1, truth_seq2, query_seq1, query_seq2
        }
    }
}
