
use crate::data_types::grouped_metrics::GroupTypeMetrics;
use crate::data_types::summary_metrics::SummaryMetrics;
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

    /// Grouped metrics, which are the bulk of tracked statistics
    group_metrics: GroupTypeMetrics,

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
            group_metrics: Default::default(),
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
        // add to the group metrics
        self.group_metrics.add_truth_zygosity(variant, expected_zygosity_count, observed_zygosity_count)?;

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
        // add to the group stats
        self.group_metrics.add_query_zygosity(variant, expected_zygosity_count, observed_zygosity_count)?;

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
        self.group_metrics.add_basepair_metrics(basepair_metrics, variant_type);
    }

    /// Adds record-basepair metrics to the current tracking
    /// # Arguments
    /// * `record_bp_metrics` - the metrics to get added into the current values
    /// * `variant_type` - optional variant type these stats are added to; if none, it goes into the "All" grouping
    pub fn add_record_bp_metrics(&mut self, record_bp_metrics: SummaryMetrics, variant_type: Option<VariantType>) {
        self.group_metrics.add_record_bp_metrics(record_bp_metrics, variant_type);
    }

    /// Sets the values for a reverse benchmark. This is primarily to set query_tp and query_fp values from a reverse comparison.
    /// # Arguments
    /// * `other` - Results from a benchmark where truth and query have been swapped.
    pub fn add_swap_benchmark(&mut self, other: &Self) -> anyhow::Result<()> {
        // add the swapped group metrics first
        self.group_metrics.add_swap_benchmark(&other.group_metrics)?;

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

    // getters
    pub fn region_id(&self) -> u64 {
        self.region_id
    }

    pub fn group_metrics(&self) -> &GroupTypeMetrics {
        &self.group_metrics
    }

    pub fn truth_variant_data(&self) -> &[VariantMetrics] {
        &self.truth_variant_data
    }

    pub fn query_variant_data(&self) -> &[VariantMetrics] {
        &self.query_variant_data
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
