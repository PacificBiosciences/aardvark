/*!
# Exact GT Optimizer
Contains the logic for optimizing the exact GT category.
This is a branch-and-bound approach that will find the minimum number of ALT -> REF changes to make the sequences exactly match.
This approach assigns no partial credit, even if there are sequences "close" to each other.

## Example usage
```rust
use aardvark_bio::data_types::coordinates::Coordinates;
use aardvark_bio::data_types::variants::Variant;
use aardvark_bio::data_types::phase_enums::Allele;
use aardvark_bio::exact_gt_optimizer::optimize_gt_alleles;

// create a simple reference string with coordinates for the region
let reference = b"ACGTACGTACGT";
let coordinates = Coordinates::new("mock".to_string(), 0, reference.len() as u64);

// create a truth variants with allele state, copy for query (exact match)
let truth_variants = vec![
    // 0-based position is 5, C>G SNV
    Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap()
];
let truth_alleles = vec![
    Allele::Alternate
];
let query_variants = truth_variants.clone();
let query_alleles = truth_alleles.clone();

// run the optimizer
let result = optimize_gt_alleles(
    reference, &coordinates,
    &truth_variants, &truth_alleles,
    &query_variants, &query_alleles
).unwrap();

// check the results
assert!(result.is_exact_match());
assert_eq!(result.num_errors(), 0);
assert_eq!(result.truth_alleles(), &truth_alleles);
assert_eq!(result.query_alleles(), &query_alleles);
*/

use anyhow::bail;
use log::{debug, trace};
use priority_queue::PriorityQueue;
use std::cmp::Reverse;
use std::collections::BTreeMap;

use crate::data_types::coordinates::Coordinates;
use crate::data_types::phase_enums::Allele;
use crate::data_types::variants::Variant;
use crate::dwfa::haplotype_dwfa::{HapDWFAError, HaplotypeDWFA};
use crate::query_optimizer::order_variants;

/// Result from our exact sequence optimization routine
#[derive(Debug, Eq, PartialEq)]
pub struct OptimizedAlleles {
    /// Finalize zygosity of truth, forcing unphased into a phase
    truth_alleles: Vec<Allele>,
    /// Final zygosity of observations, forcing unphased into a phase
    query_alleles: Vec<Allele>,
    /// Number of ALT alleles that were converted to REF to make it match; shared with truth and query
    num_errors: usize
}

impl OptimizedAlleles {
    /// Constructor
    pub fn new(truth_alleles: Vec<Allele>, query_alleles: Vec<Allele>, num_errors: usize) -> Self {
        Self {
            truth_alleles,
            query_alleles,
            num_errors
        }
    }

    /// Returns true if there were no errors
    pub fn is_exact_match(&self) -> bool {
        self.num_errors == 0
    }

    // getters
    pub fn truth_alleles(&self) -> &[Allele] {
        &self.truth_alleles
    }

    pub fn query_alleles(&self) -> &[Allele] {
        &self.query_alleles
    }

    pub fn num_errors(&self) -> usize {
        self.num_errors
    }
}

/// Given two sets of variants and alleles, this will find the minimum number of ALT -> REF changes to make the sequences exactly match.
/// This restriction means there is no partial credit, and it is (in many ways) the lower bound on accuracy.
/// # Arguments
/// * `reference` - the full reference chromosome from coordinates.chrom
/// * `coordinates` - the coordinates of the region we are comparing
/// * `truth_variants` - the provided set of truth variants
/// * `truth_alleles` - the provided set of alleles
/// * `query_variants` - the set of query variants
/// * `query_alleles` - the provided set of alleles
/// # Errors
/// * if there are errors from DWFA extension or finalization
/// # Panics
/// * if an unsupported zygosity is provided
pub fn optimize_gt_alleles(
    reference: &[u8], coordinates: &Coordinates,
    truth_variants: &[Variant], truth_alleles: &[Allele],
    query_variants: &[Variant], query_alleles: &[Allele],
) -> anyhow::Result<OptimizedAlleles> {
    // do some sanity checks
    if truth_variants.len() != truth_alleles.len() {
        bail!("truth values must have equal length");
    }

    if query_variants.len() != query_alleles.len() {
        bail!("query values must have equal length");
    }

    debug!("Truth variants:");
    for (v, a) in truth_variants.iter().zip(truth_alleles.iter()) {
        debug!("\t{a:?} => {} {:?} {:?}", v.position(), v.allele0(), v.allele1());
    }
    debug!("Query variants:");
    for (v, a) in query_variants.iter().zip(query_alleles.iter()) {
        debug!("\t{a:?} => {} {:?} {:?}", v.position(), v.allele0(), v.allele1());
    }

    // first, figure out the order we are calculating variants in
    let all_variant_order: Vec<(usize, bool)> = order_variants(truth_variants, query_variants);
    let total_variant_count = all_variant_order.len();

    // create our initial empty root node
    let mut next_node_id = 0;
    let start = coordinates.start() as usize;
    let root_node = ExactMatchNode::new(next_node_id, start);
    assert!(root_node.is_exact_match());
    next_node_id += 1;

    // put it into the queue as our seed node
    let priority = root_node.priority();
    let mut pqueue: PriorityQueue<ExactMatchNode, NodePriority> = Default::default();
    pqueue.push(root_node, priority);

    // set our best to max ED with no result
    let mut best_error_count = usize::MAX;
    let mut best_result = None;

    // I don't think we need a bucket limit for this approach, the exact match enforces some limits that are implicit
    let max_per_bucket = usize::MAX;
    let mut bucket_counts = vec![0; total_variant_count];

    // track the minimum found allele synchronize point; anything before that has a cost >= the best at that point, and can be safely ignored
    let mut min_allele_sync = 0;

    // track the number of failures at a given sync point, and auto-fail variants that break a limit
    let mut failed_counts: std::collections::BTreeMap<usize, usize> = Default::default();
    let auto_fail_threshold = 500;
    let mut auto_fail_index = 0;
    let mut auto_fail_counts = 0;

    // track the time here, we have an automatic 5 min failure, but should not hit this with other limits
    let start = std::time::Instant::now();

    // loop while our queue is not empty
    while let Some((current_node, _priority)) = pqueue.pop() {
        if current_node.num_errors() >= best_error_count {
            // the cost of this node is not better than something already found, so skip it
            continue;
        }

        if start.elapsed().as_secs() > 300 {
            bail!("300 second time limit reached");
        }

        // check which variant we are on
        let order_index = current_node.set_alleles();
        if order_index == total_variant_count {
            // we did every variant, time to finalize the node
            let mut final_node = current_node;
            final_node.finalize_dwfas(reference, coordinates.end() as usize)?;

            // if our cost is less, this is a better solution
            let final_errors = final_node.num_errors();
            if final_node.is_exact_match() && final_errors < best_error_count {
                best_error_count = final_errors;
                best_result = Some(final_node);
            }
            continue;
        }

        if order_index < min_allele_sync {
            // we have a node that is less than the minimum found allele sync, ignore it
            continue;
        }

        // make sure we didn't fill our quota for this bucket size
        if bucket_counts[order_index] >= max_per_bucket {
            continue;
        }
        bucket_counts[order_index] += 1;

        // if we know this node is on the best path, then we can clear out whatever is in the queue since it is irrelevant
        if current_node.is_synchronized() {
            assert!(order_index >= min_allele_sync);
            min_allele_sync = order_index;

            // reset the auto-fail counts and update which index will auto-fail next
            auto_fail_counts = 0;
            auto_fail_index = min_allele_sync;

            for (i, c) in failed_counts.iter() {
                debug!("\t{i} ({:?}) => {c}", all_variant_order[*i]);
            }
            failed_counts.clear();

            debug!("New min_allele_sync={min_allele_sync}");
            debug!("Sync stats: {} / {} => {} / {}",
                current_node.hap_dwfa().truth_haplotype().alleles().len(),
                current_node.hap_dwfa().query_haplotype().alleles().len(),
                current_node.hap_dwfa().truth_haplotype().sequence().len(),
                current_node.hap_dwfa().query_haplotype().sequence().len()
            );
            debug!("Current queue size: {}", pqueue.len());
            let mut queue_stats: BTreeMap<(usize, usize), usize> = Default::default();
            for (n, _) in pqueue.iter() {
                let k = (n.set_alleles(), n.num_errors());
                *queue_stats.entry(k).or_default() += 1;
            }
            for (k, v) in queue_stats.iter() {
                debug!("\t{k:?} => {v}");
            }
        }

        // now get the relevant variant info depending on truth/query
        let (variant_index, is_truth) = all_variant_order[order_index];
        let (current_variant, current_allele) = if is_truth {
            (&truth_variants[variant_index], truth_alleles[variant_index])
        } else {
            (&query_variants[variant_index], query_alleles[variant_index])
        };

        // get the sync point we can use, which is up to the start of the next variant
        let next_var_pos = if order_index == all_variant_order.len() - 1 {
            // this is the last variant
            coordinates.end() as usize
        } else {
            let (nvi, nt) = all_variant_order[order_index+1];
            if nt { truth_variants[nvi].position() as usize } else { query_variants[nvi].position() as usize }
        };
        let sync_extension = Some(next_var_pos);

        match current_allele {
            Allele::Unknown => bail!("Allele::Unknown is not supported"),
            Allele::Reference => {
                // REF allele only, so just extend by that
                let mut new_node = current_node;

                // no need to update node ID here since we are not cloning; and compatibility from extension is irrelevant
                let success = new_node.extend_variant(
                    reference, is_truth, 
                    current_variant, Allele::Reference, sync_extension, false
                )?;
                
                if success && new_node.is_exact_match() {
                    let new_priority = new_node.priority();
                    pqueue.push(new_node, new_priority);
                } else {
                    *failed_counts.entry(order_index).or_default() += 1;
                }
            },
            Allele::Alternate => {
                let extensions = [
                    (Allele::Reference, true), // expecting ALT, so if we select reference we need to add an error
                    (Allele::Alternate, false) // ALT is not an error
                ];

                // try each extension
                for (allele, is_error) in extensions.into_iter() {
                    if order_index < auto_fail_index && allele != Allele::Reference {
                        // this one need to be Reference only, no matter what; so skip this
                        continue;
                    }

                    // copy our node and update the ID
                    let mut new_node = current_node.clone();
                    new_node.set_node_id(next_node_id);
                    next_node_id += 1;

                    // try to extend it
                    let success = new_node.extend_variant(
                        reference, is_truth,
                        current_variant, allele, sync_extension, is_error
                    )?;

                    // only add 
                    if success && new_node.is_exact_match() {
                        let new_priority = new_node.priority();
                        pqueue.push(new_node, new_priority);
                    } else {
                        *failed_counts.entry(order_index).or_default() += 1;
                    }
                }
            },
        };

        // check if we need to auto-fail an index and clear the queue out a bit
        auto_fail_counts += 1;
        if auto_fail_counts >= auto_fail_threshold {
            let (auto_fail_sub_index, auto_fail_is_truth) = all_variant_order[auto_fail_index];
            let before_size = pqueue.len();
            pqueue = pqueue.into_iter()
                .filter(|(node, _priority)| {
                    let alleles = if auto_fail_is_truth {
                        node.hap_dwfa().truth_haplotype().alleles()
                    } else {
                        node.hap_dwfa().query_haplotype().alleles()
                    };

                    let allele = if auto_fail_sub_index < alleles.len() {
                        alleles[auto_fail_sub_index]
                    } else {
                        Allele::Reference
                    };

                    // only keep the ones that are Reference; i.e. we have "failed" the other
                    allele == Allele::Reference
                }).collect();

            // debug messaging
            let after_size = pqueue.len();
            trace!("Auto-fail triggered for variant {auto_fail_index}: pqueue size filtered from {before_size} to {after_size}.");

            // now set the next node up for auto-fail tracking
            auto_fail_index += 1;
            auto_fail_counts = 0;
        }
    }

    debug!("Fewest errors: {best_error_count}");

    // get the best result
    let best_node = match best_result {
        Some(bn) => bn,
        None => bail!("No result found for problem")
    };

    // bundlé and savé
    let ret = OptimizedAlleles {
        truth_alleles: best_node.hap_dwfa().truth_haplotype().alleles().to_vec(),
        query_alleles: best_node.hap_dwfa().query_haplotype().alleles().to_vec(),
        num_errors: best_node.num_errors()
    };
    Ok(ret)
}

/// Tracks the current haplotypes in our branch-and-bound exploration of the query space
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct ExactMatchNode {
    /// Unique node ID, primarily for deterministic, fixed order
    node_id: u64,
    /// Haplotype comparator, we really only care about those that are non-equal
    hap_dwfa: HaplotypeDWFA,
    /// Tracks the number of errors (flips from ALT to REF on truth or query) so far in our exact matching path
    num_errors: usize
}

/// Nice shortcut for coding; priority is (fewest flips, most set variants, smallest node index)
/// In other words: lowest cost first, then tie-break with first-come, first-serve
type NodePriority = (Reverse<usize>, usize, Reverse<u64>);
impl ExactMatchNode {
    /// Constructor
    /// * `node_id` - the unique node ID
    /// * `region_start` - the start of the region we're solving
    pub fn new(
        node_id: u64, region_start: usize,
    ) -> Self {
        let hap_dwfa = HaplotypeDWFA::new(region_start, 0);
        Self {
            node_id,
            hap_dwfa,
            num_errors: 0
        }
    }

    /// Adds variants sequence to both haplotypes. Returns true if any ALT sequences were successfully incorporated.
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `is_truth` - if true, extends the truth sequences with the variant; otherwise, the query sequences
    /// * `variant` - the variant we are traversing
    /// * `allele` - the REF/ALT indication for the haplotypes
    /// * `sync_extension` - if provided, this will further copy the reference sequence to both truth and query
    pub fn extend_variant(&mut self,
        reference: &[u8], is_truth: bool,
        variant: &Variant, allele: Allele, sync_extension: Option<usize>, is_error: bool
    ) -> Result<bool, HapDWFAError> {
        let extended = match self.hap_dwfa.extend_variant(reference, is_truth, variant, allele, sync_extension) {
            Ok(b) => b,
            Err(e) => {
                if is_allowed_error(&e) {
                    assert!(!self.is_exact_match());
                    false
                } else {
                    return Err(e);
                }
            },
        };
        if is_error {
            self.num_errors += 1;
        }
        Ok(extended)
    }

    /// Finalizes the DWFAs, nothing can get added after calling this.
    /// Returns true if this finalizing was successful and an exact match.
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `region_end` - needed to determine how far our region extends out to
    pub fn finalize_dwfas(&mut self, reference: &[u8], region_end: usize) -> Result<bool, HapDWFAError> {
        // update each DWFA, translate errors to anyhow
        match self.hap_dwfa.finalize_dwfa(reference, region_end) {
            Ok(_) => Ok(true),
            Err(e) => {
                if is_allowed_error(&e) {
                    assert!(!self.is_exact_match());
                    Ok(false)
                } else {
                    Err(e)
                }
            }
        }
    }

    /// Returns the total cost of the current DWFAs
    pub fn edit_distance(&self) -> usize {
        self.hap_dwfa.edit_distance()
    }

    /// Returns true if the two sequences are still capable of an exact match
    pub fn is_exact_match(&self) -> bool {
        self.edit_distance() == 0
    }

    /// Returns true if this node is guaranteed to be on the best path because all future changes will be independent
    pub fn is_synchronized(&self) -> bool {
        self.hap_dwfa.is_synchronized()
    }

    /// Returns the priority for this node
    pub fn priority(&self) -> NodePriority {
        (
            Reverse(self.num_errors()), // fewest number of error
            self.hap_dwfa.set_alleles() - self.num_errors(), // most correctly set alleles
            Reverse(self.node_id) // earliest node ID
        )
    }

    /// Returns the number of set alleles
    pub fn set_alleles(&self) -> usize {
        // self.haplotype1.alleles().len()
        self.hap_dwfa.set_alleles()
    }

    // getters
    pub fn hap_dwfa(&self) -> &HaplotypeDWFA {
        &self.hap_dwfa
    }

    pub fn num_errors(&self) -> usize {
        self.num_errors
    }

    // setters
    pub fn set_node_id(&mut self, node_id: u64) {
        self.node_id = node_id;
    }
}

/// Returns true if this is a "normal" error from hitting the max edit distance
fn is_allowed_error(error: &HapDWFAError) -> bool {
    matches!(
        error,
        // list of errors that are allowed
        HapDWFAError::DError { error: crate::dwfa::dynamic_wfa::DWFAError::MaxEditDistance }
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exact_match_node() {
        // simple test for the diplotype level
        let reference = b"ACGTACGTACGT";
        // let baseline1 = b"ACGTACCGTACGT"; // C insertion
        // let baseline2 = b"ACGTAGGTACGT"; // C->G SNV which overlaps insertion region
        let mut comparison_node = ExactMatchNode::new(0, 0);

        // first, add the insertion
        let ins_variant = Variant::new_insertion(0, 4, b"A".to_vec(), b"AC".to_vec()).unwrap();
        comparison_node.extend_variant(reference, false, &ins_variant, Allele::Alternate, None, false).unwrap();

        // second add the SNP variant
        let snv_variant = Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap();
        comparison_node.extend_variant(reference, false, &snv_variant, Allele::Reference, None, true).unwrap();

        // lets add a homozygous FP also
        let snv_variant2 = Variant::new_snv(0, 8, b"A".to_vec(), b"G".to_vec()).unwrap();
        comparison_node.extend_variant(reference, false, &snv_variant2, Allele::Alternate, None, false).unwrap();

        // finalize and compare
        comparison_node.finalize_dwfas(reference, reference.len()).unwrap();
        // assert_eq!(comparison_node.edit_distance(), 2); // the ED would be 2 if it did not short-circuit to 1
        assert_eq!(comparison_node.edit_distance(), 1);
        assert_eq!(comparison_node.num_errors(), 1);
        assert_eq!(comparison_node.hap_dwfa().truth_haplotype().sequence(), reference);
        assert_eq!(comparison_node.hap_dwfa().truth_haplotype().alleles(), &[]);
        assert_eq!(comparison_node.hap_dwfa().query_haplotype().sequence(), b"ACGTACCGTGCGT");
        assert_eq!(comparison_node.hap_dwfa().query_haplotype().alleles(), &[Allele::Alternate, Allele::Reference, Allele::Alternate]);
    }

    #[test]
    fn test_optimize_gt_alleles_match() {
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("mock".to_string(), 0, reference.len() as u64);
        let truth_variants = vec![
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap()
        ];
        let truth_alleles = vec![
            Allele::Alternate
        ];
        let query_variants = truth_variants.clone();
        let query_alleles = truth_alleles.clone();

        let result = optimize_gt_alleles(
            reference, &coordinates,
            &truth_variants, &truth_alleles,
            &query_variants, &query_alleles
        ).unwrap();

        assert!(result.is_exact_match());
        assert_eq!(result.num_errors(), 0);
        assert_eq!(result.truth_alleles(), &truth_alleles);
        assert_eq!(result.query_alleles(), &query_alleles);
    }

    #[test]
    fn test_optimize_gt_alleles_all_fn() {
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("mock".to_string(), 0, reference.len() as u64);
        let truth_variants = vec![
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap()
        ];
        let truth_alleles = vec![
            Allele::Alternate
        ];
        let query_variants = vec![];
        let query_alleles = vec![];

        let result = optimize_gt_alleles(
            reference, &coordinates,
            &truth_variants, &truth_alleles,
            &query_variants, &query_alleles
        ).unwrap();

        assert!(!result.is_exact_match());
        assert_eq!(result.num_errors(), 1);
        assert_eq!(result.truth_alleles(), &[Allele::Reference]);
        assert_eq!(result.query_alleles(), &query_alleles);
    }

    #[test]
    fn test_optimize_gt_alleles_all_fp() {
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("mock".to_string(), 0, reference.len() as u64);
        let truth_variants = vec![];
        let truth_alleles = vec![];
        let query_variants = vec![
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap()
        ];
        let query_alleles = vec![
            Allele::Alternate
        ];
        
        let result = optimize_gt_alleles(
            reference, &coordinates,
            &truth_variants, &truth_alleles,
            &query_variants, &query_alleles
        ).unwrap();

        assert!(!result.is_exact_match());
        assert_eq!(result.num_errors(), 1);
        assert_eq!(result.truth_alleles(), &truth_alleles);
        assert_eq!(result.query_alleles(), &[Allele::Reference]);
    }

    #[test]
    fn test_optimize_gt_alleles_diff_rep() {
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("mock".to_string(), 0, reference.len() as u64);
        let truth_variants = vec![
            Variant::new_indel(0, 5, b"CGT".to_vec(), b"GGG".to_vec()).unwrap()
        ];
        let truth_alleles = vec![
            Allele::Alternate
        ];
        let query_variants = vec![
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 7, b"T".to_vec(), b"G".to_vec()).unwrap()
        ];
        let query_alleles = vec![
            Allele::Alternate,
            Allele::Alternate
        ];
        
        let result = optimize_gt_alleles(
            reference, &coordinates,
            &truth_variants, &truth_alleles,
            &query_variants, &query_alleles
        ).unwrap();

        assert!(result.is_exact_match());
        assert_eq!(result.num_errors(), 0);
        assert_eq!(result.truth_alleles(), &truth_alleles);
        assert_eq!(result.query_alleles(), &query_alleles);
    }
}
