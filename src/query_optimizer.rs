/*!
# Query Optimizer
Contains the logic for optimizing the query sequences relative to the truth sequences in a haplotype pair.
This is implicitly looking for the best phasing of the query variants to match the truth variants.
Scoring minimizing the edit distance between the two sequences, while respecting the truth phase orientation when provided.
Conceptually, this is a branch-and-bound approach that maximizes the BASEPAIR accuracy downstream.
There may be some edge cases where this is not optimal for GT-based scoring.

## Example usage
```rust
use aardvark_bio::data_types::coordinates::Coordinates;
use aardvark_bio::data_types::phase_enums::PhasedZygosity;
use aardvark_bio::data_types::variants::Variant;
use aardvark_bio::query_optimizer::optimize_sequences;

// create a simple reference string with coordinates for the region
let reference = b"ACGTACGTACGT";
let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);

// create two truth haplotype sequences for comparison later
let truth_seq1 = b"ACGTAAGTACGT"; // C->A SNV
let truth_seq2 = b"ACGTAGGTACGT"; // C->G SNV

// these are the corresponding truth variants and zygosities for the truth sequences
let shared_variants = [
    Variant::new_snv(0, 5, b"C".to_vec(), b"A".to_vec()).unwrap(),
    Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap()
];
let truth_zygosity = [
    PhasedZygosity::PhasedHet10,
    PhasedZygosity::PhasedHet01,
];

// in this example, the variants are identical, but the query zygosities are unphased
let query_zygosity = [
    PhasedZygosity::UnphasedHeterozygous,
    PhasedZygosity::UnphasedHeterozygous,
];

// now run the optimizer
let query_sequences = optimize_sequences(
    reference, &coordinates, &shared_variants, &truth_zygosity, &shared_variants, &query_zygosity, 50
).unwrap()[0].clone();

// should be exact matches, and our phasing is resolved also
assert_eq!(query_sequences.ed1(), 0);
assert_eq!(query_sequences.ed2(), 0);
assert_eq!(query_sequences.query_seq1(), truth_seq1);
assert_eq!(query_sequences.query_seq2(), truth_seq2);
assert_eq!(query_sequences.query_zygosity(), &truth_zygosity);
```
*/

use anyhow::ensure;
use log::debug;
use priority_queue::PriorityQueue;
use std::cmp::Reverse;
use std::hash::Hash;

use crate::data_types::coordinates::Coordinates;
use crate::data_types::phase_enums::{Allele, PhasedZygosity};
use crate::data_types::variants::Variant;
use crate::dwfa::haplotype_dwfa::HaplotypeDWFA;

/// Result from our query optimization routine
#[derive(Clone, Debug)]
pub struct OptimizedHaplotypes {
    /// Final zygosity of truth, aligned with the input variants. `Excluded` was left out.
    truth_zygosity: Vec<PhasedZygosity>,
    /// Truth sequence 1, generally expected to be phased as provided
    truth_seq1: Vec<u8>,
    /// Truth sequence 2
    truth_seq2: Vec<u8>,
    /// Final zygosity of query, aligned with the input variants. `Excluded` was left out.
    query_zygosity: Vec<PhasedZygosity>,
    /// Query sequence 1, this combined with query sequence 2 should be the closest phased match to the truth haplotypes
    query_seq1: Vec<u8>,
    /// Query sequence 2
    query_seq2: Vec<u8>,
    /// Edit distance from query_seq1 to truth_seq1
    ed1: usize,
    /// Edit distance from query_seq2 to truth_seq2
    ed2: usize,
    /// Variant skip distance for truth_seq1, the number of bases that were not incorporated
    truth_vs1: usize,
    /// Variant skip distance for truth_seq2
    truth_vs2: usize,
    /// Variant skip distance for truth_seq1, the number of bases that were not incorporated
    query_vs1: usize,
    /// Variant skip distance for truth_seq2
    query_vs2: usize
}

impl OptimizedHaplotypes {
    /// Returns true if all edit distances are zero AND no variants were skipped in incorporation
    pub fn is_exact_match(&self) -> bool {
        self.ed1 + self.ed2 + self.truth_vs1 + self.truth_vs2 + self.query_vs1 + self.query_vs2 == 0
    }

    // getters
    pub fn truth_zygosity(&self) -> &[PhasedZygosity] {
        &self.truth_zygosity
    }

    pub fn truth_seq1(&self) -> &[u8] {
        &self.truth_seq1
    }

    pub fn truth_seq2(&self) -> &[u8] {
        &self.truth_seq2
    }

    pub fn query_zygosity(&self) -> &[PhasedZygosity] {
        &self.query_zygosity
    }

    pub fn query_seq1(&self) -> &[u8] {
        &self.query_seq1
    }

    pub fn query_seq2(&self) -> &[u8] {
        &self.query_seq2
    }

    pub fn ed1(&self) -> usize {
        self.ed1
    }

    pub fn ed2(&self) -> usize {
        self.ed2
    }

    pub fn truth_vs1(&self) -> usize {
        self.truth_vs1
    }

    pub fn truth_vs2(&self) -> usize {
        self.truth_vs2
    }

    pub fn query_vs1(&self) -> usize {
        self.query_vs1
    }

    pub fn query_vs2(&self) -> usize {
        self.query_vs2
    }
}

/// Given two sets of variants, this will determine the best orientations of both truth and query variants that minimizes the combined edit distance.
/// Truth variant phases are respected, whereas query is allowed to be changed (i.e., treated as an "error").
/// Variants marked optional are applied only when that choice strictly lowers the edit distance.
/// Otherwise the zygosity entry is `PhasedZygosity::Excluded`.
/// The method generates four sequences with two phase zygosities sets.
/// If two heterozygous query variants are incompatible, then they are forced onto separate haplotypes.
/// # Arguments
/// * `reference` - the full reference chromosome from coordinates.chrom
/// * `coordinates` - the coordinates of the region we are comparing
/// * `truth_variants` - the provided set of truth variants
/// * `truth_zygosity` - the provided set of zygosities; assumed to all be in the same phase block
/// * `query_variants` - the set of query variants
/// * `query_zygosity` - the corresponding zygosity of each query variant, must be heterozygous or homozygous alternate
/// * `max_branch_factor` - the maximum branch factor in the query optimizer; limits exponential blowup
/// # Errors
/// * if there are errors from DWFA extension or finalization
/// # Panics
/// * if an unsupported zygosity is provided
pub fn optimize_sequences(
    reference: &[u8], coordinates: &Coordinates,
    truth_variants: &[Variant], truth_zygosity: &[PhasedZygosity],
    query_variants: &[Variant], query_zygosity: &[PhasedZygosity],
    max_branch_factor: usize,
) -> anyhow::Result<Vec<OptimizedHaplotypes>> {
    debug!("Starting sequence optimization...");

    // do some sanity checks
    ensure!(truth_variants.len() == truth_zygosity.len(), "truth values must have equal length");
    ensure!(query_variants.len() == query_zygosity.len(), "query values must have equal length");
    ensure!(max_branch_factor > 0, "max_branch_factor must be greater than 0");

    // first, figure out the order we are calculating variants in
    let all_variant_order: Vec<(usize, bool)> = order_variants(truth_variants, query_variants);
    let total_variant_count = all_variant_order.len();

    // create our initial empty root node
    let mut next_node_id = 0;
    let start = coordinates.start() as usize;
    let root_node = ComparisonNode::new(next_node_id, start);
    next_node_id += 1;

    // put it into the queue as our seed node
    let priority = root_node.priority();
    let mut pqueue: PriorityQueue<ComparisonNode, NodePriority> = Default::default();
    pqueue.push(root_node, priority);

    // set our best to max ED with no result
    let mut best_ed = usize::MAX;
    let mut best_optional_included = usize::MAX;
    let mut best_results = vec![];

    // this tracks how many nodes of length L are encountered, and ignores any beyond our heuristic cutoff
    // in testing, this just prevent exponential blowup from purely FP variants in query
    let mut bucket_counts = vec![0; total_variant_count+1];

    // loop while our queue is not empty
    while let Some((current_node, _priority)) = pqueue.pop() {
        if current_node.total_cost() > best_ed {
            // the cost of this node is worse than something already found, so skip it
            continue;
        }

        // check which variant we are on
        let order_index = current_node.set_alleles();
        if bucket_counts[order_index] == 0 {
            debug!("Best path to {order_index} = {:?}, {} <?> {}, {} <?> {}",
                current_node.priority(),
                current_node.hap_dwfa1().truth_haplotype().sequence().len(),
                current_node.hap_dwfa1().query_haplotype().sequence().len(),
                current_node.hap_dwfa2().truth_haplotype().sequence().len(),
                current_node.hap_dwfa2().query_haplotype().sequence().len()
            );
        }

        // make sure we didn't fill our quota for this bucket size
        if bucket_counts[order_index] >= max_branch_factor {
            continue;
        }
        bucket_counts[order_index] += 1;

        if order_index == total_variant_count {
            // we did every variant, time to finalize the node
            let mut final_node = current_node;
            final_node.finalize_dwfas(reference, coordinates.end() as usize)?;

            // if our cost is less, this is a better solution
            // equal cost keeps the solution that includes fewer optional alleles
            let final_cost = final_node.total_cost();
            let optional_included = included_optional_count(truth_variants, final_node.hap_dwfa1().truth_haplotype().alleles())
                + included_optional_count(query_variants, final_node.hap_dwfa1().query_haplotype().alleles());
            match final_cost.cmp(&best_ed).then(optional_included.cmp(&best_optional_included)) {
                std::cmp::Ordering::Less => {
                    // lower edit distance, or the same distance with fewer optional alleles
                    best_ed = final_cost;
                    best_optional_included = optional_included;
                    best_results = vec![final_node];
                },
                std::cmp::Ordering::Equal => {
                    // equal to previously found result
                    best_results.push(final_node);
                },
                std::cmp::Ordering::Greater => {},
            };
            continue;
        }

        // now get the relevant variant info depending on truth/query
        let (variant_index, is_truth) = all_variant_order[order_index];
        let (current_variant, current_zyg) = if is_truth {
            (&truth_variants[variant_index], truth_zygosity[variant_index])
        } else {
            (&query_variants[variant_index], query_zygosity[variant_index])
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

        // alleles to try for this variant; optional alleles also get an omit branch
        let mut extensions: Vec<(Allele, Allele)> = Vec::with_capacity(3);
        if current_zyg.is_heterozygous() {
            if !is_truth || current_zyg == PhasedZygosity::UnphasedHeterozygous {
                // if it is unphased OR it is a query variant (phase can be wrong), we have to consider both phase orientations
                extensions.push((Allele::Reference, Allele::Alternate)); // 0|1
                extensions.push((Allele::Alternate, Allele::Reference)); // 1|0
            } else {
                // phase is fixed here, figure out which way it is and then just do that one
                let pair = match current_zyg {
                    PhasedZygosity::PhasedHet01 => (Allele::Reference, Allele::Alternate),
                    PhasedZygosity::PhasedHet10 => (Allele::Alternate, Allele::Reference),
                    _ => panic!("should not happen")
                };
                extensions.push(pair);
            }
        } else {
            // if it's not het, the only other valid kind is hom-alt
            assert_eq!(current_zyg, PhasedZygosity::HomozygousAlternate);
            extensions.push((Allele::Alternate, Allele::Alternate));
        }

        // all optional alleles get an exclude branch
        if current_variant.is_optional() {
            // leave the allele out entirely; this is not a reference call
            extensions.push((Allele::Excluded, Allele::Excluded));
        }

        // if there is only one extension, we can just extend in place
        if extensions.len() == 1 {
            let (allele1, allele2) = extensions[0];
            let mut new_node = current_node;
            new_node.extend_variant(
                reference, is_truth,
                current_variant, allele1, allele2, sync_extension
            )?;
            let new_priority = new_node.priority();
            pqueue.push(new_node, new_priority);
        } else {
            // multiple extension, so we need to create new nodes for each
            for (allele1, allele2) in extensions {
                let mut new_node = current_node.clone();
                new_node.set_node_id(next_node_id);
                next_node_id += 1;
                new_node.extend_variant(
                    reference, is_truth,
                    current_variant, allele1, allele2, sync_extension
                )?;
                let new_priority = new_node.priority();
                pqueue.push(new_node, new_priority);
            }
        }
    }
    debug!("Best ED: {best_ed}");

    ensure!(!best_results.is_empty(), "no results found");

    // convert each node into OptimizedHaplotypes
    let ret = best_results.into_iter()
        .map(|best_node| {
            // convert the allelic observations back into phased zygosities for reporting
            let hap_dwfa1 = best_node.hap_dwfa1();
            let hap_dwfa2 = best_node.hap_dwfa2();

            let truth_hap1 = hap_dwfa1.truth_haplotype();
            let truth_hap2 = hap_dwfa2.truth_haplotype();
            let query_hap1 = hap_dwfa1.query_haplotype();
            let query_hap2 = hap_dwfa2.query_haplotype();

            let truth_zygosity = convert_alleles_to_zygosity(truth_hap1.alleles(), truth_hap2.alleles());
            let query_zygosity = convert_alleles_to_zygosity(query_hap1.alleles(), query_hap2.alleles());

            // bundlé and savé
            OptimizedHaplotypes {
                truth_zygosity,
                truth_seq1: truth_hap1.sequence().to_vec(),
                truth_seq2: truth_hap2.sequence().to_vec(),
                query_zygosity,
                query_seq1: query_hap1.sequence().to_vec(),
                query_seq2: query_hap2.sequence().to_vec(),
                ed1: hap_dwfa1.edit_distance(),
                ed2: hap_dwfa2.edit_distance(),
                truth_vs1: hap_dwfa1.truth_haplotype().variant_skip_distance(),
                truth_vs2: hap_dwfa2.truth_haplotype().variant_skip_distance(),
                query_vs1: hap_dwfa1.query_haplotype().variant_skip_distance(),
                query_vs2: hap_dwfa2.query_haplotype().variant_skip_distance(),
            }
        }).collect();
    Ok(ret)
}

/// This will take a collection of truth and query variants and return the relative order.
/// Return value is a Vec of tuple (index, is_truth).
/// # Arguments
/// * `truth_variants` - the pre-sorted set of truth variants
/// * `query_variants` - the pre-sorted set of query variants
pub fn order_variants(truth_variants: &[Variant], query_variants: &[Variant]) -> Vec<(usize, bool)> {
    // add the variants
    let mut ret: Vec<(usize, bool)> = Vec::<(usize, bool)>::with_capacity(truth_variants.len()+query_variants.len());
    ret.extend((0..truth_variants.len()).zip(std::iter::repeat(true)));
    ret.extend((0..query_variants.len()).zip(std::iter::repeat(false)));

    // sort them by positions
    ret.sort_by_key(|&(i, is_truth)| if is_truth { truth_variants[i].position() } else { query_variants[i].position() });
    ret
}

/// Helper function to convert a Vec of allele pairs into PhasedZygosities, all outputs will be phased when possible.
/// TODO: do we want to rewrite this with something like `impl From<(Allele, Allele)> for PhasedZygosity`?
/// # Arguments
/// * `alleles1` - first set of alleles
/// * `alleles2` - second set of alleles
fn convert_alleles_to_zygosity(alleles1: &[Allele], alleles2: &[Allele]) -> Vec<PhasedZygosity> {
    assert_eq!(alleles1.len(), alleles2.len());
    alleles1.iter()
        .zip(alleles2.iter())
        .map(|(a1, a2)| {
            match (a1, a2) {
                // these are the only combinations we should hit at this point in the program
                (Allele::Reference, Allele::Alternate) => PhasedZygosity::PhasedHet01,
                (Allele::Alternate, Allele::Reference) => PhasedZygosity::PhasedHet10,
                (Allele::Alternate, Allele::Alternate) => PhasedZygosity::HomozygousAlternate,
                (Allele::Excluded, Allele::Excluded) => PhasedZygosity::Excluded,
                _ => panic!("no impl")
            }
        })
        .collect()
}

/// Counts optional variants that were applied rather than excluded.
fn included_optional_count(variants: &[Variant], alleles: &[Allele]) -> usize {
    variants.iter()
        .zip(alleles.iter())
        .filter(|(variant, allele)| variant.is_optional() && **allele != Allele::Excluded)
        .count()
}

/// Tracks the current haplotypes in our branch-and-bound exploration of the query space
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct ComparisonNode {
    /// Unique node ID, primarily for deterministic, fixed order
    node_id: u64,
    /// Haplotype 1 information
    hap_dwfa1: HaplotypeDWFA,
    /// Haplotype 2 information
    hap_dwfa2: HaplotypeDWFA
}

/// Nice shortcut for coding; priority is (smallest cost from ED, smallest node index)
/// In other words: lowest cost first, then tie-break with first-come, first-serve
type NodePriority = (Reverse<usize>, Reverse<u64>);

impl ComparisonNode {
    /// Constructor
    /// * `node_id` - the unique node ID
    /// * `region_start` - the start of the region we're solving
    pub fn new(
        node_id: u64, region_start: usize,
    ) -> Self {
        let hap_dwfa1 = HaplotypeDWFA::new(region_start, usize::MAX);
        let hap_dwfa2 = HaplotypeDWFA::new(region_start, usize::MAX);
        Self {
            node_id,
            hap_dwfa1,
            hap_dwfa2
        }
    }

    /// Adds variants sequence to both haplotypes. Returns true if any ALT sequences were successfully incorporated.
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `is_truth` - if true, extends the truth sequences with the variant; otherwise, the query sequences
    /// * `variant` - the variant we are traversing
    /// * `allele1` - the REF/ALT indication for hap1
    /// * `allele2` - the REF/ALT indication for hap2
    /// * `sync_extension` - if provided, this will further copy the reference sequence to both truth and query
    pub fn extend_variant(&mut self,
        reference: &[u8], is_truth: bool,
        variant: &Variant, allele1: Allele, allele2: Allele, sync_extension: Option<usize>
    ) -> anyhow::Result<bool> {
        let mut both_extended = true;
        both_extended &= self.hap_dwfa1.extend_variant(reference, is_truth, variant, allele1, sync_extension)?;
        both_extended &= self.hap_dwfa2.extend_variant(reference, is_truth, variant, allele2, sync_extension)?;
        Ok(both_extended)
    }

    /// Finalizes the DWFAs, nothing can get added after calling this
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `region_end` - needed to determine how far our region extends out to
    pub fn finalize_dwfas(&mut self, reference: &[u8], region_end: usize) -> anyhow::Result<()> {
        // update each DWFA, translate errors to anyhow
        self.hap_dwfa1.finalize_dwfa(reference, region_end)?;
        self.hap_dwfa2.finalize_dwfa(reference, region_end)?;
        Ok(())
    }

    /// Returns the total cost of the current DWFAs
    pub fn total_cost(&self) -> usize {
        self.hap_dwfa1.total_cost() + self.hap_dwfa2.total_cost()
    }

    /// Returns the priority for this node
    pub fn priority(&self) -> NodePriority {
        (
            Reverse(self.total_cost()), // lowest cost first
            Reverse(self.node_id) // then lowest node ID second
        )
    }

    /// Returns the number of set alleles
    pub fn set_alleles(&self) -> usize {
        // self.haplotype1.alleles().len()
        self.hap_dwfa1.set_alleles()
    }

    // getters
    pub fn hap_dwfa1(&self) -> &HaplotypeDWFA {
        &self.hap_dwfa1
    }
    pub fn hap_dwfa2(&self) -> &HaplotypeDWFA {
        &self.hap_dwfa2
    }

    // setters
    pub fn set_node_id(&mut self, node_id: u64) {
        self.node_id = node_id;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MAX_BRANCH_FACTOR_TEST: usize = 50;

    #[test]
    fn test_comparison_node() {
        // simple test for the diplotype level
        let reference = b"ACGTACGTACGT";
        // let baseline1 = b"ACGTACCGTACGT"; // C insertion
        // let baseline2 = b"ACGTAGGTACGT"; // C->G SNV which overlaps insertion region
        let mut comparison_node = ComparisonNode::new(0, 0);

        // first, add the insertion
        let ins_variant = Variant::new_insertion(0, 4, b"A".to_vec(), b"AC".to_vec()).unwrap();
        comparison_node.extend_variant(reference, false, &ins_variant, Allele::Alternate, Allele::Reference, None).unwrap();

        // second add the SNP variant
        let snv_variant = Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap();
        comparison_node.extend_variant(reference, false, &snv_variant, Allele::Reference, Allele::Alternate, None).unwrap();

        // lets add a homozygous FP also
        let snv_variant2 = Variant::new_snv(0, 8, b"A".to_vec(), b"G".to_vec()).unwrap();
        comparison_node.extend_variant(reference, false, &snv_variant2, Allele::Alternate, Allele::Alternate, None).unwrap();

        // finalize and compare
        comparison_node.finalize_dwfas(reference, reference.len()).unwrap();
        assert_eq!(comparison_node.total_cost(), 4); // 1 INS, 1 SNP, 1 hom SNP = 1+1+2
        assert_eq!(comparison_node.hap_dwfa1().query_haplotype().sequence(), b"ACGTACCGTGCGT");
        assert_eq!(comparison_node.hap_dwfa1().query_haplotype().alleles(), &[Allele::Alternate, Allele::Reference, Allele::Alternate]);
        assert_eq!(comparison_node.hap_dwfa2().query_haplotype().sequence(), b"ACGTAGGTGCGT");
        assert_eq!(comparison_node.hap_dwfa2().query_haplotype().alleles(), &[Allele::Reference, Allele::Alternate, Allele::Alternate]);
    }

    #[test]
    fn test_optimize_query_sequences_001() {
        // this test largely mirrors the values from test_comparison_node()
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);

        let truth_variants = [
            Variant::new_insertion(0, 4, b"A".to_vec(), b"AC".to_vec()).unwrap(),
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap(),
        ];
        let truth_zygosity = [
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet01,
        ];

        let query_variants = [
            Variant::new_insertion(0, 4, b"A".to_vec(), b"AC".to_vec()).unwrap(),
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 8, b"A".to_vec(), b"G".to_vec()).unwrap()
        ];
        let query_zygosity = [
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet01,
            PhasedZygosity::HomozygousAlternate
        ];
        
        // first, run a set of empty variants against truth, everything should be FP
        let sequences = optimize_sequences(
            reference, &coordinates, &truth_variants, &truth_zygosity, &query_variants, &query_zygosity, MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();

        assert_eq!(sequences.ed1(), 1);
        assert_eq!(sequences.ed2(), 1);
        assert_eq!(sequences.truth_seq1(), b"ACGTACCGTACGT");
        assert_eq!(sequences.truth_seq2(), b"ACGTAGGTACGT");
        assert_eq!(sequences.truth_zygosity(), &truth_zygosity);
        assert_eq!(sequences.query_seq1(), b"ACGTACCGTGCGT");
        assert_eq!(sequences.query_seq2(), b"ACGTAGGTGCGT");
        assert_eq!(sequences.query_zygosity(), &query_zygosity);
    }

    #[test]
    fn test_optimize_query_sequences_all_fn() {
        // this test has no query variants, so everything is effectively a false negative
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let truth_variants = [
            Variant::new_insertion(0, 4, b"A".to_vec(), b"AC".to_vec()).unwrap(),
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap(),
        ];
        let truth_zygosity = [
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet01,
        ];
        let query_variants = [];
        let query_zygosity = [];

        // first, run a set of empty variants against truth, everything should be FP
        let sequences = optimize_sequences(
            reference, &coordinates, &truth_variants, &truth_zygosity, &query_variants, &query_zygosity, MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();

        assert_eq!(sequences.ed1(), 1);
        assert_eq!(sequences.ed2(), 1);
        assert_eq!(sequences.query_seq1(), reference);
        assert_eq!(sequences.query_seq2(), reference);
        assert_eq!(sequences.query_zygosity(), &[]);
    }

    #[test]
    fn test_optimize_query_sequences_multiallelic() {
        // test two incompatible unphased hets in query, make sure we get 0 ED
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let truth_seq1 = b"ACGTAAGTACGT"; // C->A SNV
        let truth_seq2 = b"ACGTAGGTACGT"; // C->G SNV
        let shared_variants = [
            Variant::new_snv(0, 5, b"C".to_vec(), b"A".to_vec()).unwrap(),
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap()
        ];
        let truth_zygosity = [
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet01,
        ];
        let query_zygosity = [
            PhasedZygosity::UnphasedHeterozygous,
            PhasedZygosity::UnphasedHeterozygous,
        ];

        // first, run a set of empty variants against truth, everything should be FP
        let query_sequences = optimize_sequences(
            reference, &coordinates, &shared_variants, &truth_zygosity, &shared_variants, &query_zygosity, MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();

        // should be exact matches, and our phasing is resolved also
        assert_eq!(query_sequences.ed1(), 0);
        assert_eq!(query_sequences.ed2(), 0);
        assert_eq!(query_sequences.query_seq1(), truth_seq1);
        assert_eq!(query_sequences.query_seq2(), truth_seq2);
        assert_eq!(query_sequences.query_zygosity(), &truth_zygosity);
    }

    #[test]
    fn test_optimize_query_sequences_incompatible() {
        // test two incompatible unphased homs, make sure we get penalized
        let reference = b"ACGTACGTACGT";
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let shared_variants = [
            Variant::new_snv(0, 5, b"C".to_vec(), b"A".to_vec()).unwrap(),
            Variant::new_snv(0, 5, b"C".to_vec(), b"G".to_vec()).unwrap()
        ];
        let truth_zygosity = [
            PhasedZygosity::HomozygousAlternate,
            PhasedZygosity::HomozygousAlternate,
        ];
        let query_zygosity = [
            PhasedZygosity::HomozygousAlternate,
            PhasedZygosity::HomozygousAlternate,
        ];

        // first, run a set of empty variants against truth, everything should be FP
        let query_sequences = optimize_sequences(
            reference, &coordinates, &shared_variants, &truth_zygosity, &shared_variants, &query_zygosity, MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();

        // should be exact matches, and our phasing is resolved also
        assert_eq!(query_sequences.ed1(), 0);
        assert_eq!(query_sequences.ed2(), 0);
        assert_eq!(query_sequences.truth_vs1(), 1); // skipped in both truth and query; 1 * 2 = 2
        assert_eq!(query_sequences.truth_vs2(), 1);
        assert_eq!(query_sequences.query_vs1(), 1);
        assert_eq!(query_sequences.query_vs2(), 1);
        assert_eq!(query_sequences.query_zygosity(), &truth_zygosity); // but the zygosity is preserved
    }

    /// Two short deletions that together match one longer indel.
    /// The long allele is optional and must be kept because it drives edit distance to zero.
    fn split_deletion_fixture() -> (Vec<u8>, Vec<Variant>, Variant) {
        let mut reference = Vec::new();
        reference.extend_from_slice(b"AAAA");
        reference.extend(std::iter::repeat_n(b'C', 25));
        reference.push(b'T');
        reference.extend(std::iter::repeat_n(b'G', 26));
        reference.extend_from_slice(b"AAAA");

        let short1 = Variant::new_deletion(0, 3, reference[3..29].to_vec(), reference[3..4].to_vec()).unwrap();
        let short2 = Variant::new_deletion(0, 29, reference[29..56].to_vec(), reference[29..30].to_vec()).unwrap();
        let mut long = Variant::new_indel(0, 3, reference[3..56].to_vec(), b"AT".to_vec()).unwrap();
        long.set_optional();
        assert_eq!(short1.alt_ed().unwrap(), 25);
        assert_eq!(short2.alt_ed().unwrap(), 26);
        assert_eq!(long.alt_ed().unwrap(), 51);
        (reference, vec![short1, short2], long)
    }

    /// An optional long query allele that reproduces two required truth deletions is kept and drives edit distance to zero.
    #[test]
    fn test_optional_long_query_matches_two_short_truth() {
        let (reference, shorts, long) = split_deletion_fixture();
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let zyg = vec![PhasedZygosity::HomozygousAlternate, PhasedZygosity::HomozygousAlternate];
        let result = optimize_sequences(
            &reference, &coordinates, &shorts, &zyg, &[long], &[PhasedZygosity::HomozygousAlternate], MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();
        assert_eq!(result.ed1() + result.ed2(), 0);
        assert_eq!(result.query_zygosity(), &[PhasedZygosity::HomozygousAlternate]);
        assert_eq!(result.truth_zygosity(), &zyg);
        assert_eq!(result.query_seq1(), result.truth_seq1());
    }

    /// An optional long truth allele that reproduces two required query deletions is kept and drives edit distance to zero.
    #[test]
    fn test_optional_long_truth_matches_two_short_query() {
        let (reference, shorts, long) = split_deletion_fixture();
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let short_zyg = vec![PhasedZygosity::HomozygousAlternate, PhasedZygosity::HomozygousAlternate];
        let result = optimize_sequences(
            &reference, &coordinates, &[long], &[PhasedZygosity::HomozygousAlternate], &shorts, &short_zyg, MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();
        assert_eq!(result.ed1() + result.ed2(), 0);
        assert_eq!(result.truth_zygosity(), &[PhasedZygosity::HomozygousAlternate]);
        assert_eq!(result.query_zygosity(), &short_zyg);
        assert_eq!(result.query_seq1(), result.truth_seq1());
    }

    /// An optional query allele with no truth counterpart is left out, so the edit distance stays zero.
    #[test]
    fn test_optional_unmatched_query_is_omitted() {
        let reference = b"ACGTACGTACGTACGTACGT";
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let mut alt = b"A".to_vec();
        alt.extend(std::iter::repeat_n(b'N', 80));
        let mut insertion = Variant::new_insertion(0, 0, b"A".to_vec(), alt).unwrap();
        insertion.set_optional();
        let result = optimize_sequences(
            reference, &coordinates, &[], &[], &[insertion], &[PhasedZygosity::HomozygousAlternate], MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();
        assert_eq!(result.ed1() + result.ed2(), 0);
        assert_eq!(result.query_zygosity(), &[PhasedZygosity::Excluded]);
        assert_eq!(result.query_seq1(), reference);
    }

    /// A required query allele with no truth counterpart stays in the comparison and adds edit distance.
    #[test]
    fn test_required_unmatched_query_is_kept() {
        let reference = b"ACGTACGTACGTACGTACGT";
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let mut alt = b"C".to_vec();
        alt.extend(std::iter::repeat_n(b'N', 10));
        let insertion = Variant::new_insertion(0, 1, b"C".to_vec(), alt).unwrap();
        let result = optimize_sequences(
            reference, &coordinates, &[], &[], &[insertion], &[PhasedZygosity::HomozygousAlternate], MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();
        assert!(result.ed1() + result.ed2() > 0);
        assert_eq!(result.query_zygosity(), &[PhasedZygosity::HomozygousAlternate]);
    }

    /// An optional query allele that only partially covers a required truth allele is left out.
    #[test]
    fn test_optional_partial_overlap_is_omitted() {
        let (reference, shorts, _) = split_deletion_fixture();
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let truth = vec![shorts[0].clone()];
        let truth_zyg = vec![PhasedZygosity::HomozygousAlternate];
        let mut partial = Variant::new_deletion(0, 3, reference[3..55].to_vec(), reference[3..4].to_vec()).unwrap();
        partial.set_optional();
        assert_eq!(partial.alt_ed().unwrap(), 51);

        let with_optional = optimize_sequences(
            &reference, &coordinates, &truth, &truth_zyg, &[partial], &[PhasedZygosity::HomozygousAlternate], MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();
        let truth_only = optimize_sequences(
            &reference, &coordinates, &truth, &truth_zyg, &[], &[], MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();
        assert_eq!(with_optional.query_zygosity(), &[PhasedZygosity::Excluded]);
        assert_eq!(with_optional.ed1(), truth_only.ed1());
        assert_eq!(with_optional.ed2(), truth_only.ed2());
    }

    /// An optional query allele that fixes two truth deletions but leaves a separate truth SNV is kept, and the edit distance falls without reaching zero.
    #[test]
    fn test_optional_improvement_is_not_exact() {
        let (reference, shorts, long) = split_deletion_fixture();
        let coordinates = Coordinates::new("chrom".to_string(), 0, reference.len() as u64);
        let snv = Variant::new_snv(0, 57, b"A".to_vec(), b"G".to_vec()).unwrap();
        let mut truth = shorts;
        truth.push(snv);
        let truth_zyg = vec![PhasedZygosity::HomozygousAlternate; 3];

        let with_optional = optimize_sequences(
            &reference, &coordinates, &truth, &truth_zyg, &[long], &[PhasedZygosity::HomozygousAlternate], MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();
        let truth_only = optimize_sequences(
            &reference, &coordinates, &truth, &truth_zyg, &[], &[], MAX_BRANCH_FACTOR_TEST
        ).unwrap()[0].clone();

        assert_eq!(with_optional.query_zygosity(), &[PhasedZygosity::HomozygousAlternate]);
        assert!(with_optional.ed1() + with_optional.ed2() > 0);
        assert!(with_optional.ed1() + with_optional.ed2() < truth_only.ed1() + truth_only.ed2());
    }
}