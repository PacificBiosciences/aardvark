/*!
# Merge Solver
Contains the logic for merging variant calls from multiple inputs.
In general, the main use case is identifying regions that are identical across all inputs.
However, there are controls via the `MergeConfig` struct that allow for more flexible merging strategies when there are conflicts.

## Example usage
```rust
use aardvark_bio::data_types::coordinates::Coordinates;
use aardvark_bio::data_types::merge_benchmark::MergeClassification;
use aardvark_bio::data_types::multi_region::MultiRegion;
use aardvark_bio::data_types::phase_enums::PhasedZygosity;
use aardvark_bio::data_types::variants::Variant;
use aardvark_bio::merge_solver::{solve_merge_region, MergeConfig, MergeConfigBuilder};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;

// create a simple reference genome
let mut reference_genome = ReferenceGenome::empty_reference();
reference_genome.add_contig(
    "mock_chr1".to_string(), "ACCGTTACCAGGACTTGACAAACCG"
).unwrap();
let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 25);

// create a 3-way merge of identical variants with different zygosity
let variants = vec![
    vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
    vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
    vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
];
let zygosity = vec![
    vec![PhasedZygosity::HomozygousAlternate],
    vec![PhasedZygosity::PhasedHet01],
    vec![PhasedZygosity::HomozygousAlternate],
];

// put it all into a MultiRegion object and solve it; this one will fail because of the zygosity difference
let problem = MultiRegion::new(0, coordinates.clone(), variants.clone(), zygosity).unwrap();
let result = solve_merge_region(&problem, &reference_genome, MergeConfig::default()).unwrap();
assert_eq!(result.merge_classification(), &MergeClassification::Different);

// now lets re-try with a voting strategy; here we get the majority agree, and the ones that agree are indices 0 and 2
let merge_config = MergeConfigBuilder::default()
    .majority_voting_enabled(true)
    .build().unwrap();
let result = solve_merge_region(&problem, &reference_genome, merge_config).unwrap();
assert_eq!(result.merge_classification(), &MergeClassification::MajorityAgree { indices: vec![0, 2] });
```
*/
use anyhow::ensure;
use derive_builder::Builder;
use log::debug;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;
use std::collections::BTreeSet;

use crate::data_types::merge_benchmark::{MergeBenchmark, MergeClassification};
use crate::data_types::multi_region::MultiRegion;
use crate::data_types::phase_enums::PhasedZygosity;
use crate::data_types::variants::Variant;
use crate::query_optimizer::optimize_sequences;

/// Controls what passes our merging process
#[derive(Builder, Clone, Copy, Default)]
#[builder(default)]
pub struct MergeConfig {
    /// if true, then blocks where every either agrees or has no variants will pass
    no_conflict_enabled: bool,
    /// if true, a majority vote (all equal) will pass
    majority_voting_enabled: bool,
    /// If Some(v_index), then the corresponding VCF index will be selected if other method fails. This effectively removes failed regions.
    conflict_selection: Option<usize>
}

impl MergeConfig {
    // mostly getters
    pub fn no_conflict_enabled(&self) -> bool {
        self.no_conflict_enabled
    }

    pub fn majority_voting_enabled(&self) -> bool {
        self.majority_voting_enabled
    }

    pub fn conflict_selection(&self) -> Option<usize> {
        self.conflict_selection
    }
}

/// Entry point for comparing a region
/// # Arguments
/// * `problem` - core problem that we want to do the comparison on
/// * `reference_genome` - shared pre-loaded reference genome, intended to be provided from Arc<ReferenceGenome> reference in parallelization
/// * `max_edit_distance` - maximum edit distance between query and truth in the WFA comparison
pub fn solve_merge_region(problem: &MultiRegion, reference_genome: &ReferenceGenome, merge_strategy: MergeConfig) -> anyhow::Result<MergeBenchmark> {
    // pull out core components from the problem space
    let problem_id = problem.region_id();
    let coordinates = problem.coordinates();
    debug!("B#{problem_id} Scanning region #{problem_id}: {coordinates}");

    let reference = reference_genome.get_full_chromosome(coordinates.chrom());

    // first, build up the delta values as a shortcut
    let mut delta_scores = vec![];
    for (variants, zygs) in problem.variants().iter().zip(problem.zygosity().iter()) {
        delta_scores.push(variant_delta_length(variants, zygs)?);
        // delta_scores.push(0);
    }

    // do all the pair-wise comparisons
    let mut all_identical = true;
    let mut no_conflict = true;
    let mut match_sets: Vec<BTreeSet<usize>> = vec![Default::default(); problem.variants().len()];
    for (i, (v1, z1)) in problem.variants().iter().zip(problem.zygosity().iter()).enumerate() {
        match_sets[i].insert(i);
        for j in (i+1)..problem.variants().len() {
            let v2 = &problem.variants()[j];
            let z2 = &problem.zygosity()[j];

            let exact_match = if delta_scores[i] == delta_scores[j] {
                // check the basepair level analysis
                let all_opt_haps = optimize_sequences(
                    reference, coordinates, v1, z1, v2, z2
                )?;

                // we don't care which one is selected, we just want to check if the optimal score is ED=0
                let basepair_analysis = &all_opt_haps[0];
                basepair_analysis.is_exact_match()
            } else {
                // they have different number of inserted/deleted bases, impossible to be an exact match
                false
            };

            // handles the match category
            all_identical &= exact_match;

            // IF either variant set is empty
            no_conflict &= v1.is_empty() || v2.is_empty() ||
                // OR they exactly match, THEN there is no conflict present
                exact_match;

            if exact_match {
                match_sets[i].insert(j);
                match_sets[j].insert(i);
            }
        }
    }

    // figure out if we have a majority
    let maj_count = match_sets.len() / 2 + 1;
    let first_maj = match_sets.into_iter()
        .find(|bset| bset.len() >= maj_count)
        .map(|bset| bset.into_iter().collect::<Vec<usize>>())
        .unwrap_or_default(); // will be empty if we do not

    // figure out the summary classification
    let classification = if all_identical {
        MergeClassification::BasepairIdentical
    } else if merge_strategy.no_conflict_enabled() && no_conflict {
        // collect the non-conflicting indices and return them
        let indices = problem.variants().iter().enumerate()
            .filter_map(|(i, v)| {
                if v.is_empty() {
                    None
                } else {
                    Some(i)
                }
            })
            .collect();
        MergeClassification::NoConflict {
            indices
        }
    } else if merge_strategy.majority_voting_enabled() && !first_maj.is_empty() {
        MergeClassification::MajorityAgree { indices: first_maj }
    } else if let Some(v_index) = merge_strategy.conflict_selection() {
        // all previous have failed and we would normally return Different, but the user has said to report a particular VCF
        MergeClassification::ConflictSelection { index: v_index }
    } else {
        MergeClassification::Different
    };
    
    Ok(MergeBenchmark::new(problem_id, classification))
}

/// Calculates the total size delta for a collection of variants and zygosities.
/// For example, a 1 bp homozygous insertion would be +2, while a heterozygous 3 bp deletion would be -3.
/// SNPs have no impact on length.
/// # Arguments
/// * `variants` - the set of variants
/// * `zygosities` - the corresponding set of zygosities
/// # Errors
/// * if the number of variants and zygosities are not equal
/// * if an unknown zygosity is provided
pub fn variant_delta_length(variants: &[Variant], zygosities: &[PhasedZygosity]) -> anyhow::Result<i64> {
    let mut total_delta = 0;
    ensure!(variants.len() == zygosities.len(), "Variants and zygosities must be the same length");

    for (v, z) in variants.iter().zip(zygosities.iter()) {
        ensure!(*z != PhasedZygosity::Unknown, "Zygosity cannot be unknown");
        let v_delta = v.allele1().len() as i64 - v.allele0().len() as i64;
        let z_weight = z.to_allele_count() as i64;
        total_delta += v_delta * z_weight;
    }

    Ok(total_delta)
}

#[cfg(test)]
mod tests {
    use crate::data_types::coordinates::Coordinates;
    use crate::data_types::phase_enums::PhasedZygosity;
    use crate::data_types::variants::Variant;

    use super::*;

    /// Helper function that builds a tiny reference genome we can repeatedly use
    fn generate_simple_reference() -> ReferenceGenome {
        let mut ref_genome = ReferenceGenome::empty_reference();
        ref_genome.add_contig(
            "mock_chr1".to_string(), "ACCGTTACCAGGACTTGACAAACCG"
        ).unwrap();
        ref_genome
    }

    #[test]
    fn test_exact_mode() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 25);
        let variants = vec![
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
        ];
        let zygosity = vec![
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::HomozygousAlternate],
        ];
        let problem = MultiRegion::new(0, coordinates.clone(), variants.clone(), zygosity).unwrap();

        let result = solve_merge_region(&problem, &reference_genome, MergeConfig::default()).unwrap();
        assert_eq!(result.merge_classification(), &MergeClassification::BasepairIdentical);

        // different zygosity should get different
        let zygosity = vec![
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::PhasedHet01]
        ];
        let problem = MultiRegion::new(0, coordinates, variants, zygosity).unwrap();

        let result = solve_merge_region(&problem, &reference_genome, MergeConfig::default()).unwrap();
        assert_eq!(result.merge_classification(), &MergeClassification::Different);
    }

    #[test]
    fn test_noconflict_mode() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 25);
        let variants = vec![
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![]
        ];
        let zygosity = vec![
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::HomozygousAlternate],
            vec![],
        ];
        let problem = MultiRegion::new(0, coordinates.clone(), variants.clone(), zygosity).unwrap();
        let merge_config = MergeConfigBuilder::default().no_conflict_enabled(true).build().unwrap();

        let result = solve_merge_region(&problem, &reference_genome, merge_config).unwrap();
        assert_eq!(result.merge_classification(), &MergeClassification::NoConflict { indices: vec![0, 1] });

        // make sure this still gets different
        let variants = vec![
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()]
        ];
        let zygosity = vec![
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::PhasedHet01]
        ];
        let problem = MultiRegion::new(0, coordinates, variants, zygosity).unwrap();

        let result = solve_merge_region(&problem, &reference_genome, merge_config).unwrap();
        assert_eq!(result.merge_classification(), &MergeClassification::Different);
    }

    #[test]
    fn test_majority_mode() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 25);
        let variants = vec![
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![]
        ];
        let zygosity = vec![
            vec![PhasedZygosity::HomozygousAlternate],
            vec![PhasedZygosity::HomozygousAlternate],
            vec![],
        ];
        let problem = MultiRegion::new(0, coordinates.clone(), variants.clone(), zygosity).unwrap();
        let merge_config = MergeConfigBuilder::default().majority_voting_enabled(true).build().unwrap();

        let result = solve_merge_region(&problem, &reference_genome, merge_config).unwrap();
        assert_eq!(result.merge_classification(), &MergeClassification::MajorityAgree { indices: vec![0, 1] });

        // make sure this still gets different
        let variants = vec![
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()],
            vec![],
            vec![Variant::new_snv(0, 5, b"T".to_vec(), b"C".to_vec()).unwrap()]
        ];
        let zygosity = vec![
            vec![PhasedZygosity::HomozygousAlternate],
            vec![],
            vec![PhasedZygosity::PhasedHet01]
        ];
        let problem = MultiRegion::new(0, coordinates, variants, zygosity).unwrap();

        let result = solve_merge_region(&problem, &reference_genome, merge_config).unwrap();
        assert_eq!(result.merge_classification(), &MergeClassification::Different);
    }

    #[test]
    fn test_variant_length_delta() {
        let variants = [
            Variant::new_snv(0, 10, b"A".to_vec(), b"C".to_vec()).unwrap(),
            Variant::new_deletion(0, 12, b"ACGTACGT".to_vec(), b"A".to_vec()).unwrap(),
            Variant::new_insertion(0, 25, b"A".to_vec(), b"ACC".to_vec()).unwrap()
        ];
        let zygs = [
            PhasedZygosity::HomozygousAlternate,
            PhasedZygosity::PhasedHet01,
            PhasedZygosity::HomozygousAlternate,
        ];

        // test them individually
        assert_eq!(variant_delta_length(&variants[0..1], &zygs[0..1]).unwrap(), 0);
        assert_eq!(variant_delta_length(&variants[1..2], &zygs[1..2]).unwrap(), -7);
        assert_eq!(variant_delta_length(&variants[2..3], &zygs[2..3]).unwrap(), 4);

        // test them combined
        assert_eq!(variant_delta_length(&variants, &zygs).unwrap(), -3);
    }
}
