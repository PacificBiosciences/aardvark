
use anyhow::{bail, ensure, Context};
use derive_builder::Builder;
use log::debug;
use rust_lib_reference_genome::reference_genome::ReferenceGenome;

use crate::data_types::compare_benchmark::{CompareBenchmark, SequenceBundle};
use crate::data_types::compare_region::CompareRegion;
use crate::data_types::coordinates::Coordinates;
use crate::data_types::phase_enums::{Allele, Haplotype, PhasedZygosity};
use crate::data_types::summary_metrics::SummaryMetrics;
use crate::data_types::variants::{Variant, VariantType};
use crate::exact_gt_optimizer::optimize_gt_alleles;
use crate::parsing::stratifications::Stratifications;
use crate::query_optimizer::{optimize_sequences, OptimizedHaplotypes};
use crate::util::sequence_alignment::wfa_ed;

/// These are the variants we currently allow
const SUPPORTED_VARIANT_TYPES: [VariantType; 4] = [
    VariantType::Snv,
    VariantType::Insertion,
    VariantType::Deletion,
    VariantType::Indel
];

/// Controls what passes our merging process
#[derive(Builder, Clone, Copy, Debug)]
#[builder(default)]
pub struct CompareConfig {
    /// if True, saves the haplotype sequences which may consume more memory than normal
    enable_sequences: bool,
    /// if True, this will enable a compute shortcut for exact-matching sequences at the cost of some variant-level assessment accuracy
    enable_exact_shortcut: bool,
}

impl Default for CompareConfig {
    fn default() -> Self {
        // these settings are set to reasonable defaults for unit tests
        // main.rs will set each of them manually based on user input
        Self {
            enable_sequences: true,
            enable_exact_shortcut: false
        }
    }
}

/// Entry point for comparing a region
/// # Arguments
/// * `problem` - core problem that we want to do the comparison on
/// * `reference_genome` - shared pre-loaded reference genome, intended to be provided from Arc<ReferenceGenome> reference in parallelization
/// * `compare_config` - collection of configuration items for the compare mode
pub fn solve_compare_region(
    problem: &CompareRegion, reference_genome: &ReferenceGenome, compare_config: CompareConfig, stratifications: Option<&Stratifications>
) -> anyhow::Result<CompareBenchmark> {
    // pull out core components from the problem space
    let problem_id = problem.region_id();
    let coordinates = problem.coordinates();
    debug!("B#{problem_id} Scanning region #{problem_id}: {coordinates}");

    // pull out the relevant problem info
    let reference = reference_genome.get_full_chromosome(coordinates.chrom()); // will panic if chrom is not present
    let truth_variants = problem.truth_variants();
    let raw_truth_zygosity = problem.truth_zygosity();
    let query_variants = problem.query_variants();
    let raw_query_zygosity = problem.query_zygosity();

    // pull out the reference sequence also
    let ref_start = coordinates.start() as usize;
    let ref_end = coordinates.end() as usize;
    let ref_seq = &reference[ref_start..ref_end];

    // first, generate the best query sequences based on the truth sequences
    let all_optimized_haplotypes = optimize_sequences(
        reference, coordinates,
        truth_variants, raw_truth_zygosity,
        query_variants, raw_query_zygosity
    )?;
    debug!("B#{problem_id} Found {} equal solutions, scoring by ALT flips...", all_optimized_haplotypes.len());

    let containment_regions = if let Some(strat) = stratifications {
        // get variant coordinates and convert into 0-based inclusive
        let var_coor = problem.var_coordinates();
        let first: i32 = var_coor.start().try_into()
            .with_context(|| format!("Failed to convert coordinates: {var_coor:?}"))?;
        let last: i32 = var_coor.end().try_into()
            .with_context(|| format!("Failed to convert coordinates: {var_coor:?}"))?;
        assert!(first < last);

        let containment_regions = strat.containments(
            var_coor.chrom(), first, last - 1
        );
        Some(containment_regions)
    } else {
        None
    };

    let mut best_results = vec![];
    for optimized_haplotypes in all_optimized_haplotypes.into_iter() {
        // check if the sequences are identical with 0 errors and 0 skipped variants
        if compare_config.enable_exact_shortcut && optimized_haplotypes.is_exact_match() {
            // TODO: account for phasing if we ever add it
            debug!("B#{problem_id} exact match identified:");

            // this compute truth and query at once, under the assumption of exact match
            let mut exact_result = generate_exact_match(
                problem,
                truth_variants, raw_truth_zygosity,
                query_variants, raw_query_zygosity,
                &reference[ref_start..ref_end], &optimized_haplotypes
            )?;

            if compare_config.enable_sequences {
                let bundle = SequenceBundle::new(
                    String::from_utf8(ref_seq.to_vec())?,
                    String::from_utf8(optimized_haplotypes.truth_seq1().to_vec())?,
                    String::from_utf8(optimized_haplotypes.truth_seq2().to_vec())?,
                    String::from_utf8(optimized_haplotypes.query_seq1().to_vec())?,
                    String::from_utf8(optimized_haplotypes.query_seq2().to_vec())?,
                );
                exact_result.add_sequence_bundle(bundle);
            }

            if let Some(cr) = containment_regions {
                exact_result.add_containment_regions(cr);
            }

            return Ok(exact_result);
        }

        // if we made it here, we did not find an exact match, so we need to tease apart which variants are FP and FN
        let truth_zyg = optimized_haplotypes.truth_zygosity();
        let query_zyg = optimized_haplotypes.query_zygosity();

        // split the zygosities into alleles
        let (truth_hap1, truth_hap2): (Vec<Allele>, Vec<Allele>) = truth_zyg.iter()
            .map(|z| z.decompose_alleles())
            .unzip();
        let (query_hap1, query_hap2): (Vec<Allele>, Vec<Allele>) = query_zyg.iter()
            .map(|z| z.decompose_alleles())
            .unzip();

        // now compare the truth and query haplotypes
        let hap1_compare = optimize_gt_alleles(
            reference, coordinates,
            truth_variants, &truth_hap1,
            query_variants, &query_hap1
        )?;
        let hap2_compare = optimize_gt_alleles(
            reference, coordinates,
            truth_variants, &truth_hap2,
            query_variants, &query_hap2
        )?;

        // first iterate over all truth inputs
        let mut truth_stats  = compare_expected_observed(
            problem_id,
            optimized_haplotypes.ed1(),
            optimized_haplotypes.ed2(),
            truth_variants,
            &truth_hap1,
            hap1_compare.truth_alleles(),
            &truth_hap2,
            hap2_compare.truth_alleles()
        )?;

        if compare_config.enable_sequences {
            let bundle = SequenceBundle::new(
                String::from_utf8(ref_seq.to_vec())?,
                String::from_utf8(optimized_haplotypes.truth_seq1().to_vec())?,
                String::from_utf8(optimized_haplotypes.truth_seq2().to_vec())?,
                String::from_utf8(optimized_haplotypes.query_seq1().to_vec())?,
                String::from_utf8(optimized_haplotypes.query_seq2().to_vec())?,
            );
            truth_stats.add_sequence_bundle(bundle);
        }

        // now iterate over all the query inputs
        let query_stats = compare_expected_observed(
            problem_id,
            optimized_haplotypes.ed1(),
            optimized_haplotypes.ed2(),
            query_variants,
            &query_hap1,
            hap1_compare.query_alleles(),
            &query_hap2,
            hap2_compare.query_alleles()
        )?;

        best_results.push((optimized_haplotypes, hap1_compare, hap2_compare, truth_stats, query_stats));
    }

    // find the one with the least error by the number of ALTs flipped to REF; should all have same ED
    let (optimized_haplotypes, h1c, h2c, mut truth_stats, query_stats) = best_results.into_iter()
        .min_by_key(|(_oh, h1c, h2c, _ts, _qs)| h1c.num_errors() + h2c.num_errors()).unwrap();
    debug!("Lowest cost by ALT flips: {}", h1c.num_errors()+h2c.num_errors());

    // add query as the swap stats
    truth_stats.add_swap_benchmark(&query_stats)?;

    // now add all the basepair-level stats
    add_basepair_stats(problem, reference_genome, &mut truth_stats, &optimized_haplotypes)?;

    // add in the containment regions if we have them
    if let Some(cr) = containment_regions {
        truth_stats.add_containment_regions(cr);
    }

    // add in the
    Ok(truth_stats)
}

/// This compares the expected and observed alleles and generated the statistics in our summary output.
/// # Arguments
/// * `problem_id` - the region ID, passed into the result
/// * `ed1` - edit distance for expected and observed hap1, passed into the result
/// * `ed2` - edit distance for expected and observed hap2, passed into the result
/// * `exp_hap1` - expected alleles in hap1
/// * `obs_hap1` - observed alleles in hap1
/// * `exp_hap2` - expected alleles in hap2
/// * `obs_hap2` - observed alleles in hap2
#[allow(clippy::too_many_arguments)]
fn compare_expected_observed(
    problem_id: u64, ed1: usize, ed2: usize,
    variants: &[Variant],
    exp_hap1: &[Allele], obs_hap1: &[Allele],
    exp_hap2: &[Allele], obs_hap2: &[Allele],
) -> anyhow::Result<CompareBenchmark> {
    // sanity checks
    ensure!(variants.len() == exp_hap1.len(), "all variants and haplotype alleles must be equal length");
    ensure!(variants.len() == obs_hap1.len(), "all variants and haplotype alleles must be equal length");
    ensure!(variants.len() == exp_hap2.len(), "all variants and haplotype alleles must be equal length");
    ensure!(variants.len() == obs_hap2.len(), "all variants and haplotype alleles must be equal length");

    // initialize our return data
    let mut benchmark = CompareBenchmark::new(
        problem_id, ed1, ed2
    );

    // zip and iterate
    let zip1 = exp_hap1.iter().cloned().zip(obs_hap1.iter().cloned());
    let zip2 = exp_hap2.iter().cloned().zip(obs_hap2.iter().cloned());
    for (v, ((exp_h1, obs_h1), (exp_h2, obs_h2))) in variants.iter().zip(zip1.zip(zip2)) {
        // convert alleles into allele counts
        let exp_count = exp_h1.to_allele_count() + exp_h2.to_allele_count();
        let obs_count = obs_h1.to_allele_count() + obs_h2.to_allele_count();

        // with the new setup, we shouldn't ever find a false positive
        assert!(exp_count >= obs_count);
        benchmark.add_truth_zygosity(v, exp_count, obs_count)?;
    }

    Ok(benchmark)
}

/// Adds the basepair level statistics to our final result
/// # Arguments
/// * `problem` - the problem we're looking at
/// * `reference_genome` - full reference we use
/// * `bench_result` - a mutable result that we are adding basepair statistics to
/// * `optimized_haps` - the optimized truth/query paired haplotypes
fn add_basepair_stats(
    problem: &CompareRegion, reference_genome: &ReferenceGenome,
    bench_result: &mut CompareBenchmark, optimized_haps: &OptimizedHaplotypes
) -> anyhow::Result<()> {
    // pull out problem fields
    let problem_id = problem.region_id();
    let coordinates = problem.coordinates();
    let reference = reference_genome.get_full_chromosome(coordinates.chrom()); // will panic if chrom is not present
    let ref_start = coordinates.start() as usize;
    let ref_end = coordinates.end() as usize;

    // get our variants and zygosity from the problem and optimized solution
    let truth_variants = problem.truth_variants();
    let truth_zygosity = optimized_haps.truth_zygosity();
    let query_variants = problem.query_variants();
    let query_zygosity = optimized_haps.query_zygosity();

    for haplotype in [Haplotype::Hap1, Haplotype::Hap2] {
        // re-generate the sequences using the haplotypes
        let (truth_seq, truth_ed) = generate_haplotype_sequence(
            reference_genome, coordinates,
            truth_variants, truth_zygosity, haplotype
        )?;
        let (query_seq, query_ed) = generate_haplotype_sequence(
            reference_genome, coordinates,
            query_variants, query_zygosity, haplotype
        )?;

        // TODO: remove these assertions once we are happy with it
        let p_truth_seq = if haplotype == Haplotype::Hap1 { optimized_haps.truth_seq1() } else { optimized_haps.truth_seq2() };
        assert_eq!(&truth_seq, p_truth_seq);
        let p_query_seq = if haplotype == Haplotype::Hap1 { optimized_haps.query_seq1() } else { optimized_haps.query_seq2() };
        assert_eq!(&query_seq, p_query_seq);

        // now do a basepair level comparison
        let align_metrics = perform_basepair_compare(
            &reference[ref_start..ref_end],
            &truth_seq, &query_seq
        )?;
        debug!("B#{problem_id} Align metrics: {align_metrics:?}");
        bench_result.add_basepair_metrics(align_metrics, None);

        // also add negatives from skipping variant incorporation; remember these are doubled
        let truth_vs = 2 * truth_ed as u64;
        let query_vs = 2 * query_ed as u64;
        let skip_metrics = SummaryMetrics::new(0, truth_vs, 0, query_vs);
        bench_result.add_basepair_metrics(skip_metrics, None);

        // for each variant type, we need to filter the right variants in, and then re-compute
        for &filter_type in SUPPORTED_VARIANT_TYPES.iter() {
            // first, we filter the query sequence by variant type to get TP & FP; FNs have to be ignored since we filtered query, TP may be different
            let (filt_qv, filt_qz): (Vec<Variant>, Vec<PhasedZygosity>) = query_variants.iter().zip(query_zygosity.iter())
                .filter_map(|(qv, &qz)| if qv.variant_type() == filter_type {
                    // put the clones here to reduce overhead
                    Some((qv.clone(), qz))
                } else {
                    None
                })
                .unzip();

            let (query_tp, query_fp) = if !filt_qv.is_empty() {
                // we have at least one variant left over in the query, so we need to add these stats
                let (filtered_query, failed_ed) = generate_haplotype_sequence(
                    reference_genome, coordinates, &filt_qv, &filt_qz, haplotype
                ).with_context(|| format!("Error while generating sequence for filtered query haplotype 1: {filter_type:?}"))?;
                let failed_vs = 2 * failed_ed as u64;

                let filt_align_metrics = perform_basepair_compare(
                    &reference[ref_start..ref_end],
                    &truth_seq, &filtered_query
                )?;

                (filt_align_metrics.query_tp, filt_align_metrics.query_fp + failed_vs) // add the failed incorporations
            } else {
                (0, 0)
            };

            // second, we filter the truth sequence by variant type to get TP & FN; FPs have to be ignored since we filtered query, TP may be different
            let (filt_tv, filt_tz): (Vec<Variant>, Vec<PhasedZygosity>) = truth_variants.iter().zip(truth_zygosity.iter())
                .filter_map(|(tv, &tz)| if tv.variant_type() == filter_type {
                    // put the clones here to reduce overhead
                    Some((tv.clone(), tz))
                } else {
                    None
                })
                .unzip();

            let (truth_tp, truth_fn) = if !filt_tv.is_empty() {
                // we have at least one variant left over in the query, so we need to add these stats
                let (filtered_truth, failed_ed) = generate_haplotype_sequence(
                    reference_genome, coordinates, &filt_tv, &filt_tz, haplotype
                ).with_context(|| format!("Error while generating sequence for filtered truth haplotype 1: {filter_type:?}"))?;
                let failed_vs = 2 * failed_ed as u64;

                let filt_align_metrics = perform_basepair_compare(
                    &reference[ref_start..ref_end],
                    &filtered_truth, &query_seq
                )?;

                (filt_align_metrics.truth_tp, filt_align_metrics.truth_fn + failed_vs)
            } else {
                (0, 0)
            };

            // we pull a piece from each comparison
            let vt_align_metrics = SummaryMetrics::new(
                truth_tp, truth_fn, query_tp, query_fp
            );
            debug!("\t{filter_type:?} => {vt_align_metrics:?}");
            bench_result.add_basepair_metrics(vt_align_metrics, Some(filter_type));
        }
    }

    Ok(())
}

/// Shortcut function for when we find a result that matches EXACTLY at the basepair level.
/// # Arguments
/// * `problem` - the original problem
/// * `truth_variants` - variants that are considered truth
/// * `truth_zygosity` - the zygosity of those variants
/// * `query_variants` - variants that are considered truth
/// * `query_zygosity` - the zygosity of those variants
/// * `reference_sequence` - the reference sequence for the region (not full chrom)
/// * `optimized_haps` - the haplotype sequences that were discovered
#[allow(clippy::too_many_arguments)]
fn generate_exact_match(
    problem: &CompareRegion,
    truth_variants: &[Variant], truth_zygosity: &[PhasedZygosity],
    query_variants: &[Variant], query_zygosity: &[PhasedZygosity],
    reference_sequence: &[u8], optimized_haps: &OptimizedHaplotypes
) -> anyhow::Result<CompareBenchmark> {
    // exact match, so ED = 0
    let mut bench_result = CompareBenchmark::new(problem.region_id(), 0, 0);

    // generate the expected zygosity counts and then add all of those exactly
    let expected_zyg_counts = generate_expected_zyg_counts(truth_zygosity);
    for (&ev, variant) in expected_zyg_counts.iter().zip(truth_variants.iter()) {
        debug!("\t{ev}, ={ev}, {variant:?}");
        bench_result.add_truth_zygosity(variant, ev, ev)
            .with_context(|| format!("Failed to add summary stats for variant {variant:?}"))?;
    }

    // do the same 
    let expected_query_counts = generate_expected_zyg_counts(query_zygosity);
    for (&ev, variant) in expected_query_counts.iter().zip(query_variants.iter()) {
        debug!("\t{ev}, ={ev}, {variant:?}");
        bench_result.add_query_zygosity(variant, ev, ev)
            .with_context(|| format!("Failed to add summary stats for variant {variant:?}"))?;
    }

    // REMEMBER: all basepair metrics are doubled to keep floating point away
    // first, the full length metrics should always match, and effectively mask redundant variants
    assert_eq!(optimized_haps.truth_seq1(), optimized_haps.query_seq1());
    assert_eq!(optimized_haps.truth_seq2(), optimized_haps.query_seq2());

    // figure out the base-pair level differences
    let ed1 = wfa_ed(reference_sequence, optimized_haps.truth_seq1())?;
    let ed2 = wfa_ed(reference_sequence, optimized_haps.truth_seq2())?;
    let joint_ed_score = 2 * (ed1 + ed2) as u64;

    // add them to both truth.TP and query.TP
    let shared_metrics = SummaryMetrics::new(joint_ed_score, 0, joint_ed_score, 0);
    bench_result.add_basepair_metrics(shared_metrics, None);

    // now add metrics for the truth variants by type
    for (&ev, variant) in expected_zyg_counts.iter().zip(truth_variants.iter())  {
        // figure out how many bases are correctly added, which is the different in REF and ALT sequences for the variant
        // we double everything here to avoid floating-point
        let ref_alt_dist = 2 * variant.alt_ed()? as u64;
        // create the metrics, scaling by expected zygosity
        let var_metrics = SummaryMetrics::new((ev as u64)*ref_alt_dist, 0, 0, 0);

        // add to the variant type
        let variant_type = variant.variant_type();
        bench_result.add_basepair_metrics(var_metrics, Some(variant_type));
    }

    // now add the query form metrics; in most cases these match, but different representations can create alternate variant-level metrics
    let expected_query_counts = generate_expected_zyg_counts(query_zygosity);
    for (&ev, variant) in expected_query_counts.iter().zip(query_variants.iter()) {
        // figure out how many bases are correctly added, which is the different in REF and ALT sequences for the variant
        // we double everything here to avoid floating-point
        let ref_alt_dist = 2 * variant.alt_ed()? as u64;
        // create the metrics, scaling by expected zygosity
        let var_metrics = SummaryMetrics::new(0, 0, (ev as u64)*ref_alt_dist, 0);

        // add to the variant type
        let variant_type = variant.variant_type();
        bench_result.add_basepair_metrics(var_metrics, Some(variant_type));
    }

    Ok(bench_result)
}

/// This will compare two sequences and flag any missing or different bases from the truth sequence.
/// Critically, there is no information indicating if a base is a "variant" or just buffer.
/// The values returned in the SummaryMetrics are doubled to account for half-correctness.
/// You get the same recall/prec/f1 since everything is doubled, but the numbers are inflated.
/// # Arguments
/// * `reference_sequence` - the reference sequence or baseline
/// * `truth_sequence` - the truth or baseline sequence; output TP + TN should equal the length of this sequence
/// * `query_sequence` - the query sequence, which is aligned to truth
fn perform_basepair_compare(
    reference_sequence: &[u8],
    truth_sequence: &[u8], query_sequence: &[u8]
) -> anyhow::Result<SummaryMetrics> {
    /*
    Core problem:
    - let truth_total = the number of bases in a region that are changed, inserted, or deleted in the truth set
        - one SNV = 1 truth total
        - pure insertion / deletion = L, where L is the number of bases added/deleted
        - indels are more complicated, but some combination of the above
        - more generally, edit_distance(REF allele v. ALT allele)
    - Find the number of bases that are correctly changed in the query also; this is truth.TP
    - Alternatively, find the number of bases that are NOT correctly changed in the query sequence; this is truth.FN
    - If we find one, we get the other for free
    - If we solve this operation, we can swap truth and query to get query.TP and query.FP

    Math form:
    - still need to figure out best way to do this per-variant; start with REF and progressively add variants to either truth or query
    - conceptually, make the reference-truth-query triangle and use a system of equations:
        - X = ED_ref_truth = TP+FN
        - Y = ED_ref_query = TP+FP
        - Z = ED_truth_query = FP+FN
        - Solve that to get:
        - TP = (X+Y-Z) / 2 - yes, you can get partial corrects; how to handle? floor it? ceil it? actual have partial correctness?
        - FN = X - TP
        - FP = Y - TP
    */
    // calculate the baseline differences; we double everything here to avoid floating-point; all results should *technically* be halved on the backend
    let ed_ref_truth = 2 * wfa_ed(reference_sequence, truth_sequence)? as u64;
    let ed_ref_query = 2 * wfa_ed(reference_sequence, query_sequence)? as u64;
    let ed_truth_query = 2 * wfa_ed(truth_sequence, query_sequence)? as u64;

    // at the basepair level, truth.TP and query.TP are strictly equal
    let tp_shared = (ed_ref_truth + ed_ref_query - ed_truth_query) / 2;

    // now derive these
    let truth_fn = ed_ref_truth - tp_shared;
    let query_fp = ed_ref_query - tp_shared;

    // set and return
    let summary_metrics = SummaryMetrics {
        truth_tp: tp_shared,
        truth_fn,
        query_tp: tp_shared,
        query_fp
    };
    Ok(summary_metrics)
}

/*
Another idea:
- One issue with some of the above approaches is that we are allowing FP/FN variants to "disappear" in the graph traversal
- We are flagging these correctly in the "best match" GT and Hap scoring; but when we compare at the basepair level, these are partially missing
- What if:
    - we assume that truth is phased, or force it to be phased (which we already do); these are two "baseline" haplotypes
    - we then do a dynamic WFA traversal of both, intelligently forcing the query GT to behave as called; i.e. 0|1 -> one get REF, one gets ALT
    - we put these into a lowest cost queue, where we are adding one query variant at a time, allowing for flips when possible
    - if we traversing both at once, we have a combined cost that we can traverse in lowest-ED order
    - the lowest ED full path is our "best" path where the two best formations of the query sequence are assigned to the two truth sequences
    - at that point, we can do direct sequence comparisons using `perform_basepair_compare`, and the output should be trustworthy
*/

/// Given a set of variants with phased zygosities, this will put together the selected haplotype sequence.
/// # Arguments
/// * `reference_genome` - the pre-loaded reference
/// * `region` - the region we are looking at in the reference
/// * `variants` - set of variants that might get added
/// * `zygosities` - the zygosity of the calls
/// * `haplotype` - tells us which alleles in the variant/zygosity combination to add to the sequence
/// # Errors
/// * if an Unknown zygosity is encountered
/// * if we encountered two non-reference variants on the same haplotype that conflict with each other
/// # Panics
/// * if the region specified is not in the reference genome
fn generate_haplotype_sequence(
    reference_genome: &ReferenceGenome, region: &Coordinates, variants: &[Variant], zygosities: &[PhasedZygosity], haplotype: Haplotype
) -> anyhow::Result<(Vec<u8>, usize)> {
    // convert the zygosity + haplotype phase into an allele series
    let mut alleles = Vec::<Allele>::with_capacity(zygosities.len());
    
    for zyg in zygosities.iter() {
        // turn zygosity into an allele pair
        let (a1, a2) = match *zyg {
            PhasedZygosity::Unknown => bail!("Unknown zygosity encountered in variant"),
            PhasedZygosity::HomozygousReference => (Allele::Reference, Allele::Reference),
            PhasedZygosity::UnphasedHeterozygous => (Allele::Reference, Allele::Alternate),
            PhasedZygosity::PhasedHet01 => (Allele::Reference, Allele::Alternate),
            PhasedZygosity::PhasedHet10 => (Allele::Alternate, Allele::Reference),
            PhasedZygosity::HomozygousAlternate => (Allele::Alternate, Allele::Alternate),
        };
    
        // select which allele based on the haplotype
        let allele = match haplotype {
            Haplotype::Hap1 => a1,
            Haplotype::Hap2 => a2,
        };
        alleles.push(allele);
    }

    // Now just call it for the alleles we found
    generate_allele_sequence(reference_genome, region, variants, &alleles)
}

/// Given a reference genome region and a set of variants with alleles, this will combine them into the full haplotype sequence.
/// It will also return the total edit distance of any variants that could not be incorporated due to variant incompatibility (e.g., overlaps).
/// # Arguments
/// * `reference_genome` - the pre-loaded reference
/// * `region` - the region we are looking at in the reference
/// * `variants` - set of variants that might get added
/// * `alleles` - the allele chain to put together
/// # Errors
/// * if an Unknown allele is encountered
/// * if we encountered two non-reference variants on the same haplotype that conflict with each other
/// # Panics
/// * if the region specified is not in the reference genome
fn generate_allele_sequence(
    reference_genome: &ReferenceGenome, region: &Coordinates, variants: &[Variant], alleles: &[Allele]
) -> anyhow::Result<(Vec<u8>, usize)> {
    if variants.len() != alleles.len() {
        bail!("variants and alleles must be the same length");
    }

    let chrom = region.chrom();
    let mut current_ref_position = region.start() as usize;
    let mut current_sequence = vec![];
    let mut failed_ed = 0;
    for (variant, &allele) in variants.iter().zip(alleles.iter()) {
        if allele == Allele::Reference && variant.convert_index(0) == 0 {
            // indicates that this allele matches the reference genome, lets skip it entirely
            continue;
        }
        
        // first, extend the current sequence up to the variant position
        let vpos = variant.position() as usize;
        if vpos < current_ref_position {
            // this indicates that the previous variant extended past the start of this variant; usually happens when you have an indel on one allele and something else on the other
            // we decided to handle conflicts by simply ignoring the second variant that triggers the conflict
            // that would be this variant, so just treat it as REF; will get converted into a FN downstream
            // we do want to track how many bases this impacts though
            assert_eq!(variant.convert_index(0), 0); // verify REF is in allele0
            failed_ed += variant.alt_ed()?;
            continue;
        }

        // fill in from the previous position to here
        current_sequence.extend_from_slice(
            reference_genome.get_slice(chrom, current_ref_position, vpos)
        );
        current_ref_position = vpos;

        // now add the appropriate sequence from the variant and increase our reference position to match the variant reference length
        let variant_sequence = match allele {
            Allele::Unknown => bail!("Allele::Unknown encountered"),
            Allele::Reference => variant.allele0(),
            Allele::Alternate => variant.allele1(),
        };
        current_sequence.extend_from_slice(variant_sequence);
        current_ref_position += variant.ref_len();
    }

    // fill in to the end now
    let end_pos = region.end() as usize;
    current_sequence.extend_from_slice(
        reference_genome.get_slice(chrom, current_ref_position, end_pos)
    );

    Ok((current_sequence, failed_ed))
}

/// For a set of zygosities, this will generate a flattened count of the alternate allele.
/// # Arguments
/// * `zygosities` - the zygosities
fn generate_expected_zyg_counts(zygosities: &[PhasedZygosity]) -> Vec<u8> {
    zygosities.iter()
        .map(|z| {
            match *z {
                PhasedZygosity::Unknown => todo!("handle Unknown zygosity"),
                PhasedZygosity::HomozygousReference => 0,
                PhasedZygosity::UnphasedHeterozygous |
                PhasedZygosity::PhasedHet01 |
                PhasedZygosity::PhasedHet10 => 1,
                PhasedZygosity::HomozygousAlternate => 2
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use crate::data_types::summary_metrics::{SummaryGtMetrics, SummaryMetrics};
    use crate::data_types::variant_metrics::{VariantMetrics, VariantSource};

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
    fn test_generate_haplotype_sequence_001() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();

        // simple SNP test
        let region = Coordinates::new("mock_chr1".to_string(), 5, 15); // TACCAGGACT
        let variants = [
            Variant::new_snv(0, 10, b"G".to_vec(), b"C".to_vec()).unwrap()
        ];
        let zygosities = [
            PhasedZygosity::PhasedHet10
        ];
        let (hap1_test, ed1) = generate_haplotype_sequence(
            &reference_genome, &region, &variants, &zygosities, Haplotype::Hap1
        ).unwrap();
        assert_eq!(&hap1_test, b"TACCACGACT"); // SNP at pos 10 in 5-15: TACCAGGACT -> TACCACGACT
        assert_eq!(ed1, 0);

        let (hap2_test, ed2) = generate_haplotype_sequence(
            &reference_genome, &region, &variants, &zygosities, Haplotype::Hap2
        ).unwrap();
        assert_eq!(&hap2_test, b"TACCAGGACT"); // no SNP, should match reference
        assert_eq!(ed2, 0);
    }

    #[test]
    fn test_generate_haplotype_sequence_002() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();

        // test with an indel and an overlapping SNP; these are on opposite haplotypes though so it should be okay
        let region = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let variants = [
            Variant::new_deletion(0, 3, b"GT".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 4, b"T".to_vec(), b"G".to_vec()).unwrap()
        ];
        let zygosities = [
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet01
        ];
        let (hap1_test, ed1) = generate_haplotype_sequence(
            &reference_genome, &region, &variants, &zygosities, Haplotype::Hap1
        ).unwrap();
        assert_eq!(&hap1_test, b"ACCGTACCA"); // ACCGTTACCA -> ACCGTACCA
        assert_eq!(ed1, 0);

        let (hap2_test, ed2) = generate_haplotype_sequence(
            &reference_genome, &region, &variants, &zygosities, Haplotype::Hap2
        ).unwrap();
        assert_eq!(&hap2_test, b"ACCGGTACCA"); // ACCGTTACCA -> ACCGGTACCA
        assert_eq!(ed2, 0);
    }

    #[test]
    fn test_generate_haplotype_sequence_conflict() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();

        // test with an indel and an overlapping SNP; if only one is present (hap1), it's fine; if both are present, we should get an error
        let region = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let variants = [
            Variant::new_deletion(0, 3, b"GT".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 4, b"T".to_vec(), b"G".to_vec()).unwrap()
        ];
        let zygosities = [
            PhasedZygosity::PhasedHet01,
            PhasedZygosity::HomozygousAlternate
        ];
        let (hap1_test, ed1) = generate_haplotype_sequence(
            &reference_genome, &region, &variants, &zygosities, Haplotype::Hap1
        ).unwrap();
        assert_eq!(&hap1_test, b"ACCGGTACCA"); // ACCGTTACCA -> ACCGGTACCA
        assert_eq!(ed1, 0);

        let (hap2_test, ed2) = generate_haplotype_sequence(
            &reference_genome, &region, &variants, &zygosities, Haplotype::Hap2
        ).unwrap();
        assert_eq!(&hap2_test, b"ACCGTACCA"); // ACCGTTACCA -> ACCGTACCA
        assert_eq!(ed2, 1);
    }

    /// Simple test with just a single SNV. The calls are an exact clone() of the truth inputs.
    #[test]
    fn test_solve_compare_region_simple_snv() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let truth_variants = vec![
            Variant::new_snv(0, 2, b"C".to_vec(), b"G".to_vec()).unwrap()
        ];
        let truth_zygosity = vec![
            PhasedZygosity::PhasedHet10
        ];
        let query_variants = truth_variants.clone();
        let query_zygosity = truth_zygosity.clone();
        let problem = CompareRegion::new(
            0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
        ).unwrap();

        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();
        assert_eq!(result.total_ed(), 0); // should be exact paths
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(1, 0, 1, 0, 0, 0));
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(1, 0, 1, 0));
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(2*1, 0, 2*1, 0));
        assert_eq!(result.truth_variant_data(), &[VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap()]);
        assert_eq!(result.query_variant_data(), &[VariantMetrics::new(VariantSource::Query, 1, 1).unwrap()]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACGGTTACCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACGGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq2, "ACCGTTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACCGTTACCA");
    }

    /// Test where the truth has two SNVs side-by-side, but the query set has one "indel" entry and the phase orientation is flipped.
    /// Should still be identical despite the representation differences.
    #[test]
    fn test_solve_compare_region_double_single() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let truth_variants = vec![
            Variant::new_snv(0, 2, b"C".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 3, b"G".to_vec(), b"T".to_vec()).unwrap()
        ];
        let truth_zygosity = vec![
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet10
        ];

        // has been collapsed into a single variant here, and also the phase orientation is swapped; should still be an exact match though
        let query_variants = vec![
            Variant::new_indel(0, 2, b"CG".to_vec(), b"GT".to_vec()).unwrap()
        ];
        let query_zygosity = vec![
            PhasedZygosity::PhasedHet01
        ];
        let problem = CompareRegion::new(
            0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
        ).unwrap();

        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();
        assert_eq!(result.total_ed(), 0); // should be exact paths
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(2, 0, 1, 0, 0, 0)); // two variants here
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(2, 0, 1, 0));
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(2*2, 0, 2*2, 0)); // 2bp difference is shared
        assert_eq!(result.truth_variant_data(), &[
            VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap(),
            VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap()
        ]);
        assert_eq!(result.query_variant_data(), &[VariantMetrics::new(VariantSource::Query, 1, 1).unwrap()]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACGTTTACCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACGTTTACCA");
        assert_eq!(&sequence_bundle.truth_seq2, "ACCGTTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACCGTTACCA");
    }

    /// Opposite of the previous test: truth has one variant entry, query has two.  Should still be identical.
    #[test]
    fn test_solve_compare_region_single_double() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let truth_variants = vec![
            Variant::new_indel(0, 2, b"CG".to_vec(), b"GT".to_vec()).unwrap()
        ];
        let truth_zygosity = vec![
            PhasedZygosity::PhasedHet01
        ];

        let query_variants = vec![
            Variant::new_snv(0, 2, b"C".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 3, b"G".to_vec(), b"T".to_vec()).unwrap()
        ];
        let query_zygosity = vec![
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet10
        ];

        let problem = CompareRegion::new(
            0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
        ).unwrap();

        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();
        assert_eq!(result.total_ed(), 0); // should be exact paths
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(1, 0, 2, 0, 0, 0)); // only one variant in truth set
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(1, 0, 2, 0));
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(2*2, 0, 2*2, 0)); // 2bp difference is shared
        assert_eq!(result.truth_variant_data(), &[VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap()]);
        assert_eq!(result.query_variant_data(), &[
            VariantMetrics::new(VariantSource::Query, 1, 1).unwrap(),
            VariantMetrics::new(VariantSource::Query, 1, 1).unwrap()
        ]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACCGTTACCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq2, "ACGTTTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACGTTTACCA");
    }

    /// Insertion of a C. In truth, it is left shifted, in query it is right shifted
    #[test]
    fn test_solve_compare_region_insertion_shift() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let truth_variants = vec![
            Variant::new_insertion(0, 0, b"A".to_vec(), b"AC".to_vec()).unwrap()
        ];
        let truth_zygosity = vec![
            PhasedZygosity::HomozygousAlternate
        ];

        let query_variants = vec![
            Variant::new_insertion(0, 2, b"C".to_vec(), b"CC".to_vec()).unwrap()
        ];
        let query_zygosity = vec![
            PhasedZygosity::HomozygousAlternate,
        ];

        let problem = CompareRegion::new(
            0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
        ).unwrap();

        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();
        assert_eq!(result.total_ed(), 0); // should be exact paths
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(1, 0, 1, 0, 0, 0)); // only one variant in truth set
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(2, 0, 2, 0)); // but two haplotype variants
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(2*2, 0, 2*2, 0)); // 1bp difference on both haps
        assert_eq!(result.truth_variant_data(), &[VariantMetrics::new(VariantSource::Truth, 2, 2).unwrap()]);
        assert_eq!(result.query_variant_data(), &[VariantMetrics::new(VariantSource::Query, 2, 2).unwrap()]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACCCGTTACCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACCCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq2, "ACCCGTTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACCCGTTACCA");
    }

    /// Simple test with just a single SNV. The calls are an exact clone() of the truth inputs.
    #[test]
    fn test_solve_compare_region_simple_fn() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let truth_variants = vec![
            Variant::new_snv(0, 2, b"C".to_vec(), b"G".to_vec()).unwrap()
        ];
        let truth_zygosity = vec![
            PhasedZygosity::PhasedHet10
        ];
        let query_variants = vec![];
        let query_zygosity = vec![];
        let problem = CompareRegion::new(
            0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
        ).unwrap();

        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();
        assert_eq!(result.total_ed(), 1); // 1 BP ed
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(0, 1, 0, 0, 0, 0));
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(0, 1, 0, 0));
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(0, 2*1, 0, 0)); // 1bp FN
        assert_eq!(result.truth_variant_data(), &[VariantMetrics::new(VariantSource::Truth, 1, 0).unwrap()]);
        assert_eq!(result.query_variant_data(), &[]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACGGTTACCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq2, "ACCGTTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACCGTTACCA");
    }

    #[test]
    fn test_solve_compare_region_complex_001() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let truth_variants = vec![
            Variant::new_snv(0, 2, b"C".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 4, b"T".to_vec(), b"C".to_vec()).unwrap(),
            Variant::new_snv(0, 6, b"A".to_vec(), b"C".to_vec()).unwrap()
        ];
        let truth_zygosity = vec![
            PhasedZygosity::PhasedHet10,
            PhasedZygosity::PhasedHet01,
            PhasedZygosity::PhasedHet10
        ];
        let query_variants = vec![
            Variant::new_snv(0, 2, b"C".to_vec(), b"G".to_vec()).unwrap(),
            // Variant::new_snv(0, 4, b"T".to_vec(), b"C".to_vec()).unwrap(),
            Variant::new_snv(0, 6, b"A".to_vec(), b"C".to_vec()).unwrap()
        ];
        let query_zygosity = vec![
            PhasedZygosity::HomozygousAlternate,
            PhasedZygosity::PhasedHet01
        ];
        let problem = CompareRegion::new(
            0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
        ).unwrap();

        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();
        assert_eq!(result.total_ed(), 2); // two BP delta
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(2, 1, 1, 1, 0, 1));
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(2, 1, 2, 1));
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(2*2, 2*1, 2*2, 2*1)); // 2 matching bp; 1 FN, 1 FP
        assert_eq!(result.truth_variant_data(), &[
            VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap(), // FP in truth space
            VariantMetrics::new(VariantSource::Truth, 1, 0).unwrap(), // FN in truth space
            VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap()
        ]);
        assert_eq!(result.query_variant_data(), &[
            VariantMetrics::new(VariantSource::Query, 1, 2).unwrap(), // FP in query space
            VariantMetrics::new(VariantSource::Query, 1, 1).unwrap()
        ]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACGGTTCCCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACGGTTCCCA");
        assert_eq!(&sequence_bundle.truth_seq2, "ACCGCTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACGGTTACCA");
    }

    #[test]
    fn test_solve_compare_region_ambiguous() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let coordinates = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let truth_variants = vec![
            Variant::new_snv(0, 4, b"T".to_vec(), b"C".to_vec()).unwrap()
        ];
        let truth_zygosity = vec![
            PhasedZygosity::HomozygousAlternate
        ];
        let query_variants = vec![
            Variant::new_snv(0, 4, b"T".to_vec(), b"C".to_vec()).unwrap(),
            Variant::new_snv(0, 4, b"T".to_vec(), b"A".to_vec()).unwrap()
        ];
        let query_zygosity = vec![
            PhasedZygosity::PhasedHet01,
            PhasedZygosity::PhasedHet10
        ];
        let problem = CompareRegion::new(
            0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity
        ).unwrap();

        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();
        assert_eq!(result.total_ed(), 1); // path through graph is an exact match still
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(0, 1, 1, 1, 1, 0));
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(1, 1, 1, 1));
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(3, 1, 3, 1)); // one allele is on track, the other is half-right, half-wrong
        assert_eq!(result.truth_variant_data(), &[
            VariantMetrics::new(VariantSource::Truth, 2, 1).unwrap() // expected homozygous, but found something else
        ]);
        assert_eq!(result.query_variant_data(), &[
            VariantMetrics::new(VariantSource::Query, 1, 1).unwrap(), // TP in query space since it started with 1
            VariantMetrics::new(VariantSource::Query, 0, 1).unwrap() // FP in query space
        ]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACCGCTACCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACCGATACCA"); // it puts A into q1
        assert_eq!(&sequence_bundle.truth_seq2, "ACCGCTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACCGCTACCA"); // it puts C into q2; these can technically swap
    }

    #[test]
    fn test_solve_compare_region_skip_variants() {
        // fixed simple reference we can use for these tests
        let reference_genome = generate_simple_reference();
        let region = Coordinates::new("mock_chr1".to_string(), 0, 10); // ACCGTTACCA
        let variants = vec![
            Variant::new_deletion(0, 3, b"GT".to_vec(), b"G".to_vec()).unwrap(),
            Variant::new_snv(0, 4, b"T".to_vec(), b"G".to_vec()).unwrap()
        ];
        let zygosities = vec![
            PhasedZygosity::PhasedHet01,
            PhasedZygosity::HomozygousAlternate
        ];

        let problem = CompareRegion::new(
            0, region, variants.clone(), zygosities.clone(), variants, zygosities
        ).unwrap();
        let result = solve_compare_region(&problem, &reference_genome, Default::default(), None).unwrap();

        // our results should have two of the SNVs skipped, one in truth and one in query; this is because they can't get incorporated
        assert_eq!(result.total_ed(), 0); // add ED comes from the variant skips
        let group_metrics = result.group_metrics();
        assert_eq!(*group_metrics.gt(), SummaryGtMetrics::new(1, 1, 1, 1, 1, 1)); // one right, one wrong in each
        assert_eq!(*group_metrics.hap(), SummaryMetrics::new(2, 1, 2, 1)); // two correct alleles, one skipped
        assert_eq!(*group_metrics.basepair(), SummaryMetrics::new(4, 2, 4, 2)); // same, but doubled for basepair
        assert_eq!(group_metrics.variant_basepair().get(&VariantType::Snv).unwrap(), &SummaryMetrics::new(3, 1, 3, 1)); // check that one SNV is missed, but half right in this situation
        assert_eq!(group_metrics.variant_basepair().get(&VariantType::Deletion).unwrap(), &SummaryMetrics::new(2, 0, 2, 0)); // check that one SNV is missed, but half right in this situation
        assert_eq!(result.truth_variant_data(), &[
            VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap(),
            VariantMetrics::new(VariantSource::Truth, 2, 1).unwrap()
        ]);
        assert_eq!(result.query_variant_data(), &[
            VariantMetrics::toggle_source(&VariantMetrics::new(VariantSource::Truth, 1, 1).unwrap()),
            VariantMetrics::toggle_source(&VariantMetrics::new(VariantSource::Truth, 2, 1).unwrap())
        ]);

        // check the sequences also
        let sequence_bundle = result.sequence_bundle().unwrap();
        assert_eq!(&sequence_bundle.ref_seq,    "ACCGTTACCA");
        assert_eq!(&sequence_bundle.truth_seq1, "ACCGGTACCA");
        assert_eq!(&sequence_bundle.query_seq1, "ACCGGTACCA");
        assert_eq!(&sequence_bundle.truth_seq2, "ACCGTACCA");
        assert_eq!(&sequence_bundle.query_seq2, "ACCGTACCA");
    }

    #[test]
    fn test_perform_basepair_compare_simple() {
        let result = perform_basepair_compare(
            b"ACGTACGT",
            b"ACGTACGT", // exact match to REF
            b"ACCTAGGT" // two false positives
        ).unwrap();
        assert_eq!(result, SummaryMetrics::new(0, 0, 0, 2 * 2));

        let result = perform_basepair_compare(
            b"ACGTACGT",
            b"ACCTAGGT", // two false negatives
            b"ACGTACGT" // exact match to REF
        ).unwrap();
        assert_eq!(result, SummaryMetrics::new(0, 2 * 2, 0, 0));

        let result = perform_basepair_compare(
            b"ACGTACGT",
            b"ACCTAGGT", // two TP
            b"ACCTAGGT" // exact match to truth, so two TP
        ).unwrap();
        assert_eq!(result, SummaryMetrics::new(2 * 2, 0, 2 * 2, 0));
    }

    #[test]
    fn test_perform_basepair_compare_mixed() {
        let result = perform_basepair_compare(
            b"TAT",
            b"TCT", // SNP
            b"TACT" // looks like neither, it's half-right
        ).unwrap();
        assert_eq!(result, SummaryMetrics::new(1, 1, 1, 1)); // the one situation where you don't get doubles, and thats because it's half-right, half-wrong
    }

    #[test]
    fn test_perform_basepair_compare_indels() {
        let result = perform_basepair_compare(
            b"TAT",
            b"TAAAAT", // big insertion
            b"TAAAT" // mostly correctly, but one FN base
        ).unwrap();
        assert_eq!(result, SummaryMetrics::new(2*2, 2*1, 2*2, 0));

        let result = perform_basepair_compare(
            b"TAT",
            b"TAAAAT", // big insertion
            b"TAAAAAT" // mostly correctly, but one FP base
        ).unwrap();
        assert_eq!(result, SummaryMetrics::new(2*3, 0, 2*3, 2*1));
    }
}