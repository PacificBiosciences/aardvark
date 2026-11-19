//! Times `solve_compare_region` on synthetic regions.
//! Each case is built before the timed loop.
//! Haplotype sequences are left off so storing them does not dominate the measurement.

use aardvark_bio::data_types::compare_region::CompareRegion;
use aardvark_bio::data_types::coordinates::Coordinates;
use aardvark_bio::data_types::phase_enums::PhasedZygosity;
use aardvark_bio::data_types::variants::Variant;
use aardvark_bio::waffle_solver::{solve_compare_region, CompareConfig, CompareConfigBuilder};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rust_lib_reference_genome::reference_genome::ReferenceGenome;

const CHROM: &str = "chr1";
const BASES: &[u8] = b"ACGT";
/// Share of variants dropped from each side in the random-drop case.
const DROP_PERCENT: u64 = 10;

/// One prepared comparison, ready to time.
struct CompareCase {
    label: String,
    genome: ReferenceGenome,
    problem: CompareRegion,
}

/// Deterministic xorshift64, so the benches need no extra crate.
struct Rng {
    state: u64,
}

impl Rng {
    /// Seeds the generator.
    /// Zero is a fixed point of xorshift, so the low bit is forced on.
    /// # Arguments
    /// * `seed` - any start state
    fn new(seed: u64) -> Self {
        Self { state: seed | 1 }
    }

    /// Next raw value.
    fn next_u64(&mut self) -> u64 {
        let mut state = self.state;
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        self.state = state;
        state
    }

    /// One random reference base.
    fn base(&mut self) -> u8 {
        BASES[(self.next_u64() % 4) as usize]
    }

    /// A base other than `current`.
    /// # Arguments
    /// * `current` - an ACGT base that must not be returned
    fn other_base(&mut self, current: u8) -> u8 {
        let choices: &[u8] = match current {
            b'A' => b"CGT",
            b'C' => b"AGT",
            b'G' => b"ACT",
            b'T' => b"ACG",
            _ => panic!("expected an ACGT base, got {current}"),
        };
        choices[(self.next_u64() % 3) as usize]
    }

    /// `len` random bases.
    /// # Arguments
    /// * `len` - number of bases to draw
    fn sequence(&mut self, len: usize) -> Vec<u8> {
        let mut sequence = Vec::with_capacity(len);
        for _ in 0..len {
            sequence.push(self.base());
        }
        sequence
    }

    /// True for about `percent` out of every 100 calls.
    /// # Arguments
    /// * `percent` - chance, from 0 to 100
    fn chance(&mut self, percent: u64) -> bool {
        self.next_u64() % 100 < percent
    }
}

/// Loads `sequence` as the only contig.
/// # Arguments
/// * `sequence` - uppercase ACGT reference
fn load_genome(sequence: &[u8]) -> ReferenceGenome {
    let mut genome = ReferenceGenome::empty_reference();
    let text = std::str::from_utf8(sequence).expect("reference sequence is ASCII");
    genome.add_contig(CHROM.to_string(), text).expect("contig");
    genome
}

/// Packages a region that covers the whole synthetic contig.
/// Every variant is a phased heterozygote on hap1.
/// # Arguments
/// * `label` - criterion function name
/// * `sequence` - reference the variants were cut from
/// * `truth_variants` - truth records, in position order
/// * `query_variants` - query records, in position order
fn compare_case(
    label: &str,
    sequence: &[u8],
    truth_variants: Vec<Variant>,
    query_variants: Vec<Variant>,
) -> CompareCase {
    let coordinates = Coordinates::new(CHROM.to_string(), 0, sequence.len() as u64);
    let truth_zygosity = vec![PhasedZygosity::PhasedHet10; truth_variants.len()];
    let query_zygosity = vec![PhasedZygosity::PhasedHet10; query_variants.len()];
    let problem = CompareRegion::new(
        0, coordinates, truth_variants, truth_zygosity, query_variants, query_zygosity,
    ).expect(label);
    CompareCase {
        label: label.to_string(),
        genome: load_genome(sequence),
        problem,
    }
}

/// SNV whose reference allele is the base at `position`.
/// # Arguments
/// * `sequence` - reference the allele is read from
/// * `position` - 0-based locus
/// * `rng` - source of the alternate base
fn snv_at(sequence: &[u8], position: u64, rng: &mut Rng) -> Variant {
    let reference = sequence[position as usize];
    let alternate = rng.other_base(reference);
    Variant::new_snv(0, position, vec![reference], vec![alternate]).expect("snv")
}

/// Builds a random reference and a list of small variants on it.
/// Mixes SNVs, 1-5 bp insertions, and 2-5 bp deletions, spaced so they do not overlap.
/// Returns the reference and the variants in position order.
/// # Arguments
/// * `count` - number of variants to place
/// * `seed` - makes the reference and variants repeatable
fn generate_small_variants(count: usize, seed: u64) -> (Vec<u8>, Vec<Variant>) {
    // one event every 16 bp, with room for a 5 bp deletion
    let stride = 16;
    let margin = 8;
    let mut rng = Rng::new(seed);
    let sequence = rng.sequence(margin + count * stride + margin);

    let mut variants = Vec::with_capacity(count);
    for index in 0..count {
        let position = (margin + index * stride) as u64;
        let pos = position as usize;
        let variant = match index % 3 {
            0 => snv_at(&sequence, position, &mut rng),
            1 => {
                // 1 to 5 inserted bases after the anchor
                let extra = 1 + (index % 5);
                let anchor = sequence[pos];
                let mut alternate = Vec::with_capacity(1 + extra);
                alternate.push(anchor);
                for _ in 0..extra {
                    alternate.push(rng.base());
                }
                Variant::new_insertion(0, position, vec![anchor], alternate).expect("insertion")
            },
            _ => {
                // 2 to 5 reference bases, anchored on the first
                let len = 2 + (index % 4);
                let reference = sequence[pos..pos + len].to_vec();
                let anchor = vec![reference[0]];
                Variant::new_deletion(0, position, reference, anchor).expect("deletion")
            },
        };
        variants.push(variant);
    }

    (sequence, variants)
}

/// Copy of `variants` that leaves each record out with probability `percent`.
/// # Arguments
/// * `variants` - the full list, in position order
/// * `percent` - chance each record is dropped, from 0 to 100
/// * `rng` - source of the drop decisions
fn drop_randomly(variants: &[Variant], percent: u64, rng: &mut Rng) -> Vec<Variant> {
    variants.iter()
        .filter(|_| !rng.chance(percent))
        .cloned()
        .collect()
}

/// Truth and query both carry every variant.
/// # Arguments
/// * `count` - number of small variants
fn exact(count: usize) -> CompareCase {
    let (sequence, variants) = generate_small_variants(count, 0xE5AC_0000 + count as u64);
    compare_case(&format!("exact_{count}"), &sequence, variants.clone(), variants)
}

/// Truth has every variant and query has none, so every record is a false negative.
/// # Arguments
/// * `count` - number of small variants
fn all_fn(count: usize) -> CompareCase {
    let (sequence, variants) = generate_small_variants(count, 0xFA11_0000 + count as u64);
    compare_case(&format!("all_fn_{count}"), &sequence, variants, vec![])
}

/// Query has every variant and truth has none, so every record is a false positive.
/// # Arguments
/// * `count` - number of small variants
fn all_fp(count: usize) -> CompareCase {
    let (sequence, variants) = generate_small_variants(count, 0xFB11_0000 + count as u64);
    compare_case(&format!("all_fp_{count}"), &sequence, vec![], variants)
}

/// Truth and query start from the same variants, then each independently loses `DROP_PERCENT` of them.
/// Records dropped from only one side become false negatives or false positives.
/// # Arguments
/// * `count` - number of small variants before dropping
fn random_drop(count: usize) -> CompareCase {
    let (sequence, variants) = generate_small_variants(count, 0xD50B_0000 + count as u64);
    let mut rng = Rng::new(0xD50C_0000 + count as u64);
    let truth_variants = drop_randomly(&variants, DROP_PERCENT, &mut rng);
    let query_variants = drop_randomly(&variants, DROP_PERCENT, &mut rng);
    compare_case(&format!("random_drop_{count}"), &sequence, truth_variants, query_variants)
}

/// Truth has one large deletion and query has nothing.
/// The two haplotypes differ by `len` bases, well past the edit distance where the wavefront gives way to the bit-parallel scan.
/// # Arguments
/// * `len` - reference length of the deleted stretch
fn missing_sv(len: usize) -> CompareCase {
    let position = 100usize;
    let mut rng = Rng::new(0xB16_0000 + len as u64);
    let sequence = rng.sequence(position + len + 100);

    let reference = sequence[position..position + len].to_vec();
    let anchor = vec![reference[0]];
    let deletion = Variant::new_sv_deletion(
        0, position as u64, reference, anchor,
    ).expect("sv deletion");

    compare_case(&format!("missing_sv_{len}"), &sequence, vec![deletion], vec![])
}

/// Every case the pull-request comparison reports.
fn benchmark_cases() -> Vec<CompareCase> {
    let mut cases = Vec::new();

    // identical small variants
    for count in [10, 100] {
        cases.push(exact(count));
    }

    // truth only
    for count in [10, 100] {
        cases.push(all_fn(count));
    }

    // query only; either phase of a lone query variant costs the same
    for count in [10, 100] {
        cases.push(all_fp(count));
    }

    // a tenth of each side dropped at random
    for count in [10, 100] {
        cases.push(random_drop(count));
    }

    // large truth deletion missing from the query
    for len in [10_000, 100_000] {
        cases.push(missing_sv(len));
    }

    cases
}

/// Compare configuration used for every case.
/// Sequences are disabled.
/// The branch cap is the same 50 that `CompareConfig::default` uses.
/// Unset builder fields would fall back to `Default` for the field type, which is 0 for the cap.
fn compare_config() -> CompareConfig {
    CompareConfigBuilder::default()
        .enable_sequences(false)
        .max_branch_factor(50)
        .enable_hpc_basepair_scoring(false)
        .build()
        .expect("compare config")
}

/// Times each prepared region.
fn bench_compare(c: &mut Criterion) {
    let config = compare_config();
    let cases = benchmark_cases();
    let mut group = c.benchmark_group("compare-group");

    for case in &cases {
        group.bench_function(&case.label, |b| {
            b.iter(|| {
                black_box(solve_compare_region(
                    &case.problem,
                    &case.genome,
                    config,
                    None,
                ).expect("solve_compare_region"));
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_compare);
criterion_main!(benches);
