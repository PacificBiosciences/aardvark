/*!
# DWFA module
Contains Dynamic WaveFront Alignment implementations that enable fast quantification of edit distance between two dynamic sequences.
The core implementation is the `DWFALite` struct, which is a wrapper around the original DWFA from [waffle_con](https://github.com/PacificBiosciences/waffle_con/blob/main/src/dynamic_wfa.rs).
Additionally, the `HaplotypeDWFA` struct is a wrapper that facilitates the construction of two comparator haplotypes at once, which is core to the Aardvark use case of optimized variant phasing for comparison.

## Example DWFALite usage
```rust
use aardvark_bio::dwfa::dynamic_wfa::DWFALite;

// initialize with the same short sequence
let mut dwfa = DWFALite::default();
dwfa.update(b"ACGT", b"ACGT").unwrap();
assert_eq!(dwfa.edit_distance(), 0);

// simulate a single mismatch (SNV) difference between the two sequences
dwfa.update(b"ACGTA", b"ACGTC").unwrap();
assert_eq!(dwfa.edit_distance(), 1);

// finalize with the remaining sequence
dwfa.finalize(b"ACGTACGT", b"ACGTCCGT").unwrap();
assert_eq!(dwfa.edit_distance(), 1);
```
*/
/// Tweak of the original DWFA from waffle_con, we need some minor, but non-standard, adjustments for our use case
pub mod dynamic_wfa;
/// Wrapper that facilitates the construction of two comparator haplotypes at once. Also wraps the DWFALite for them.
pub mod haplotype_dwfa;
