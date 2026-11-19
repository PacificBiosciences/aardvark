
/*!
# Sequence alignment
The bit-parallel edit distance follows Myers (1999) in the form explained by Hyyrö (2001).
Adjacent scores differ by +1, 0, or −1, and one machine word advances 64 rows of the grid.
Hyyrö (2003) places those same deltas in diagonal blocks to test a distance threshold.
`advance_block` scans a full column instead.

Myers, G. 1999. A fast bit-vector algorithm for approximate string matching based on dynamic programming. Journal of the ACM 46(3): 395–415.
<https://doi.org/10.1145/316542.316550>
<https://www.gersteinlab.org/courses/452/09-spring/pdf/Myers.pdf>

Hyyrö, H. 2001. Explaining and extending the bit-parallel approximate string matching algorithm of Myers. University of Tampere, Department of Computer Science, report A-2001-10.
<http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.19.7158&rep=rep1&type=pdf>

Hyyrö, H. 2003. A bit-vector algorithm for computing Levenshtein and Damerau edit distances. Nordic Journal of Computing 10(1): 29–39.
<https://gwern.net/doc/cs/algorithm/2003-hyyro.pdf>
*/

/// Simply DWFA wrapper for edit distance on two known sequences.
/// This tends to be *much* faster than the `edit_distance(...)` approach.
/// # Arguments
/// * `v1` - the first sequence
/// * `v2` - the second sequence
/// # Errors
/// * if the underlying DWFALite has an issue
pub fn wfa_ed(v1: &[u8], v2: &[u8]) -> anyhow::Result<usize> {
    let mut dwfa = crate::dwfa::dynamic_wfa::DWFALite::default();
    dwfa.finalize(v1, v2)?;
    Ok(dwfa.edit_distance())
}

/// Same result as `wfa_ed`, given a lower bound on the edit distance.
/// The wavefront hands over to the bit-parallel scan once it would cost more, so a distance known to be past that point goes straight to the scan.
/// Otherwise this runs `wfa_ed` unchanged, so it is never slower.
/// # Arguments
/// * `v1` - the first sequence
/// * `v2` - the second sequence
/// * `lower_bound` - a value no larger than the edit distance between `v1` and `v2`
/// # Errors
/// * if the underlying DWFALite has an issue
pub fn wfa_ed_with_lower_bound(v1: &[u8], v2: &[u8], lower_bound: usize) -> anyhow::Result<usize> {
    // every unit of length difference costs one edit, so it is a lower bound too
    let lower_bound = lower_bound.max(v1.len().abs_diff(v2.len()));
    if crate::dwfa::dynamic_wfa::wavefront_exceeds_bitparallel(lower_bound, v1, v2) {
        // lower bound is too high, so go straight to the bit-parallel algorithm
        Ok(bitparallel_edit_distance(v1, v2))
    } else {
        // lower bound is low enough, so use the dynamic WFA
        wfa_ed(v1, v2)
    }
}

/// Global edit distance with unit costs, by the bit-parallel algorithm of Myers (1999) in the form explained by Hyyrö (2001).
/// The result matches `edit_distance`.
/// The shorter sequence is the pattern, so the cost is O(min(n, m) * max(n, m) / 64).
/// # Arguments
/// * `v1` - the first sequence
/// * `v2` - the second sequence
pub fn bitparallel_edit_distance(v1: &[u8], v2: &[u8]) -> usize {
    // the shorter sequence is the pattern (rows) and the longer one is the text (columns)
    let (pattern, text) = if v1.len() <= v2.len() { (v1, v2) } else { (v2, v1) };
    let pattern_len = pattern.len();
    if pattern_len == 0 {
        return text.len();
    }
    let blocks = pattern_len.div_ceil(64);

    // symbols that occur in the pattern get rows 1..; row 0 stays zero and stands for every other symbol
    let mut symbol = [0u16; 256];
    let mut symbols = 0u16;
    for &base in pattern.iter() {
        if symbol[base as usize] == 0 {
            symbols += 1;
            symbol[base as usize] = symbols;
        }
    }
    let mut pattern_eq = vec![0u64; (symbols as usize + 1) * blocks];
    for (index, &base) in pattern.iter().enumerate() {
        let row = symbol[base as usize] as usize;
        pattern_eq[row * blocks + index / 64] |= 1u64 << (index % 64);
    }

    // vertical deltas of column 0 are all +1, since D[i][0] = i
    let mut positive_vertical = vec![u64::MAX; blocks];
    let mut negative_vertical = vec![0u64; blocks];
    let last_block = blocks - 1;
    let last_bit_shift = (pattern_len - 1) % 64;
    let mut score = pattern_len as i64;
    for &base in text.iter() {
        let row = symbol[base as usize] as usize;
        let equality = &pattern_eq[row * blocks..][..blocks];

        // the top row grows by 1 per column (D[0][j] = j); the horizontal carry enters as (positive, negative)
        let (mut carry_positive, mut carry_negative) = (1u64, 0u64);
        for block in 0..last_block {
            let (carry_out_positive, carry_out_negative) = advance_block(
                equality[block],
                &mut positive_vertical[block],
                &mut negative_vertical[block],
                carry_positive,
                carry_negative,
            );
            carry_positive = carry_out_positive;
            carry_negative = carry_out_negative;
        }

        // the last block holds the score row at bit (pattern_len - 1) % 64
        let equality_last = equality[last_block];
        let positive = positive_vertical[last_block];
        let negative = negative_vertical[last_block];
        let horizontal = horizontal_delta(equality_last | carry_negative, positive, negative);
        score += ((horizontal.0 >> last_bit_shift) & 1) as i64 - ((horizontal.1 >> last_bit_shift) & 1) as i64;
        let _ = advance_block(
            equality_last,
            &mut positive_vertical[last_block],
            &mut negative_vertical[last_block],
            carry_positive,
            carry_negative,
        );
    }
    debug_assert!(score >= 0);
    score as usize
}

/// One 64-row step of the bit-parallel recurrence.
/// Returns the horizontal delta's outgoing carry as `(positive, negative)`.
/// # Arguments
/// * `pattern_match` - equality bits of this block for the current text symbol
/// * `positive_vertical` - positive vertical delta bits, updated in place
/// * `negative_vertical` - negative vertical delta bits, updated in place
/// * `carry_positive` - incoming positive horizontal carry (0 or 1)
/// * `carry_negative` - incoming negative horizontal carry (0 or 1)
pub fn advance_block(
    pattern_match: u64,
    positive_vertical: &mut u64,
    negative_vertical: &mut u64,
    carry_positive: u64,
    carry_negative: u64,
) -> (u64, u64) {
    let (positive_horizontal, negative_horizontal) = horizontal_delta(
        pattern_match | carry_negative,
        *positive_vertical,
        *negative_vertical,
    );
    let carry_out = (positive_horizontal >> 63, negative_horizontal >> 63);
    let positive_horizontal = (positive_horizontal << 1) | carry_positive;
    let negative_horizontal = (negative_horizontal << 1) | carry_negative;
    let match_or_negative = pattern_match | *negative_vertical;
    *positive_vertical = negative_horizontal | !(match_or_negative | positive_horizontal);
    *negative_vertical = positive_horizontal & match_or_negative;
    carry_out
}

/// Horizontal delta bits `(positive, negative)` before the incoming carry is shifted in.
/// # Arguments
/// * `match_or_carry` - equality bits with the incoming negative carry already folded in
/// * `positive_vertical` - positive vertical delta bits
/// * `negative_vertical` - negative vertical delta bits
pub fn horizontal_delta(match_or_carry: u64, positive_vertical: u64, negative_vertical: u64) -> (u64, u64) {
    let horizontal_xor = (((match_or_carry & positive_vertical).wrapping_add(positive_vertical)) ^ positive_vertical) | match_or_carry;
    let positive = negative_vertical | !(horizontal_xor | positive_vertical);
    let negative = positive_vertical & horizontal_xor;
    (positive, negative)
}

/// Returns the edit distance between two u8 Vecs by doing the full grid calculation.
/// This version is row-based (rows are length of v1) for the main loop.
/// # Arguments
/// * `v1` - the first sequence
/// * `v2` - the second sequence
pub fn edit_distance(v1: &[u8], v2: &[u8]) -> usize {
    // structured such that each "row" is the length of v1 (i.e. v1 is conceptually on the x-axis)
    let l1: usize = v1.len();
    let mut row: Vec<usize> = vec![0; l1+1];
    let mut prev_row: Vec<usize> = (0..l1+1).collect();
    
    // go through each row
    for (i, &c2) in v2.iter().enumerate() {
        row[0] = i+1;
        for (j, &c1) in v1.iter().enumerate() {
            row[j+1] = [
                // skip a character in v2
                prev_row[j+1]+1,
                // skip a character in v1 
                row[j]+1,
                // diagonal match/mismatch
                prev_row[j]+({
                    if c1 == c2 {
                        0
                    } else {
                        1
                    }
                })
            ].into_iter().min().unwrap();
        }

        // swap the rows at the end of each iteration
        std::mem::swap(&mut row, &mut prev_row);
    }

    prev_row[l1]
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_edit_distance() {
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
        let v3: Vec<u8> = vec![1, 2, 3, 5];
        let v4: Vec<u8> = vec![];

        assert_eq!(edit_distance(&v1, &v1), 0);
        assert_eq!(edit_distance(&v1, &v2), 1);
        assert_eq!(edit_distance(&v1, &v3), 2);
        assert_eq!(edit_distance(&v1, &v4), 5);

        assert_eq!(edit_distance(&v2, &v2), 0);
        assert_eq!(edit_distance(&v2, &v3), 3);
        assert_eq!(edit_distance(&v2, &v4), 5);
        
        assert_eq!(edit_distance(&v3, &v3), 0);
        assert_eq!(edit_distance(&v3, &v4), 4);

        assert_eq!(edit_distance(&v4, &v4), 0);
    }

    /// Any valid lower bound gives the same distance as the full grid, including bounds past the bit-parallel hand-over.
    #[test]
    fn test_wfa_ed_with_lower_bound_matches_full_grid() {
        let mut state: u64 = 0x2545F4914F6CDD1D;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        for case in 0..300 {
            let length_1 = (next() % 400) as usize;
            let sequence_1: Vec<u8> = (0..length_1).map(|_| b"ACGT"[(next() % 4) as usize]).collect();

            // a related sequence in most cases, and an unrelated one in the rest
            let sequence_2: Vec<u8> = if case % 3 == 0 {
                (0..(next() % 400) as usize).map(|_| b"ACGT"[(next() % 4) as usize]).collect()
            } else {
                sequence_1.iter().filter(|_| next() % 10 != 0).copied().collect()
            };

            let expected = edit_distance(&sequence_1, &sequence_2);
            for lower_bound in [0, expected / 2, expected] {
                assert_eq!(wfa_ed_with_lower_bound(&sequence_1, &sequence_2, lower_bound).unwrap(), expected, "case {case} bound {lower_bound}");
            }
        }
    }

    #[test]
    fn test_wfa_ed() {
        let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
        let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
        let v3: Vec<u8> = vec![1, 2, 3, 5];
        let v4: Vec<u8> = vec![];

        assert_eq!(wfa_ed(&v1, &v1).unwrap(), 0);
        assert_eq!(wfa_ed(&v1, &v2).unwrap(), 1);
        assert_eq!(wfa_ed(&v1, &v3).unwrap(), 2);
        assert_eq!(wfa_ed(&v1, &v4).unwrap(), 5);

        assert_eq!(wfa_ed(&v2, &v2).unwrap(), 0);
        assert_eq!(wfa_ed(&v2, &v3).unwrap(), 3);
        assert_eq!(wfa_ed(&v2, &v4).unwrap(), 5);

        assert_eq!(wfa_ed(&v3, &v3).unwrap(), 0);
        assert_eq!(wfa_ed(&v3, &v4).unwrap(), 4);

        assert_eq!(wfa_ed(&v4, &v4).unwrap(), 0);
    }

    #[test]
    fn test_edit_error_001() {
        let v1 = [65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 67, 65, 65, 65];
        let v2 = [65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 67, 65, 65, 65, 65, 65, 65, 67, 65, 65, 65];
        let v3 = [65, 65, 65, 65, 65, 65, 65, 65, 65, 65,     65, 65, 65, 65, 65, 65, 67, 65, 65, 65];

        assert_eq!(edit_distance(&v1, &v3), 1);
        assert_eq!(edit_distance(&v2, &v3), 1);
        assert_eq!(edit_distance(&v3, &v1), 1);
        assert_eq!(edit_distance(&v3, &v2), 1);

        assert_eq!(wfa_ed(&v1, &v3).unwrap(), 1);
        assert_eq!(wfa_ed(&v2, &v3).unwrap(), 1);
        assert_eq!(wfa_ed(&v3, &v1).unwrap(), 1);
        assert_eq!(wfa_ed(&v3, &v2).unwrap(), 1);
    }

    /// The bit-parallel distance matches the full grid, including empty sequences, block borders, and bytes outside ACGT.
    #[test]
    fn test_bitparallel_edit_distance_matches_full_grid() {
        assert_eq!(bitparallel_edit_distance(b"", b""), 0);
        assert_eq!(bitparallel_edit_distance(b"", b"ACGT"), 4);
        assert_eq!(bitparallel_edit_distance(b"ACGT", b""), 4);

        // lengths on either side of the 64-bit block size, both orders
        for length in [1usize, 63, 64, 65, 127, 128, 129] {
            let pattern: Vec<u8> = (0..length).map(|i| b"ACGTN"[i % 5]).collect();
            let mut text = pattern.clone();
            text.push(b'A');
            assert_eq!(bitparallel_edit_distance(&pattern, &text), edit_distance(&pattern, &text), "length {length}");
            assert_eq!(bitparallel_edit_distance(&text, &pattern), edit_distance(&pattern, &text), "length {length} swapped");
        }

        // every byte value, so the symbol table has to hold all 256
        let all_bytes: Vec<u8> = (0..=255).collect();
        let mut reversed = all_bytes.clone();
        reversed.reverse();
        assert_eq!(bitparallel_edit_distance(&all_bytes, &reversed), edit_distance(&all_bytes, &reversed));

        // deterministic pseudo-random sequences over ACGT and N
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        let alphabet = b"ACGTN";
        for case in 0..400 {
            let length_1 = (next() % 200) as usize;
            let length_2 = if case % 4 == 0 { (next() % 3) as usize } else { (next() % 200) as usize };
            let alphabet_size = if case % 2 == 0 { 2 } else { 5 };
            let sequence_1: Vec<u8> = (0..length_1).map(|_| alphabet[(next() % alphabet_size) as usize]).collect();
            let sequence_2: Vec<u8> = (0..length_2).map(|_| alphabet[(next() % alphabet_size) as usize]).collect();
            assert_eq!(bitparallel_edit_distance(&sequence_1, &sequence_2), edit_distance(&sequence_1, &sequence_2), "case {case}");
        }
    }
}