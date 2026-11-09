
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

/// Global (end-to-end) edit distance with unit costs, by the bit-parallel algorithm of Myers (1999) in the
/// block form of Hyyro (2003), as edlib computes its NW mode: O(|v1| * |v2| / 64). Same value as `edit_distance(...)`.
/// # Arguments
/// * `v1` - the first sequence
/// * `v2` - the second sequence
pub fn bitparallel_edit_distance(v1: &[u8], v2: &[u8]) -> usize {
    // the shorter sequence is the pattern (rows), the longer one the text (columns)
    let (pattern, text) = if v1.len() <= v2.len() { (v1, v2) } else { (v2, v1) };
    let m = pattern.len();
    if m == 0 {
        return text.len();
    }
    let blocks = m.div_ceil(64);
    // a compact alphabet: the symbols of the pattern get rows 1..; row 0 (all zero) is for every other symbol
    let mut symbol = [0u8; 256];
    let mut symbols = 0usize;
    for &c in pattern.iter() {
        if symbol[c as usize] == 0 {
            symbols += 1;
            symbol[c as usize] = symbols as u8;
        }
    }
    let mut peq: Vec<u64> = vec![0; (symbols + 1) * blocks];
    for (i, &c) in pattern.iter().enumerate() {
        peq[symbol[c as usize] as usize * blocks + i / 64] |= 1u64 << (i % 64);
    }
    // vertical deltas of column 0 are all +1 (D[i][0] = i)
    let mut pv: Vec<u64> = vec![u64::MAX; blocks];
    let mut mv: Vec<u64> = vec![0; blocks];
    let last = blocks - 1;
    let last_bit_shift = (m - 1) % 64;
    // the score of the last pattern row, D[m][j]
    let mut score = m as i64;
    for &c in text.iter() {
        let eq_row = &peq[symbol[c as usize] as usize * blocks..][..blocks];
        // global mode: the top row grows by +1 per column (D[0][j] = j); hin as (positive, negative) bits
        let (mut hin_p, mut hin_m) = (1u64, 0u64);
        for b in 0..last {
            let (eq0, p, mm) = (eq_row[b], pv[b], mv[b]);
            let xv = eq0 | mm;
            let eq = eq0 | hin_m;
            let xh = (((eq & p).wrapping_add(p)) ^ p) | eq;
            let ph = mm | !(xh | p);
            let mh = p & xh;
            let (ph_out, mh_out) = (ph >> 63, mh >> 63);
            let ph = (ph << 1) | hin_p;
            let mh = (mh << 1) | hin_m;
            pv[b] = mh | !(xv | ph);
            mv[b] = ph & xv;
            hin_p = ph_out;
            hin_m = mh_out;
        }
        // the last block: its score row is bit (m - 1) % 64
        let (eq0, p, mm) = (eq_row[last], pv[last], mv[last]);
        let xv = eq0 | mm;
        let eq = eq0 | hin_m;
        let xh = (((eq & p).wrapping_add(p)) ^ p) | eq;
        let ph = mm | !(xh | p);
        let mh = p & xh;
        score += ((ph >> last_bit_shift) & 1) as i64 - ((mh >> last_bit_shift) & 1) as i64;
        let ph = (ph << 1) | hin_p;
        let mh = (mh << 1) | hin_m;
        pv[last] = mh | !(xv | ph);
        mv[last] = ph & xv;
    }
    score as usize
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

    #[test]
    fn test_bitparallel_edit_distance_equals_full_grid() {
        // deterministic pseudo-random sequences over ACGT (and N), lengths across the 64-bit block borders
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        let alphabet = b"ACGTN";
        for case in 0..400 {
            let l1 = (next() % 200) as usize;
            let l2 = if case % 4 == 0 { (next() % 3) as usize } else { (next() % 200) as usize };
            let k = if case % 2 == 0 { 2 } else { 5 };
            let v1: Vec<u8> = (0..l1).map(|_| alphabet[(next() % k) as usize]).collect();
            let mut v2: Vec<u8> = v1.iter().filter(|_| next() % 10 != 0).copied().collect();
            v2.truncate(l2.max(v2.len().min(l2 + 50)));
            if case % 3 == 0 {
                v2 = (0..l2).map(|_| alphabet[(next() % k) as usize]).collect();
            }
            assert_eq!(bitparallel_edit_distance(&v1, &v2), edit_distance(&v1, &v2), "case {case}");
            assert_eq!(bitparallel_edit_distance(&v2, &v1), edit_distance(&v1, &v2), "case {case} swapped");
        }
    }
}