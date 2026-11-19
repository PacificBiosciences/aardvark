
use rustc_hash::FxHashMap as HashMap;

use crate::util::boundary_alignment::BoundaryTracker;

#[derive(thiserror::Error, Debug)]
pub enum DWFAError {
    #[error("maximum edit distance exceeded")]
    MaxEditDistance,
    #[error("DWFA is already finalized, it cannot be changed further")]
    AlreadyFinalized,
    // this is a pass-through error from the boundary alignment code
    #[error(transparent)]
    Boundary(anyhow::Error),
}

/// The core dynamic WFA structure for a lite implementation.
/// It is structured such that all the sequences being built are maintained **outside** of this struct (hence the "lite").
/// Essentially, it is keep the wavefront information, but all sequences are updated and tracked elsewhere.
/// As a result, all updates must pass the sequences so the information can get tracked.
/// If the outside sequences are changed in any way other than appending, then this may become desynchronized.
/// Conceptually, if this is a 2D grid, the baseline sequence will go from top to bottom (y-axis) and the other sequence will go from left to right (x-axis).
/// This means each character we add to `other_seq` will add a new _column_ to the grid, and each character added to `baseline_seq` will add a new _row_.
/// 
/// This specialized version will run until the end of both baseline and other, knowing they're both dynamic.
/// We removed some of the specialized functionality from the main library to keep this simple.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct DWFALite {
    /// The minimum edit distance from `other_seq` to `baseline_seq` allowing for a large port of the tail of the `baseline_seq` to get ignored.
    edit_distance: usize,
    /// This is the core wavefront for the current `edit_distance`.
    /// It should always be length 2*`edit_distance`-1.
    /// For our mental model, the median diagonal means we have used the same number of bases in our baseline and other.
    /// Anything below the median diagonal means we have used more baseline than other.
    /// Anything above the median diagonal means we have used more other than baseline.
    /// The value stored will always be the number of bases consumed in `other_seq`.
    /// After any character addition, the wavefront should be touching one (or both) of the final column/rows in the grid.
    /// Finalizing, will force it into the final corner, touching both.
    wavefront: Vec<usize>,
    /// Farthest offset reached into the baseline sequence.
    max_baseline_distance: usize,
    /// Farthest offset reached into the other sequence.
    max_other_distance: usize,
    /// If true, this DWFA has been finalized, meaning it cannot be extended anymore.
    is_finalized: bool,
    /// Once the wavefront would be too expensive to grow, the exact grid boundaries take over from it.
    /// The wavefront fields above are no longer updated after that.
    boundary: Option<BoundaryTracker>,
    /// The maximum allowed edit distance before throwing an error
    max_edit_distance: usize
}

/// This is the *minimum* edit distance at which we *may* switch to the grid boundaries.
/// This is not a guaranteed switch at this edit distance.
const MIN_BOUNDARY_EDIT_DISTANCE: usize = 200;

impl Default for DWFALite {
    fn default() -> Self {
        DWFALite {
            edit_distance: 0,
            wavefront: vec![0],
            max_baseline_distance: 0,
            max_other_distance: 0,
            is_finalized: false,
            boundary: None,
            max_edit_distance: usize::MAX
        }
    }
}

/// Returns true if growing the wavefront to `next` costs more than one bit-parallel scan of both sequences.
/// The wavefront costs O(ED^2) and the bit-parallel distance costs O(n * m / 64).
/// One bit-parallel word is worth about four diagonal steps, so the switch is at ED^2 > n * m / 256.
/// # Arguments
/// * `next` - the edit distance the wavefront would grow to
/// * `baseline_seq` - the baseline sequence
/// * `other_seq` - the other sequence
pub fn wavefront_exceeds_bitparallel(next: usize, baseline_seq: &[u8], other_seq: &[u8]) -> bool {
    // WFA cost is ED^2
    let wfa_cost = next.saturating_mul(next);

    // bit-parallel cost is O(n * m / 256)
    let l1 = baseline_seq.len().max(1);
    let l2 = other_seq.len().max(1);
    let bitparallel_cost = l1.saturating_mul(l2) / 256;

    // switch when the WFA cost is more than one bit-parallel word
    wfa_cost > bitparallel_cost
}
    

impl DWFALite {
    pub fn with_max_edit_distance(max_edit_distance: usize) -> Self {
        Self {
            max_edit_distance,
            ..Default::default()
        }
    }

    /// This is the main function to extend the WFA with a new symbol.
    /// This function will handle all of the extending and potential edit distance increases as well.
    /// Once the wavefront would cost more than a bit-parallel scan, the exact grid boundaries take over and the result is unchanged.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If this function is called after `finalize()` has been called.
    /// * If we exceed the maximum edit distance
    /// * If either sequence is shorter than what the grid boundaries already consumed
    pub fn update(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, DWFAError> {
        if self.is_finalized {
            return Err(DWFAError::AlreadyFinalized);
        }

        // if we are already using the grid boundaries, direct straight to it
        if self.boundary.is_some() {
            return self.update_boundary(baseline_seq, other_seq);
        }

        // maximally extend everything along the current diagonals
        self.extend(baseline_seq, other_seq)?;
        
        // check how it looks
        while !self.reached_baseline_end(baseline_seq) && !self.reached_other_end(other_seq) {
            // check if we should switch to the grid boundaries
            let next = self.edit_distance.saturating_add(1);
            if next > MIN_BOUNDARY_EDIT_DISTANCE && wavefront_exceeds_bitparallel(next, baseline_seq, other_seq) {
                // the wavefront is too expensive to grow, so switch to the grid boundaries
                return self.enter_boundary(baseline_seq, other_seq);
            }

            // increase the edit distance, re-extension happens automatically
            self.increase_edit_distance(baseline_seq, other_seq)?;
        }


        Ok(self.edit_distance)
    }

    /// Builds the grid boundaries for the current sequences and uses them from now on.
    /// The distance returned is the same one the wavefront would have reached.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    /// # Errors
    /// * If we exceed the maximum edit distance
    /// * If either sequence is shorter than what the grid boundaries already consumed
    fn enter_boundary(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, DWFAError> {
        self.boundary = Some(BoundaryTracker::default());
        self.update_boundary(baseline_seq, other_seq)
    }

    /// Extends the grid boundaries to the current sequences.
    /// The distance is the smallest value on the last row or last column, which is when the wavefront reaches an end of either sequence.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    /// # Errors
    /// * If we exceed the maximum edit distance
    /// * If either sequence is shorter than what was already consumed
    fn update_boundary(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, DWFAError> {
        let tracker = self.boundary.as_mut()
            .expect("boundary tracker is active");
        tracker.extend(baseline_seq, other_seq)
            .map_err(DWFAError::Boundary)?;
        let distance = tracker.boundary_minimum();
        debug_assert!(distance >= self.edit_distance);

        // the wavefront reports a cap error with the distance left at cap + 1
        if distance > self.max_edit_distance {
            self.edit_distance = self.max_edit_distance + 1;
            return Err(DWFAError::MaxEditDistance);
        }
        self.edit_distance = distance;
        Ok(distance)
    }

    /// Extends one diagonal until it mismatches or reaches either sequence end.
    /// The farthest offsets reached on that diagonal are folded into the cached maxima.
    /// # Arguments
    /// * `index` - the wavefront diagonal to extend
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    fn extend_diagonal(&mut self, index: usize, baseline_seq: &[u8], other_seq: &[u8]) {
        // example at edit distance 1:
        // baseline offset would equal: [d+1, d, d-1] which corresponds to
        // a skipped base in primary, equal usage, and a skipped base in other

        // if we truncate a wavefront (e.g., a bounded size); then `i` below will be too small for every wavefront from the start
        //     that gets truncated; this means we want the below formula to be `*d + self.edit_distance - (i + truncate_prefix)`
        //     where `truncate_prefix` is the number we have cut off the front (in total).
        //     IMPORTANT: this formula will need to get updated everywhere that `baseline_offset` is calculated

        let offset = self.wavefront[index];
        let mut baseline_offset = offset + self.edit_distance - index;
        let mut other_offset = offset;

        // an offset already past either end stops here, including one that wrapped below zero
        if baseline_offset < baseline_seq.len() && other_offset < other_seq.len() {
            // compare 8 bytes at a time; the first differing byte ends the extension, the same as comparing byte by byte
            while baseline_offset + 8 <= baseline_seq.len() && other_offset + 8 <= other_seq.len() {
                // convert them into u64s for 8 symbols at a time
                let baseline_word = u64::from_le_bytes(baseline_seq[baseline_offset..baseline_offset + 8].try_into().unwrap());
                let other_word = u64::from_le_bytes(other_seq[other_offset..other_offset + 8].try_into().unwrap());

                // compare the two words
                let diff = baseline_word ^ other_word;
                if diff != 0 {
                    // the lowest set bit is the first differing byte on a little-endian read
                    let matching = (diff.trailing_zeros() / 8) as usize;
                    baseline_offset += matching;
                    other_offset += matching;
                    break;
                }

                // both matched, so move on to the next 8 bytes
                baseline_offset += 8;
                other_offset += 8;
            }

            // finish the tail, and a word that ended at a mismatch, byte by byte
            while baseline_offset < baseline_seq.len() &&
                other_offset < other_seq.len() &&
                baseline_seq[baseline_offset] == other_seq[other_offset] {
                baseline_offset += 1;
                other_offset += 1;
            }
        }

        // we are past the end of either sequence OR the sequences are not equal at this position, so we are done extending
        self.wavefront[index] = other_offset;
        self.max_baseline_distance = self.max_baseline_distance.max(baseline_offset);
        self.max_other_distance = self.max_other_distance.max(other_offset);
    }

    /// This will take the current wavefront and try to extend each one along it's current diagonal.
    /// In most cases, this will not do much work.
    /// This should also be stable such that calling it after finalizing has no impact.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * None so far
    fn extend(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<(), DWFAError> {
        // this is easier logic than trying to handle the option syntax below
        let len = self.wavefront.len();
        for i in 0..len {
            // `i` is the index in the wavefront
            // `i // 2` is always the middle diagonal
            // anything less than that has used more baseline symbols
            // anything more than that has used more other symbols
            // anything above or below needs to shift the comparison
            self.extend_diagonal(i, baseline_seq, other_seq);
        }
        Ok(())
    }

    /// This will increase the edit distance for this DWFA and create a new larger wavefront.
    /// The existing allocation grows in place, and each diagonal is extended once its new offset is final.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    /// # Errors
    /// * If the DWFA is already finalized
    /// * If we exceed the maximum edit distance
    fn increase_edit_distance(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<(), DWFAError> {
        if self.is_finalized {
            return Err(DWFAError::AlreadyFinalized);
        }

        // first, increase the distance we're at
        self.edit_distance += 1;
        if self.edit_distance > self.max_edit_distance {
            return Err(DWFAError::MaxEditDistance);
        }

        // grow the existing allocation; the two new diagonals start at zero
        let old_len = self.wavefront.len();
        self.wavefront.resize(old_len + 2, 0);

        // scatter from the end so each old offset is still visible
        for i in (0..old_len).rev() {
            let d = self.wavefront[i];

            // deletion (skipping) of a base in `baseline_seq`; this does not change the distance into `other_seq`
            self.wavefront[i] = self.wavefront[i].max(d);

            // mismatch progresses both sequences
            self.wavefront[i + 1] = self.wavefront[i + 1].max(d + 1);

            // insertion of a base into `baseline_seq`; this progresses `other_seq`
            self.wavefront[i + 2] = self.wavefront[i + 2].max(d + 1);

            // the diagonal two steps ahead has now received every contribution
            self.extend_diagonal(i + 2, baseline_seq, other_seq);
        }

        // the first two diagonals are final only after the lowest old offset is scattered
        self.extend_diagonal(0, baseline_seq, other_seq);
        self.extend_diagonal(1, baseline_seq, other_seq);
        Ok(())
    }

    /// This function is similar to `update(...)`, but it will perform all the final steps for an end-to-end calculation.
    /// This will trigger the algorithm to make sure we have reached the end of both the `baseline_seq` and `other_seq`, potentially increasing edit distance further.
    /// After calling this, no further extensions can be made, it is a "finalized" calculation.
    /// Only `edit_distance()` is meaningful afterwards.
    /// When the wavefront would be expensive the distance is computed directly, and the wavefront is left where it stopped.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    /// # Errors
    /// * If the edit distance cannot be increased further and it needs to be.
    /// * If either sequence is shorter than what the grid boundaries already consumed
    pub fn finalize(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<(), DWFAError> {
        if self.is_finalized {
            return Err(DWFAError::AlreadyFinalized);
        }

        // the grid boundaries already hold the end-to-end distance in their corner
        if let Some(tracker) = self.boundary.as_mut() {
            tracker.extend(baseline_seq, other_seq)
                .map_err(DWFAError::Boundary)?;
            let global_distance = tracker.corner();
            if global_distance > self.max_edit_distance {
                self.edit_distance = self.max_edit_distance + 1;
                return Err(DWFAError::MaxEditDistance);
            }
            self.edit_distance = global_distance;
            self.is_finalized = true;
            return Ok(());
        }

        // always try to extend before finalizing
        self.extend(baseline_seq, other_seq)?;

        // a length difference above the cap cannot be aligned, so stop with the same error the wavefront would reach
        let lower_bound = baseline_seq.len().abs_diff(other_seq.len());
        if lower_bound > self.max_edit_distance {
            self.edit_distance = self.max_edit_distance + 1;
            return Err(DWFAError::MaxEditDistance);
        }

        // a small cap aborts the wavefront before a bit-parallel scan of the whole sequences would finish
        let cap_is_cheaper = self.max_edit_distance != usize::MAX
            && !wavefront_exceeds_bitparallel(self.max_edit_distance, baseline_seq, other_seq);

        // keep increasing the edit distance until we reach the full diagonal or we hit the WFA cap
        while !self.reached_full_diagonal(baseline_seq, other_seq) {
            let next = self.edit_distance.saturating_add(1);
            if !cap_is_cheaper && wavefront_exceeds_bitparallel(next, baseline_seq, other_seq) {
                // break this loop if it becomes too expensive to grow the wavefront
                break;
            }
            self.increase_edit_distance(baseline_seq, other_seq)?;
        }

        // if we reached the full diagonal, we are done
        if !self.reached_full_diagonal(baseline_seq, other_seq) {
            // the wavefront stopped early; the end-to-end distance does not depend on it
            let global_distance = crate::util::sequence_alignment::bitparallel_edit_distance(baseline_seq, other_seq);
            debug_assert!(global_distance >= self.edit_distance);

            // if we hit the WFA cap, we are done
            if global_distance > self.max_edit_distance {
                self.edit_distance = self.max_edit_distance + 1;
                return Err(DWFAError::MaxEditDistance);
            }

            // otherwise, we set the edit distance to the global distance and mark it as finalized
            self.edit_distance = global_distance;
        }

        // mark it as finalized
        self.is_finalized = true;
        Ok(())
    }

    /// Helper function that will determine the farthest distance reached into the `baseline_seq` so far.
    pub fn maximum_baseline_distance(&self) -> usize {
        self.max_baseline_distance
    }

    /// Helper function that will determine the farthest distance reached into the `other_seq` so far.
    /// After `update()`, this is the `other_seq` length unless the grid boundaries took over from the wavefront.
    /// After `finalize()`, this is the `other_seq` length only when the wavefront was driven to the end.
    pub fn maximum_other_distance(&self) -> usize {
        self.max_other_distance
    }

    /// Returns true if at least one wavefront is touching the end of the baseline sequence
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    pub fn reached_baseline_end(&self, baseline_seq: &[u8]) -> bool {
        // under "normal", the == will work, but >= is required at finalizing
        self.maximum_baseline_distance() >= baseline_seq.len()
    }

    /// Returns true if at least one wavefront is touching the end of the other sequence
    /// # Arguments
    /// * `other_seq` - the other sequence
    pub fn reached_other_end(&self, other_seq: &[u8]) -> bool {
        // under "normal", the == will work, but >= is required at finalizing
        self.maximum_other_distance() >= other_seq.len()
    }

    /// Returns true if at least one wavefront is at the final diagonal position; i.e. both sequences are fully aligned end-to-end
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    pub fn reached_full_diagonal(&self, baseline_seq: &[u8], other_seq: &[u8]) -> bool {
        self.wavefront.iter().enumerate()
            .any(|(i, &d)| {
                // we want to check if ANY of the wavefronts are at the end of both diagonals
                let base_dist = d + self.edit_distance - i;
                let other_dist = d;
                base_dist >= baseline_seq.len() && other_dist >= other_seq.len()
            })
    }

    /// This will return the set of candidate extensions that do not require increasing the edit distance.
    /// It also includes how many times that character was counted in the event of multiple possible extension points.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence, typically getting updates
    pub fn get_extension_candidates(&self, baseline_seq: &[u8], other_seq: &[u8]) -> HashMap<u8, usize> {
        let mut ret: HashMap<u8, usize> = Default::default();
        for (i, &d) in self.wavefront.iter().enumerate() {
            let other_offset = d;
            if other_offset == other_seq.len() {
                let offset = d + self.edit_distance - i;
                if offset < baseline_seq.len() {
                    // ret.insert(baseline_seq[offset]);
                    let entry = ret.entry(baseline_seq[offset]).or_insert(0);
                    *entry += 1;
                }
            }
        }
        ret
    }

    // Getters below
    pub fn edit_distance(&self) -> usize {
        self.edit_distance
    }

    pub fn wavefront(&self) -> &[usize] {
        &self.wavefront
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cached maxima match a full wavefront scan after extension and edit-distance growth.
    #[test]
    fn test_cached_wavefront_maxima() {
        let baseline = b"ACGTACGTACGT";
        let other = b"ACGTACCTACGT";
        let mut dwfa = DWFALite::default();
        dwfa.update(baseline, other).unwrap();
        dwfa.finalize(baseline, other).unwrap();

        let scanned_baseline = dwfa.wavefront().iter().enumerate()
            .map(|(i, &d)| d + dwfa.edit_distance() - i)
            .max()
            .unwrap();
        let scanned_other = *dwfa.wavefront().iter().max().unwrap();
        assert_eq!(dwfa.maximum_baseline_distance(), scanned_baseline);
        assert_eq!(dwfa.maximum_other_distance(), scanned_other);
    }

    /// Growing the wavefront in place keeps one cell per diagonal.
    #[test]
    fn test_in_place_wavefront_length() {
        let baseline = b"AAAA";
        let other = b"CCCC";
        let mut dwfa = DWFALite::default();
        dwfa.update(baseline, other).unwrap();
        dwfa.finalize(baseline, other).unwrap();
        let distance = dwfa.edit_distance();
        assert_eq!(distance, baseline.len());
        assert_eq!(dwfa.wavefront().len(), distance * 2 + 1);
    }

    #[test]
    fn test_new() {
        let dwfa = DWFALite::default();
        assert_eq!(dwfa.edit_distance(), 0);
        assert_eq!(dwfa.wavefront(), &[0]);
    }

    /// An empty sequence has an edit distance equal to the other sequence's length.
    #[test]
    fn test_empty() {
        let mut dwfa = DWFALite::default();
        let sequence = b"";
        let other = b"ACGT";
        dwfa.finalize(sequence, other).unwrap();
        assert_eq!(dwfa.edit_distance(), 4);

        let mut dwfa = DWFALite::default();
        let sequence = b"ACGT";
        let other = b"";
        dwfa.finalize(sequence, other).unwrap();
        assert_eq!(dwfa.edit_distance(), 4);
    }

    #[test]
    fn test_exact_match() {
        let sequence = b"ACGTACGTACGT";
        let mut other_seq = vec![];
        let mut dwfa = DWFALite::default();
        for &c in sequence.iter() {
            other_seq.push(c);
            assert_eq!(dwfa.update(sequence, &other_seq).unwrap(), 0);
        }
    }

    #[test]
    fn test_simple_mismatch() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACCTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_simple_insertion() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACIGTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_simple_deletion() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_complex_001() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACTACGCACGGGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        dwfa.finalize(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), 4);
    }

    #[test]
    fn test_complex_002() {
        //modified_seq has 2 separate deletions, 1 2bp insertion, and 1 mismatch
        let sequence     = b"AACGGATCAAGCTTACCAGTATTTACGT";
        let alt_sequence = b"AACGGACAAAAGCTTACCTGTATTACGT";
        
        let mut dwfa = DWFALite::default();
        /*
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        */
        // test a shortcut
        dwfa.update(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), 5);
    }

    #[test]
    fn test_big_insertion() {
        // one big insertion in the middle
        //let sequence     = b"AACGGATTTTACGT";
        //let alt_sequence = b"AACGGATAAAAGCTTACCTGTTTTACGT";
        let sequence     = b"AA";
        let alt_sequence = b"ATA";
        
        let mut dwfa = DWFALite::default();
        dwfa.finalize(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), alt_sequence.len() - sequence.len());
    }

    #[test]
    fn test_big_deletion() {
        // one big deletion in the middle
        let sequence     = b"ATTTTTTTTTTAAAAAAAAAA";
        let alt_sequence = b"AAAAAAAAAAA";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), sequence.len() - alt_sequence.len());
    }

    #[test]
    fn test_required_finalize() {
        // the ALT is a lot smaller than the baseline, so we need to run a finalize to get the full ED
        let sequence     = b"ATTTTTTTTTTA";
        let alt_sequence = b"AA";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }

        // here it has only compared AT to AA, so ED=1
        assert_eq!(dwfa.edit_distance(), 1);

        // after finalizing, it will do end-to-end comparison
        dwfa.finalize(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), sequence.len() - alt_sequence.len());
    }

    #[test]
    fn test_cloning() {
        let sequence     = b"AAAAAAA";
        let alt_sequence = b"AAACAAA";

        let mut dwfa = DWFALite::default();
        let mut dwfa2 = dwfa.clone();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &sequence[..(l+1)]).unwrap();
            dwfa2.update(sequence, &alt_sequence[..(l+1)]).unwrap();

            if sequence[l] == alt_sequence[l] {
                // same sequence still
                assert_eq!(dwfa, dwfa2);
            } else {
                // should have different sequences now
                assert_ne!(dwfa, dwfa2);

                // re-clone the first one
                dwfa2 = dwfa.clone();
            }
        }

        // in the end, both should exactly match due to cloning
        assert_eq!(dwfa.edit_distance(), 0);
        assert_eq!(dwfa2.edit_distance(), 0);
    }

    #[test]
    fn test_big_early_termination() {
        let c1 =     "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACACGGGATTATGGGCAAGGGTAACTGACATCTGCTCAGCCTCAACGTACCCGTCTCAAATGCGGCCAGGCGGTGGGGTAAGCAGGAATGAGGCAGGGGTGGGGTTGCCCTGAGGAGGATGATCCCAACGAGGGCGTGAGCAGGGGACCCGAGTTGGAACTACCACATTGCTTTATTGTACATTAGAGCCTCTGGCTAGGGAGCAGGCTGGGGACTAGGTACCCCATTCTAGCGGGGCACAGCACAAAGCTCATAGGGGGATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCTGGGCCGGGGCTGTCCGGCGGCCACGGAGAAGCTGAAGTGCTGCAGCAGGGAGGTGAAGAAGAGGAAGAGCTCCATGCGGGCCAGGGGCTCCCCGAGGCATGCACGGCGGCCTGTGGGGAGGGGAGGGGCGTCAGTGAGCCTGGCTCCTGGGTGATACCCCTGCAAGACTCCACGGAAGGGGACAGGGAGCCGGGCTCCCCACAGGCACCTGCTGAGAAAGGCAGGAAGGCCTCCGGCTTCACAAAGTGGCCCTGGGCATCCAGGAAGTGTTCGGGGTGGAAGCGGAAGGGCTTCTCCCAGACGGCCTCATCCTTCAGCACCGATGACAGGTTGGTGATGAGTGTCGTTCCCTGGGCAGGAGATGCAGGGTGAGAGTGGGGACTGGACTCTAGGATGCTGGGACCCCTGCCACCAAACACACGGGGGACACACACTGCCTGGCACACAGCTGGACTCTGTCAACTAGTCCTGCGCCCGAGAAGCTCCACAGTACCCTCTCCGACCCCACAGCAGGGCGCAGTCACACCTCTCAGAGGCACCCACACTGCCCCCTCTCCCTGCAGGCGCTGGGTCCTCCAACATTCTGGCAGGTCCTGGTTTGTCTCCCCACTAGACGGGGGCTCTGGATGGACAGGCCAGCCCTGCCTATACTCTGGACCCCCCACCCAAGTGGGGACAGTCAGTGTGGTGGCATTGAGGACTAGGTGGCCAGGGTTCCTAGAGTGGGCCCACCTGGCAGTAGCCATGCTGGGGCTATCACCAGGGGCTGGTGCTGAGCTGGGGTGAGGAGGGCGCCAGGCCTACCTTAGGGATGCGGAAGCCCTGTACTTCGATGTCACGGGATGTCATATGGGTCACACCCAGGGGGACGATGTCCCCAAAGCGCTGCACCTCATGAATCACGGCAGTGGTGTAGGGCATGTGAGCCTGGTCACCCATCTCTGGTCGCCGCACCTGCCCTATCACGTCGTCGATCTCCTGTTGGACACGGACTGGACAGACATGCGTCCCCACAATGGGTCAGCACCCAGGGGACACTCTCCTTCCTCCTGTGTTGGAGGAAGTTAGGCTTACAGGAGCCTGGCCACGCCTGTGCTGGAAGCCCCGGGTGTCCCAGCTAAGCCCAGGGGCCCCCAGCTGTACCCTTCCTCCCTCAGTCCCTGCCTTGGGCCCCAGCTGGGCTCACGCTGCACATCCAGGTGTAGGATCATGAGCAGGAGGCCCCAGGCCAGCGTGGTCAAGGTGGTCACCATCCCGGCAAGGAACAGGTTACCCACCACTATGCGCAGGTTCTCATCATTGAAGCTGCTCTCAGGGCTCCCCTTGGCCTGAGCAGGGCCGAGAGGATACTCAGGGGATAGAACGGGGTAGCCCCCAAATGACCTCCAATTCTGCACCTGTCAGCCCAGATGCGGCTCGCCGGGTGATGCACTGGTCCAACCTTTTGCCCAGCCTCCCCTCATTCCTCCTGGGACGTTCAACCCACCACCCTTGCCCCCCACCGTGGCAGCCACTCTCACCTTCTCCTTCTTTGCCAGGAAGGCCTCAGTCAGGTCTCGGGGTGGCTGGGCTGGGTCCCAGGTCATCCTGTGCTCAGTTAGCAGCTCATCCAGCTGGGTCAGGAAAGCCTTTTGGAAGCGTAGGACCTTGCCAGCCAGCGCTGGGATGTGCGGGAGGACGGGGACAGCATTCAGCACCTACACCAGACAGAACCGGGTCTCAATCCTTCCTGTGCTCTGCGTTCATCTGGACCAGTCTCAGGCCCCAGCCATCTCCAGGAAGACCCAGGGCCTGCCTGTCCTTACCACTGACCTCACCAAGTCCCTCCCCAAGTGCCAGCCTCCACCCTCTCTCTCCTTGCCCAGAGGAGAAACCTAAAATCGAAATCTCCAACGTGGACGGGGGTACAGAGTCCTTGGCCTCTCCTGGTGCCCCCTGACCCGGGCACACCTCTCCCACGACCATGTCTGAGATGTCCCCTCCTCCTCCAGGCCCTTCTTACAGTGGGGTCTCCTGGAATGTCCTTTCCCAAACCCATCTACGCAAATCCTGCCCTTCGGAGGCCCCAGTCCAGCCCCGGCACCTCTCAGGAGCTCGCCCTGCAAAGACCCTTGCTCCGCACCTCGCGCAGGAAGCCCGACTCCTCCTTCGATCCCTCCCTGAGCTAGGTCCAGCAGCCTGAGGAAGCGAGGGTCGTCGTACTCGAAGCGGCGCCCGCAGGTGAGGGAGGCGATCACGTTGCTCACGGCTTTGTCCAAGAGACCGTTGGGGCGAAAGGGGCGTCCTGGGGGTGGGAGATGCGGGTAAGGGGTCGCCTTCTCCGTCCCCCGCCTTCCCAGTTCCCGCTGTGTGCCCTTCTGCCCATCACCCACCGGCTTGGTCGGCGAAGGCGGCACAAAGGCAGGCGGCCTCCTCGGTCACCCACTGCTCCAGCGACTTCTTGCCCAGGCCCAAGTTGCGCAAGGTGGACACGGAGAAGCGCCTCTGCTCGCGCCACGCGGGCCCATAGCGCGACAGGATCACCCCTGTGGGCGGGACGGACACGTGGGCGTTGCCATGAAGGCCTTGGCCCCACCCTCCGCCACCCACTCCAACCCTGGCGCTCCACAAGGTCTCCCGCAGTCCCTAGCCCGGTCCAGCTGGGCACAGGGCCCACTCTTTGCTCACCCACATTGCTCCCCTGCCTGGGGCGGGGTTTGGCCCCACCTCGTCTCTGCCCACCCTGACCACCTTTCCACTCAAGGAAGATCCCGCCCGTCCCGCCCACACTGAGCCCGCAGCATAGGCGCGGTCCCCGCCACCGCCACTTCGACGCATCAGCCTCGCCCACCGGGCTTCTGGCGGGTCTGGGCAGTAGCCCCGCCCCCTCCCAGCCCACAGACTCGCACCTCCCCCGTGCAGGTGGTTTCCTGGCCCACTGTCCTCAGCCCACTCGCTGGCCTTTATCTCTGTTTCACGTCCAGGACCCCACGCCCTGTCGGCGCTGCTTGGGCTACGGTCACTGTCCACCCGGGGCCCACGGAAACGCGGTCTCTGTCCCCCACCGCCGCTTGCCTTGGGAACGCGGCCCGAAGCCCAGGACCTGGTAGATGGGCGCAGGCGGGCGGTCGGCCGTGTCCTCGCCGCGGGTCACCATCGCCTCGCGCACGGCCGCCAGCCCATTGAGCACGACCACCGGCGTCCAGGCCAGCTGCAGGCTGAACACGTCCCCGAAGCGGCGCCGCAACTGCAGAGGGAGGGTCAGGGCCTCTTGTCAAGCCAGGATCCCCCCAGACTACAGGTCCTAGTCCTATTTGAACCTTGGACGACCCCCGGGGCTACCAGGAGTGAGCAGGTGGAAGGAGGAGACCCAGCCTCCTGATCCTGGGGCGGGGGTGGGGGTCACACCTTCTGTGATGGAGGAACTCAGTTTGGATGCGTCACCCAGGTATGACCTTGCAAGAGTCACCAAAATTGCCGAGAGGCCCCAGTTAGCATCCCATTCCCAGATGATGGTCCATGCCGGTGAGCAGTGAGGCCCGAGGACCCACAGTGCAAAAGGTTTGAACCGGGTCACTGCACCCCCTTCATCCTCGATTTCGTGATTTAAACGGCACTCAGGACTAACTCATCTTCCATTCCCAAGGCCTTTCCTTCTGGTGTCAGCAGAAGGGACTTTGTACTCCATAACATATGTTGCCCAATGGGCTTGCATGCCCACTGCCAAGTCCAGCTCCACCTCCAGGCCCTTGCCCTACTCTTCCTTGGCCTTTGGAAAATCCAGTCCTTCATGCCATGTATAAATGTCCTTCCCCAGGACGTCCCCCAAACCTGCTTCCCCTTCTCAGCCTGGCTTCTGATCCAGCCTGTGGTTTAACCCACCACCCATGTTTGCTGGTGGTGGGGCATCCTCAGGACCTCTGCCGCCCTCCAGGACCTCCTCCCTCACCTGGTCGAAGCAGTATGGTGTGTTCTGGAAGTCCACATGCAGCAAGGTTGCCCAGCCCGGGCAGTGGCAGGGGACCTGGCGGGTAGCGTGCAGCCCAGCGTTGGTGCCGGTGCATCAGGTCCACCAGGAGCAGGAAGATGGCCACTATCATGGCCAGGGGCACCAGTGCTTCTAGCCCCATGGCTGCCTCACTACCAACTGGGCTCCTCTGGACACACCTGGCACCCCCACCCCACCAGGCACAGAGGACCAGGCAGGACACTCTCAGCACACCGAGCGCGTGACCCTTCCCTTATAAAGGGAGCTGATGATGGCCTTCGCCCTCTGCTGTGAGTGAACCTGCTGTGTTGACTGTGCTGCCAGTGGCAGAGTCAGGCCAGGGTGGGTATGGGCTGCTCCAGAGGTCCTTGCCGCTGCTTCCTGCTCCAGGCCCTTACCCAGGGTAGGGTGGTAGAAAGGCCTGGTCGGAGAAGTCACCCCCTCTCCCCACTCCAAGCTCCCCAAGCCCACACAGGCTTCTGGGATAACCAGGGTCTCAGTGGACCCGGCCATCCACCTCCCAGCTAGGCTCATACACCGTAATGTAGTCACAACCCCTCCTCCAGAACATGGCCTTGCCCTTTCCCTACCCCCACCTGCCCACTCCAGAGTGACCTTCAGCACCCTTATCTGTCACTGGCACTTACCTGGGGCCTTAGAGCTCCTGATGATGAGTGGCATCATGGGCCTGGTCCCTTCACTTCACCTTGCACTCTTGACATGCACAGACGCTATGCACACACCTGATGGTGCACAGATCTCTTGTCCACTCCCAGACACTTGTCCACTTGTTCACACTTGCAGGGACACGATTACACATGCAGAAAATCACCCACACAAAGACAATATTCACACATACACAGACTCACACTGACACTCAGGGCACACATTCTCTCTCACACACACCAGTCACACACACATACAGACCCGGCACCAAGTACCCCACTTCCCAGCCATGCCCAAGGTTTCCTGGATGGGACCTCTCCTGTCCAGAGGCTGCTCCCAGTGAGCCTCAAAGCTGTCACGTGGATCCCAGCTCAGCCCACATTCTGGGCTCTGGCCGGGCCATGGCTTCCTGTTTGCAACAGGGCTGTTCCCAGAGCTCCCAGTTGGTAGCCTGAAGGCCCTTGCCCCAGCCTGTGACAGCATCCTCCAGGGCTGCCTGAGGGTCGTCATTCTCCACTGCTTCCTGGCCTCCATGTTTCTGATTAGAAATCTGGTGGAAACATTATGGAGGATCCTTTATTTAGGATATGTTGCTTTTTTATTTTTATTTTTTCTTTAGACAGGGTCTCACTCTGTTGCCCGGGCCGGAGTGCAGTGGCAGGATCACGGCTCACTGCAATCTCAACATCAAGTGGACCTCCTGCCTCCCAAGTAGCTGGGACTACAGGCACCACCGAGCCCAAATAATTTTTTTTTTGAGACGGAGTTTTGCTCTGTCGCCCAGGTGGGAGTGCAATGATGCGATCTCGGCTCACTGCAACCTCCACCTCCAGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTACAGGTGCCCACCACCATGCCTGGCTGATTTTTTGTA";
        let seq_23 = "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAACTGACATCTGCTCAGCCTCAACGTACCCGTCTCAAATGCGGCCAGGCGGTGGGGTAAGCAGGAATGAGGCAGGGGTGGGGTTGCCCTGAGGAGGATGATCCCAACGAGGGCGTGAGCAGGGGACCCGAGTTGGAACTACCACATTGCTTTATTGTACATTAGAGCCTCTGGCTAGGGAGCAGGCTGGGGACTAGGTACCCCATTCTAGCGGGGCACAGCACAAAGCTCGTAGGGGGATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCTGGGCCGGGGCTGTCCGGCGGCCACGGAGAAGCTGAAGTGCTGCAGCAGGGAGGTGAAGAAGAGGAAGAGCTCCATGCGGGCCAGGGGCTCCCCGAGGCATGCACGGCGGCCTGTGGGGAGGGGAGGGGCGTCAGTGAGCCTGGCTCCTGGGTGATACCCCTGCAAGACTCCACGGAAGGGGACAGGGAGCCGGGCTCCCCACAGGCACCTGCTGAGAAAGGCAGGAAGGCCTCCGGCTTCACAAAGTGGCCCTGGGCATCCAGGAAGTGT";

        // iterate, making sure everything is fine even as we go well beyond seq_23
        let mut dwfa = DWFALite::default();
        for i in 0..c1.len() {
            dwfa.update(seq_23.as_bytes(), &c1.as_bytes()[0..(i+1)]).unwrap();
            assert!(dwfa.edit_distance() <= 2);
        }
        assert_eq!(dwfa.edit_distance(), 2);

        // when we finalize, now we detect the big ED
        dwfa.finalize(seq_23.as_bytes(), c1.as_bytes()).unwrap();
        assert_eq!(dwfa.edit_distance(), 5278);
    }

    /// A large edit distance is computed directly and still matches the full grid, without driving the wavefront to the end.
    #[test]
    fn test_finalize_large_distance_matches_full_grid() {
        let baseline = vec![b'A'; 300];
        let other = vec![b'C'; 300];
        let mut dwfa = DWFALite::default();
        dwfa.finalize(&baseline, &other).unwrap();
        assert_eq!(dwfa.edit_distance(), crate::util::sequence_alignment::edit_distance(&baseline, &other));
        assert!(dwfa.wavefront().len() < dwfa.edit_distance() * 2 + 1);
        assert!(dwfa.is_finalized);
    }

    /// A large distance hands over to the grid boundaries during `update()`, with the same distance the wavefront reaches.
    /// Later updates and `finalize()` stay exact as the sequences grow.
    #[test]
    fn test_update_uses_boundary_for_large_distance() {
        let baseline = vec![b'A'; 3_000];
        let mut other = vec![b'C'; 3_000];
        let mut dwfa = DWFALite::default();
        assert_eq!(dwfa.update(&baseline, &other).unwrap(), 3_000);
        assert!(dwfa.boundary.is_some());

        // the other sequence grows, which the wavefront would also have to follow
        other.extend_from_slice(&[b'C'; 100]);
        let mut wavefront_only = DWFALite::default();
        wavefront_only.update(&baseline, &other[..3_000]).unwrap();
        assert_eq!(dwfa.update(&baseline, &other).unwrap(), wavefront_only.update(&baseline, &other).unwrap());

        dwfa.finalize(&baseline, &other).unwrap();
        assert_eq!(dwfa.edit_distance(), crate::util::sequence_alignment::edit_distance(&baseline, &other));
        assert!(dwfa.is_finalized);
    }

    /// Distances below the floor stay on the wavefront, even when the bit-parallel scan would be cheaper.
    #[test]
    fn test_update_keeps_wavefront_for_small_distance() {
        let baseline = vec![b'A'; MIN_BOUNDARY_EDIT_DISTANCE];
        let other = vec![b'C'; MIN_BOUNDARY_EDIT_DISTANCE];
        let mut dwfa = DWFALite::default();
        assert_eq!(dwfa.update(&baseline, &other).unwrap(), MIN_BOUNDARY_EDIT_DISTANCE);
        assert!(dwfa.boundary.is_none());
    }

    /// Grows `sequence` by `count` symbols that follow `template`, with occasional deletions, insertions, and substitutions.
    /// # Arguments
    /// * `sequence` - the sequence to append to
    /// * `template` - the sequence being followed, repeated if needed
    /// * `count` - how many template positions to consume
    /// * `next` - source of pseudo-random numbers
    fn grow_like_template(sequence: &mut Vec<u8>, template: &[u8], count: usize, next: &mut impl FnMut() -> u64) {
        for _ in 0..count {
            let base = template[sequence.len() % template.len()];
            match next() % 25 {
                0 => {},
                1 => sequence.push(b"ACGT"[(next() % 4) as usize]),
                2 => { sequence.push(base); sequence.push(base); },
                _ => sequence.push(base),
            }
        }
    }

    /// Once the grid boundaries take over, every `update()` and the final `finalize()` match the pure wavefront.
    /// Both sequences grow in random chunks, and the hand-over happens at a different step in each case.
    #[test]
    fn test_boundary_matches_wavefront_while_growing() {
        let mut state: u64 = 0xA0761D6478BD642F;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        for case in 0..300 {
            // two sequences that mostly agree in some cases, and are unrelated in others
            let template: Vec<u8> = (0..700).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let unrelated: Vec<u8> = (0..700).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let other_template = if case % 3 == 0 { &unrelated } else { &template };

            let mut baseline: Vec<u8> = vec![];
            let mut other: Vec<u8> = vec![];
            let mut wavefront_only = DWFALite::default();
            let mut subject = DWFALite::default();
            let switch_step = (next() % 8) as usize;
            for step in 0..8 {
                // chunk sizes straddle the 64-symbol block size, and either side may not grow
                let chunk = |value: u64| match value % 6 {
                    0 => 0,
                    1 => 1,
                    2 => 63,
                    3 => 64,
                    4 => 65,
                    _ => (value / 6 % 150) as usize,
                };
                let baseline_chunk = chunk(next());
                let other_chunk = chunk(next());
                grow_like_template(&mut baseline, &template, baseline_chunk, &mut next);
                grow_like_template(&mut other, other_template, other_chunk, &mut next);

                let expected = wavefront_only.update(&baseline, &other).unwrap();
                let actual = if step == switch_step {
                    subject.enter_boundary(&baseline, &other).unwrap()
                } else {
                    subject.update(&baseline, &other).unwrap()
                };
                assert_eq!(actual, expected, "case {case} step {step} lengths {} {}", baseline.len(), other.len());
                assert_eq!(subject.edit_distance(), wavefront_only.edit_distance(), "case {case} step {step}");
            }

            wavefront_only.finalize(&baseline, &other).unwrap();
            subject.finalize(&baseline, &other).unwrap();
            assert!(subject.boundary.is_some());
            assert_eq!(subject.edit_distance(), wavefront_only.edit_distance(), "case {case} final");
            assert_eq!(subject.edit_distance(), crate::util::sequence_alignment::edit_distance(&baseline, &other), "case {case} final grid");
        }
    }

    /// A cap is reported the same way by the grid boundaries as by the wavefront, with the distance left at cap + 1.
    #[test]
    fn test_boundary_respects_max_edit_distance() {
        let mut state: u64 = 0xE7037ED1A0B428DB;
        let mut next = || { state ^= state << 13; state ^= state >> 7; state ^= state << 17; state };
        for case in 0..200 {
            let cap = (next() % 40) as usize;
            let baseline: Vec<u8> = (0..(next() % 150) as usize).map(|_| b"ACGT"[(next() % 4) as usize]).collect();
            let other: Vec<u8> = (0..(next() % 150) as usize).map(|_| b"ACGT"[(next() % 4) as usize]).collect();

            // update, then finalize if the update was allowed
            let mut wavefront_only = DWFALite::with_max_edit_distance(cap);
            let mut subject = DWFALite::with_max_edit_distance(cap);
            let expected = wavefront_only.update(&baseline, &other);
            let actual = subject.enter_boundary(&baseline, &other);
            assert_eq!(actual.is_ok(), expected.is_ok(), "case {case} cap {cap} update");
            assert_eq!(subject.edit_distance(), wavefront_only.edit_distance(), "case {case} cap {cap} update distance");
            if expected.is_err() {
                continue;
            }
            let expected = wavefront_only.finalize(&baseline, &other);
            let actual = subject.finalize(&baseline, &other);
            assert_eq!(actual.is_ok(), expected.is_ok(), "case {case} cap {cap} finalize");
            assert_eq!(subject.edit_distance(), wavefront_only.edit_distance(), "case {case} cap {cap} finalize distance");
        }
    }

    /// A small cap reports `MaxEditDistance` with the distance left at cap + 1, and the calculation is not finalized.
    #[test]
    fn test_finalize_respects_max_edit_distance() {
        // different lengths: the lower bound rejects before either algorithm runs
        let mut dwfa = DWFALite::with_max_edit_distance(0);
        let error = dwfa.finalize(b"AAAA", b"AA").unwrap_err();
        assert!(matches!(error, DWFAError::MaxEditDistance));
        assert_eq!(dwfa.edit_distance(), 1);
        assert!(!dwfa.is_finalized);

        // same length, so the distance itself is above the cap
        let baseline = [b'A'; 200];
        let other = [b'C'; 200];
        let mut dwfa = DWFALite::with_max_edit_distance(0);
        let error = dwfa.finalize(&baseline, &other).unwrap_err();
        assert!(matches!(error, DWFAError::MaxEditDistance));
        assert_eq!(dwfa.edit_distance(), 1);
        assert!(!dwfa.is_finalized);

        // the exact distance is within the cap
        let mut dwfa = DWFALite::with_max_edit_distance(4);
        dwfa.finalize(b"AAAA", b"CCCC").unwrap();
        assert_eq!(dwfa.edit_distance(), 4);
        assert!(dwfa.is_finalized);
    }

    /// Extends from a start offset one byte at a time, as the original loop did.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    /// * `baseline_offset` - where the extension starts in the baseline
    /// * `other_offset` - where the extension starts in the other sequence
    /// # Returns
    /// The `(baseline, other)` offsets where the extension stops.
    fn extend_bytewise(baseline_seq: &[u8], other_seq: &[u8], mut baseline_offset: usize, mut other_offset: usize) -> (usize, usize) {
        if baseline_offset < baseline_seq.len() && other_offset < other_seq.len() {
            while baseline_offset < baseline_seq.len() &&
                other_offset < other_seq.len() &&
                baseline_seq[baseline_offset] == other_seq[other_offset] {
                baseline_offset += 1;
                other_offset += 1;
            }
        }
        (baseline_offset, other_offset)
    }

    /// The word-wise extension stops at the same offsets as a byte-wise loop for every mismatch position around the 8-byte word borders.
    #[test]
    fn test_extend_diagonal_matches_bytewise() {
        // sequence lengths straddle the word size so the tail loop and the word loop both run
        for length in [0usize, 1, 7, 8, 9, 15, 16, 17, 31, 32, 33, 70] {
            let baseline: Vec<u8> = (0..length).map(|i| b"ACGT"[(i * 7 + i / 3) % 4]).collect();

            // no mismatch, and a mismatch at every position
            let mut others: Vec<Vec<u8>> = vec![baseline.clone()];
            for pos in 0..length {
                let mut other = baseline.clone();
                other[pos] = if other[pos] == b'A' { b'C' } else { b'A' };
                others.push(other);
            }

            // a shorter and a longer other sequence
            if length > 0 {
                others.push(baseline[..length - 1].to_vec());
            }
            let mut longer = baseline.clone();
            longer.extend_from_slice(b"ACGTACGTACGT");
            others.push(longer);

            for other in others.iter() {
                // each wavefront index at edit distance 2 shifts the baseline offset differently
                for start in 2..=5usize {
                    let mut dwfa = DWFALite { edit_distance: 2, wavefront: vec![start; 5], ..Default::default() };
                    for index in 0..5 {
                        let baseline_start = start + 2 - index;
                        let (expected_baseline, expected_other) = extend_bytewise(&baseline, other, baseline_start, start);
                        let before_baseline = dwfa.max_baseline_distance;
                        let before_other = dwfa.max_other_distance;
                        dwfa.extend_diagonal(index, &baseline, other);
                        assert_eq!(dwfa.wavefront[index], expected_other, "length {length} start {start} index {index}");
                        assert_eq!(dwfa.max_baseline_distance, before_baseline.max(expected_baseline), "length {length} start {start} index {index}");
                        assert_eq!(dwfa.max_other_distance, before_other.max(expected_other), "length {length} start {start} index {index}");
                    }
                }
            }
        }
    }
}
