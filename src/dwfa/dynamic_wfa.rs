
use rustc_hash::FxHashMap as HashMap;

#[derive(thiserror::Error, Debug)]
pub enum DWFAError {
    #[error("maximum edit distance exceeded")]
    MaxEditDistance,
    #[error("DWFA is already finalized, it cannot be changed further")]
    AlreadyFinalized,
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
    /// If true, this DWFA has been finalized, meaning it cannot be extended anymore.
    is_finalized: bool,
    /// The maximum allowed edit distance before throwing an error
    max_edit_distance: usize
}

impl Default for DWFALite {
    fn default() -> Self {
        DWFALite {
            edit_distance: 0,
            wavefront: vec![0],
            is_finalized: false,
            max_edit_distance: usize::MAX
        }
    }
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
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If this function is called after `finalize()` has been called.
    pub fn update(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, DWFAError> {
        if self.is_finalized {
            return Err(DWFAError::AlreadyFinalized);
        }

        // maximally extend everything along the current diagonals
        self.extend(baseline_seq, other_seq)?;
        
        // check how it looks
        while !self.reached_baseline_end(baseline_seq) && !self.reached_other_end(other_seq) {
            // increase the edit distance, re-extension happens automatically
            self.increase_edit_distance(baseline_seq, other_seq)?;
        }


        Ok(self.edit_distance)
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
        for (i, d) in self.wavefront.iter_mut().enumerate() {
            // `i` is the index in the wavefront
            // `i // 2` is always the middle diagonal
            // anything less than that has used more baseline symbols
            // anything more than that has used more other symbols
            // anything above or below needs to shift the comparison

            // for this particular wavefront, extend as far as possible
            loop {
                // example at edit distance 1:
                // baseline offset would equal: [d+1, d, d-1] which corresponds to
                // a skipped base in primary, equal usage, and a skipped base in other

                // if we truncate a wavefront (e.g., a bounded size); then `i` below will be too small for every wavefront from the start
                //     that gets truncated; this means we want the below formula to be `*d + self.edit_distance - (i + truncate_prefix)`
                //     where `truncate_prefix` is the number we have cut off the front (in total).
                //     IMPORTANT: this formula will need to get updated everywhere that `baseline_offset` is calculated

                let baseline_offset = *d + self.edit_distance - i;
                let other_offset = *d;
                if baseline_offset >= baseline_seq.len() ||
                    other_offset >= other_seq.len() ||
                    baseline_seq[baseline_offset] != other_seq[other_offset] {
                    // if we are past the end of either sequence OR
                    // the sequences are not equal at this position THEN
                    // we are done extending
                    break;
                }

                // we are not done, so add one to this wavefront
                *d += 1;
            }
        }
        Ok(())
    }

    /// This will increase the edit distance for this DWFA and create a new larger wavefront.
    /// This function will automatically call `extend()` to enforce the assumptions.
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

        // create an empty new wavefront
        let new_wf_len = self.wavefront.len() + 2;
        let mut new_wavefront: Vec<usize> = vec![0; new_wf_len];

        // now we need to populate the wavefront
        for (i, &d) in self.wavefront.iter().enumerate() {
            // deletion (skipping) of a base in `baseline_seq`; this does not change the distance into `other_seq`
            new_wavefront[i] = new_wavefront[i].max(d);

            // mismatch progresses both sequences
            new_wavefront[i+1] = new_wavefront[i+1].max(d + 1);

            // insertion of a base into `baseline_seq`; this progresses `other_seq`
            new_wavefront[i+2] = new_wavefront[i+2].max(d + 1);
        }

        // finally, save the new wavefront
        self.wavefront = new_wavefront;

        // re-extend
        self.extend(baseline_seq, other_seq)?;
        Ok(())
    }

    /// This function is similar to `update(...)`, but it will perform all the final steps for an end-to-end calculation.
    /// This will trigger the algorithm to make sure we have reached the end of both the `baseline_seq` and `other_seq`, potentially increasing edit distance further.
    /// After calling this, no further extensions can be made, it is a "finalized" calculation.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence
    /// * `other_seq` - the other sequence
    /// # Errors
    /// * If the edit distance cannot be increased further and it needs to be.
    pub fn finalize(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<(), DWFAError> {
        if self.is_finalized {
            return Err(DWFAError::AlreadyFinalized);
        }

        // always try to extend before finalizing
        self.extend(baseline_seq, other_seq)?;

        // now check where we are
        while !self.reached_full_diagonal(baseline_seq, other_seq) {
            // increase ED as long as we are not at the end of BOTH sequences
            self.increase_edit_distance(baseline_seq, other_seq)?;
        }
        self.is_finalized = true;
        Ok(())
    }

    /// Helper function that will determine the farthest distance reached into the `baseline_seq` so far.
    pub fn maximum_baseline_distance(&self) -> usize {
        // baseline distance requires some compute
        // the 0-index corresponds to deleting `edit_distance` bases in `baseline`, so it has the largest offset
        // each additional iteration pushes the diagonal closer to inserting bases into `baseline`, so the shift gets progressively smaller
        self.wavefront.iter().enumerate()
            .map(|(i, &d)| d + self.edit_distance - i)
            .max().unwrap()
    }

    /// Helper function that will determine the farthest distance reached into the `other_seq` so far.
    /// After a public function call, this should _always_ be the `other_seq` length.
    pub fn maximum_other_distance(&self) -> usize {
        // other distance is directly tracked in our wavefront
        *self.wavefront.iter().max().unwrap()
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

    #[test]
    fn test_new() {
        let dwfa = DWFALite::default();
        assert_eq!(dwfa.edit_distance(), 0);
        assert_eq!(dwfa.wavefront(), &[0]);
    }

    #[test]
    fn test_empty() {
        let mut dwfa = DWFALite::default();
        let sequence = b"";
        let other = b"ACGT";
        dwfa.finalize(sequence, other).unwrap();
        assert_eq!(dwfa.edit_distance(), 4);
        assert_eq!(dwfa.wavefront(), &[0, 1, 2, 3, 4, 4, 4, 4, 4]);

        let mut dwfa = DWFALite::default();
        let sequence = b"ACGT";
        let other = b"";
        dwfa.finalize(sequence, other).unwrap();
        assert_eq!(dwfa.edit_distance(), 4);
        assert_eq!(dwfa.wavefront(), &[0, 1, 2, 3, 4, 4, 4, 4, 4]);
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
            dwfa.update(seq_23.as_bytes(), c1[0..(i+1)].as_bytes()).unwrap();
            assert!(dwfa.edit_distance() <= 2);
        }
        assert_eq!(dwfa.edit_distance(), 2);

        // when we finalize, now we detect the big ED
        dwfa.finalize(seq_23.as_bytes(), c1.as_bytes()).unwrap();
        assert_eq!(dwfa.edit_distance(), 5278);
    }
}
