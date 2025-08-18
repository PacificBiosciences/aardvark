
use crate::data_types::phase_enums::Allele;
use crate::data_types::variants::Variant;
use crate::dwfa::dynamic_wfa::{DWFAError, DWFALite};
use crate::util::sequence_alignment::edit_distance;

#[derive(thiserror::Error, Debug)]
pub enum HapDWFAError {
    #[error("error from DWFA: {error}")]
    DError { error: DWFAError },
    #[error("cannot extend Allele::Unknown")]
    UnknownAllele,
}

/// Haplotype-specific details. This node does not track the comparators for efficiency, so they must be provided with each call.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct HaplotypeDWFA {
    /// Tracks the constructed truth haplotype
    truth_haplotype: HaplotypeTracker,
    /// Tracks the constructed query haplotype
    query_haplotype: HaplotypeTracker,
    /// The DWFA that tracks our edit distance to the baseline/comparator sequence
    dwfa: DWFALite
}

impl HaplotypeDWFA {
    /// Constructor
    /// # Arguments
    /// * `region_start` - reference start of this DWFA
    pub fn new(region_start: usize, max_edit_distance: usize) -> Self {
        Self {
            truth_haplotype: HaplotypeTracker::new(region_start),
            query_haplotype: HaplotypeTracker::new(region_start),
            // no wildcards, and we want end-to-end without early termination
            dwfa: DWFALite::with_max_edit_distance(max_edit_distance)
        }
    }

    /// Adds variants sequence to both haplotypes. Returns true if any ALT components were successfully incorporated.
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `is_truth` - if true, add the variant to the truth sequence; otherwise, to the query sequence
    /// * `variant` - the variant we are traversing
    /// * `allele` - the REF/ALT indication
    /// * `sync_extension` - if provided, this will further copy the reference sequence to both truth and query
    pub fn extend_variant(&mut self, reference: &[u8], is_truth: bool, variant: &Variant, allele: Allele, sync_extension: Option<usize>) -> Result<bool, HapDWFAError> {
        // extend the appropriate haplotype
        let success = if is_truth {
            // extend our query by the reference up to the sync point, passing the sync point also
            if let Some(sync_pos) = sync_extension {
                self.query_haplotype.copy_reference(reference, sync_pos)?;
            }
            // now add the variant to truth
            self.truth_haplotype.extend_variant(reference, variant, allele, sync_extension)?
        } else {
            // extend our truth by the reference up to the sync point
            if let Some(sync_pos) = sync_extension {
                self.truth_haplotype.copy_reference(reference, sync_pos)?;
            }
            // now add the variant to query, passing the sync point also
            self.query_haplotype.extend_variant(reference, variant, allele, sync_extension)?
        };

        // now update the DWFA
        self.update_dwfa()?;
        Ok(success)
    }

    /// Updates the contained DWFA after a change to the internally tracked sequence.
    /// # Arguments
    /// * `baseline` - the baseline sequence we compare against
    fn update_dwfa(&mut self) -> Result<(), HapDWFAError> {
        // update each DWFA, translate errors to anyhow
        match self.dwfa.update(self.truth_haplotype.sequence(), self.query_haplotype.sequence()) {
            Ok(_) => Ok(()),
            Err(error) => Err(HapDWFAError::DError { error })
        }
    }

    /// Fills out the rest of the region with reference sequence and then finalizes the DWFA.
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `region_end` - the end position of the region
    pub fn finalize_dwfa(&mut self, reference: &[u8], region_end: usize) -> Result<(), HapDWFAError> {
        // add any remaining reference sequence
        self.truth_haplotype.copy_reference(reference, region_end)?;
        self.query_haplotype.copy_reference(reference, region_end)?;
        self.update_dwfa()?;

        // finalize the DWFA, translate errors to anyhow
        match self.dwfa.finalize(self.truth_haplotype.sequence(), self.query_haplotype.sequence()) {
            Ok(_) => Ok(()),
            Err(error) => Err(HapDWFAError::DError { error })
        }
    }

    /// Returns true if the sequences are identical AND they are waiting at the same reference position.
    /// In this situation, they are "synchronized", in the sense that all following changes are independent of the previous changes.
    pub fn is_synchronized(&self) -> bool {
        /*
        Key insight comes from https://www.biorxiv.org/content/10.1101/023754v2.full.pdf
        Page 3, part (c) of the reducing paths explanation, paraphrased:
        - IF the two paths so far are completely identical (same length also, so not just compatible)
        - AND they are at the same "reference" position
        - THEN any changes happening after this point are 100% independent of what has happened up to this point

        In other words, we *know* this node is on the best path, so all other nodes can be pruned from the search space.
         */
        self.edit_distance() == 0 &&
            self.truth_haplotype.sequence().len() == self.query_haplotype.sequence().len() &&
            self.truth_haplotype.ref_pos() == self.query_haplotype.ref_pos()
    }

    /// Edit distance from the underlying DWFA.
    pub fn edit_distance(&self) -> usize {
        self.dwfa.edit_distance()
    }

    /// Returns the total number of set alleles across both truth and query
    pub fn set_alleles(&self) -> usize {
        self.truth_haplotype.alleles().len() + self.query_haplotype.alleles().len()
    }

    /// Returns the total skip distance of truth + query
    pub fn total_variant_skip_distance(&self) -> usize {
        self.truth_haplotype.variant_skip_distance() + self.query_haplotype.variant_skip_distance()
    }

    /// Returns the total cost of this comparison, which includes both the edit distance from the DWFA and also the cost of ignored, incompatible variants.
    pub fn total_cost(&self) -> usize {
        self.edit_distance() + self.total_variant_skip_distance()
    }

    // getters
    pub fn truth_haplotype(&self) -> &HaplotypeTracker {
        &self.truth_haplotype
    }

    pub fn query_haplotype(&self) -> &HaplotypeTracker {
        &self.query_haplotype
    }
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct HaplotypeTracker {
    /// Current position in the reference sequence
    ref_pos: usize,
    /// The alleles associated with this node
    alleles: Vec<Allele>,
    /// The generated haplotype sequence so far
    sequence: Vec<u8>,
    /// The total ED of skipped variants
    variant_skip_distance: usize
}

impl HaplotypeTracker {
    /// Constructor
    /// # Arguments
    /// * `region_start` - tracks the current reference position
    pub fn new(region_start: usize) -> Self {
        Self {
            ref_pos: region_start,
            alleles: vec![],
            sequence: vec![],
            variant_skip_distance: 0
        }
    }

    /// Adds variants sequence to both haplotypes. Returns true if any ALT components were successfully incorporated.
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `variant` - the variant we are traversing
    /// * `allele` - the REF/ALT indication
    /// * `ref_extension` - if provided, this will further copy the reference sequence after the allele is added through the provided position
    pub fn extend_variant(&mut self, reference: &[u8], variant: &Variant, allele: Allele, ref_extension: Option<usize>) -> Result<bool, HapDWFAError> {
        // no matter what, we can always add reference sequence up to the variant start
        let variant_start = variant.position() as usize;
        self.copy_reference(reference, variant_start)?;

        // extend hap1 by the corresponding allele
        let success = match allele {
            Allele::Unknown => return Err(HapDWFAError::UnknownAllele),
            Allele::Reference => {
                // reference allele, so nothing to add
                true
            },
            Allele::Alternate => {
                // check for compatible variants
                if self.ref_pos <= variant_start {
                    // now add the variant sequence
                    let seq = variant.allele1();
                    self.sequence.extend_from_slice(seq);

                    // update our ref pos for this haplotype
                    self.ref_pos = variant_start + variant.ref_len();
                    true
                } else {
                    // incompatible, so just ignore it; but we need to add the ED as incompatible cost
                    self.variant_skip_distance += edit_distance(variant.allele0(), variant.allele1());
                    false
                }
            }
        };
        self.alleles.push(allele);

        // if provided with more extension power, then do it
        if let Some(p) = ref_extension {
            self.copy_reference(reference, p)?;
        }

        Ok(success)
    }

    /// Copies reference sequence up to `region_end`.
    /// # Arguments
    /// * `reference` - the full reference sequence
    /// * `region_end` - the end position of the region
    pub fn copy_reference(&mut self, reference: &[u8], region_end: usize) -> Result<(), HapDWFAError> {
        // add any remaining reference sequence
        if self.ref_pos < region_end {
            // first update the reference sequence to this point
            let reference_update = &reference[self.ref_pos..region_end];
            self.sequence.extend_from_slice(reference_update);
            self.ref_pos = region_end;
        }
        Ok(())
    }

    // getters
    pub fn ref_pos(&self) -> usize {
        self.ref_pos
    }

    pub fn alleles(&self) -> &[Allele] {
        &self.alleles
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn variant_skip_distance(&self) -> usize {
        self.variant_skip_distance
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_haplotype_node() {
        // simple test just to make sure the extension works at the haplotype level
        let reference = b"ACGTACGTACGT";
        let region_start = 3;
        let region_end = 11;

        // add the variant to both truth and query
        let variant = Variant::new_snv(0, 6, b"G".to_vec(), b"A".to_vec()).unwrap();
        let mut hap_node = HaplotypeDWFA::new(region_start, usize::MAX); // start with first 'T'
        hap_node.extend_variant(
            reference, true,
            &variant, Allele::Alternate, None
        ).unwrap();
        hap_node.extend_variant(
            reference, false,
            &variant, Allele::Alternate, None
        ).unwrap();

        // add an extra variant that is a FP and generates ED
        let variant2 = Variant::new_snv(0, 7, b"T".to_vec(), b"A".to_vec()).unwrap();
        hap_node.extend_variant(
            reference, false,
            &variant2, Allele::Alternate, None
        ).unwrap();
        assert_eq!(hap_node.edit_distance(), 0); // ED won't get generated just yet

        // add a reference variant, so no change
        let variant3 = Variant::new_insertion(0, 8, b"A".to_vec(), b"TTTTTTTT".to_vec()).unwrap();
        hap_node.extend_variant(
            reference, false,
            &variant3, Allele::Reference, None
        ).unwrap();
        assert_eq!(hap_node.edit_distance(), 0);

        // finalize it, cutting off the last 'T'
        hap_node.finalize_dwfa(reference, region_end).unwrap();

        // check all these values
        assert_eq!(hap_node.edit_distance(), 1); // one FP variant that gets detected on finalizing
        assert_eq!(hap_node.truth_haplotype().alleles(), &[Allele::Alternate]);
        assert_eq!(hap_node.query_haplotype().alleles(), &[Allele::Alternate, Allele::Alternate, Allele::Reference]);
        assert_eq!(hap_node.truth_haplotype().sequence(), b"TACATACG"); // G->A
        assert_eq!(hap_node.query_haplotype().sequence(), b"TACAAACG"); // G->A and FP T->A
    }

    #[test]
    fn test_haplotype_node_incompatible() {
        // simple test just to make sure the extension works at the haplotype level
        let reference = b"ACGTACGTACGT";
        let region_start = 0;
        let region_end = reference.len();

        // add the variant to both truth and query
        let variant = Variant::new_snv(0, 6, b"G".to_vec(), b"A".to_vec()).unwrap();
        let mut hap_node = HaplotypeDWFA::new(region_start, usize::MAX); // start with first 'T'
        hap_node.extend_variant(
            reference, true,
            &variant, Allele::Alternate, None
        ).unwrap();
        hap_node.extend_variant(
            reference, false,
            &variant, Allele::Alternate, None
        ).unwrap();

        // add an extra variant that is a FP and incompatible
        let variant2 = Variant::new_snv(0, 6, b"G".to_vec(), b"T".to_vec()).unwrap();
        hap_node.extend_variant(
            reference, false,
            &variant2, Allele::Alternate, None
        ).unwrap();
        assert_eq!(hap_node.edit_distance(), 0); // ED won't get generated
        assert_eq!(hap_node.total_variant_skip_distance(), 1); // but VD does

        // finalize it, cutting off the last 'T'
        hap_node.finalize_dwfa(reference, region_end).unwrap();

        // check all these values
        assert_eq!(hap_node.edit_distance(), 0); // no ED
        assert_eq!(hap_node.total_variant_skip_distance(), 1); // one VD
        assert_eq!(hap_node.total_cost(), 1); // 0 + 1
    }
}