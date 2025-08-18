
use crate::util::sequence_alignment::wfa_ed;

/// All the variant types we are currently allowing
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord)]
pub enum VariantType {
    /// REF and ALT are both length = 1
    Snv=0,
    /// REF length = 1, ALT length > 1
    Insertion,
    /// REF length > 1, ALT length = 1
    Deletion,
    /// REF and ALT lengths > 1
    Indel,
    /// Must have two alleles and be tagged with SVTYPE=INS; ALT >= REF
    SvInsertion,
    /// Must have two alleles and be tagged with SVTYPE=DEL; ALT <= REF
    SvDeletion,
    /// Must have two alleles and be tagged with SVTYPE=DUP
    SvDuplication,
    /// Must have two alleles and be tagged with SVTYPE=INV
    SvInversion,
    /// Must have two alleles and be tagged with SVTYPE=BND
    SvBreakend,
    /// Must have two alleles and be tagged with TRID=####, ALT < REF
    TrContraction,
    // Must have two alleles and be tagged with TRID=####, ALT >= REF
    TrExpansion,
    /// Something that doesn't match the above criteria, must be 1 or 2 alleles
    Unknown // make sure Unknown is always the last one in the list
}

/// Zygosity definitions, mostly used elsewhere
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord)]
pub enum Zygosity {
    HomozygousReference=0,
    Heterozygous,
    HomozygousAlternate,
    Unknown // make sure Unknown is always the last one in the list
}

#[derive(thiserror::Error, Debug)]
pub enum VariantError {
    #[error("allele0 length must match reference length when index_allele0=0")]
    Allele0RefLen,
    #[error("allele{index} must be length 1")]
    AlleleLen1{ index: usize },
    #[error("allele{index} is empty (length = 0)")]
    EmptyAllele{ index: usize },
    #[error("index_allele0 must be < index_allele1")]
    IndexAlleleOrder,
    #[error("{variant_type:?} does not support multi-allelic sites")]
    MultiAllelicNotAllowed { variant_type: VariantType },
    #[error("reference must have length > 1")]
    RefLenGT1,
    #[error("alternate must have length > 1")]
    AltLenGT1,
    #[error("SV deletion ALT length must be <= REF length")]
    SvDeletionLen,
    #[error("SV insertion ALT length must be >= REF length")]
    SvInsertionLen,
    #[error("TR contraction ALT length must be < REF length")]
    TrContractionLen,
    #[error("TR expansion ALT length must be >= REF length")]
    TrExpansionLen
}

/// A variant definition structure.
/// It currently assumes that chromosome is fixed and that the variant is a SNP.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Variant {
    /// The vcf index from the input datasets, use 0 if unneeded
    vcf_index: usize,
    /// The type of variant represented by this entry
    variant_type: VariantType,
    /// The coordinate of the event in the VCF file, 0-based
    position: u64,
    /// the first allele value
    allele0: Vec<u8>,
    /// the second allele value
    allele1: Vec<u8>,

    // auxiliary booleans
    /// if true, flags this as a variant to ignore for _some_ reason
    is_ignored: bool
}

impl Variant {
    /// Creates a new single-nucleotide variant (SNV).
    /// For SNV variants, all alleles must be exactly 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match a single-nucleotide variant
    pub fn new_snv(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // SNV alleles must be length 1
        if allele0.len() != 1 {
            return Err(VariantError::AlleleLen1 { index: 0 });
        }
        if allele1.len() != 1 {
            return Err(VariantError::AlleleLen1 { index: 1 });
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Snv,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    /// Creates a new deletion variant.
    /// Deletions must have a REF allele long than 1 bp, and all ALT alleles must be exactly 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match a deletion variant
    pub fn new_deletion(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // reference length must be greater than 1 to be a deletion
        if allele0.len() <= 1 {
            return Err(VariantError::RefLenGT1);
        }

        // this one must always be length 1
        if allele1.len() != 1 {
            return Err(VariantError::AlleleLen1 { index: 1 });
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Deletion,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    /// Creates a new insertion variant.
    /// Insertions must have a REF allele exactly 1 bp long, and all ALT alleles must be longer than 1 bp.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match an insertion variant
    pub fn new_insertion(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // if reference allele is present, it must be length 1 for this type
        if allele0.len() != 1 {
            return Err(VariantError::AlleleLen1 { index: 0 });
        }

        if allele1.len() <= 1 {
            return Err(VariantError::AltLenGT1);
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Insertion,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    /// Creates a new indel variant.
    /// All indels alleles must be more than 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match an indel variant
    pub fn new_indel(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // REF and ALT lengths must be greater than 1 to be an indel
        if allele0.len() <= 1 {
            return Err(VariantError::RefLenGT1);
        }

        if allele1.len() <= 1 {
            return Err(VariantError::AltLenGT1);
        }
        
        // there's no real reason to believe in any shared sequence between alleles
        // we've seen it not work above, not worth even trying to codify warning here IMO
        // assert!(???)
        
        Ok(Variant {
            vcf_index,
            variant_type: VariantType::Indel,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    /// Creates a new SV deletion variant.
    /// SV deletions must have a REF allele long than 1 bp, and all ALT alleles must be exactly 1 bp long.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match a deletion variant
    /// * if the reference allele is passed in and it does not have the same length as `ref_len`
    pub fn new_sv_deletion(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // reference length must be greater than 1 to be a deletion
        if allele0.len() <= 1 {
            return Err(VariantError::RefLenGT1);
        }

        // restriction lifted such that now allele0 must be >= allele1
        if allele1.len() > allele0.len() {
            return Err(VariantError::SvDeletionLen);
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::SvDeletion,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    /// Creates a new SV insertion variant.
    /// SV insertions must have a REF allele exactly 1 bp long, and all ALT alleles must be longer than 1 bp.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match an insertion variant
    pub fn new_sv_insertion(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // restriction lifted such that now allele1 must be >= allele0
        if allele1.len() < allele0.len() {
            return Err(VariantError::SvInsertionLen);
        }
        
        // allele0 is always <= allele1, so just make sure it's not empty
        if allele0.is_empty() {
            return Err(VariantError::EmptyAllele { index: 0 });
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::SvInsertion,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    /// Creates a new tandem repeat contraction variant, functionally these act very similar to indel types.
    /// All tandem repeat alleles must be at least 1 bp long by VCF definition.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match a tandem repeat variant
    pub fn new_tr_contraction(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // all alleles must be >= 1 for tandem repeats, most are longer though
        if allele0.is_empty() {
            return Err(VariantError::EmptyAllele { index: 0 });
        }
        if allele1.is_empty() {
            return Err(VariantError::EmptyAllele { index: 1 });
        }

        // allele0 must be > allele1
        if allele1.len() >= allele0.len() {
            return Err(VariantError::TrContractionLen);
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::TrContraction,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    /// Creates a new tandem repeat expansion variant, functionally these act very similar to indel types.
    /// All tandem repeat alleles must be at least 1 bp long by VCF definition.
    /// # Arguments
    /// * `vcf_index` - the index of the source VCF file
    /// * `position` - the coordinate of the variant in a contig
    /// * `allele0` - the first allele (usually REF)
    /// * `allele1` - the second allele (usually ALT[0])
    /// # Errors
    /// * if the provided sequences do not match a tandem repeat variant
    pub fn new_tr_expansion(vcf_index: usize, position: u64, allele0: Vec<u8>, allele1: Vec<u8>) -> Result<Variant, VariantError> {
        // all alleles must be >= 1 for tandem repeats, most are longer though
        if allele0.is_empty() {
            return Err(VariantError::EmptyAllele { index: 0 });
        }
        if allele1.is_empty() {
            return Err(VariantError::EmptyAllele { index: 1 });
        }

        // allele0 must be <= allele1
        if allele1.len() < allele0.len() {
            return Err(VariantError::TrExpansionLen);
        }

        Ok(Variant {
            vcf_index,
            variant_type: VariantType::TrExpansion,
            position,
            allele0,
            allele1,
            is_ignored: false
        })
    }

    pub fn set_ignored(&mut self) {
        self.is_ignored = true;
    }

    /// This will determine the best matching allele (0 or 1) or return 2 if neither match.
    /// Primary purpose of this is to convert all variant observations into a 0/1 scheme.
    /// This method requires an exact match of the allele.
    /// # Arguments
    /// * `allele` - the allele that needs to get converted to a 0 or 1 (or 2 if neither match)
    pub fn match_allele(&self, allele: &[u8]) -> u8 {
        if allele == &self.allele0[..] {
            0
        } else if allele == &self.allele1[..] {
            1
        } else {
            2
        }
    }

    /// This will return the index allele for a given haplotype index.
    /// Input must always be 0 or 1, but it might get converted to something else at multi-allelic sites.
    /// # Arguments
    /// * `index` - must be 0, 1, or 2 (unknown)
    /// # Panics
    /// * if anything other than 0, 1, or 2 is provided
    pub fn convert_index(&self, index: u8) -> u8 {
        // This is a legacy function from when we had multi-allelics; if we ever allow that again, we will need to patch this
        if index == 0 {
            0
        } else if index == 1 {
            1
        } else if index == 2 {
            // we just need some indicator that it's undetermined, this will work for now
            u8::MAX  
        } else {
            panic!("index must be 0, 1, or 2");
        }
    }

    /// Calculates the REF/ALT edit distance for this variant
    pub fn alt_ed(&self) -> anyhow::Result<usize> {
        wfa_ed(&self.allele0, &self.allele1)
    }

    // getters
    pub fn vcf_index(&self) -> usize {
        self.vcf_index
    }

    pub fn variant_type(&self) -> VariantType {
        self.variant_type
    }

    pub fn position(&self) -> u64 {
        self.position
    }

    pub fn ref_len(&self) -> usize {
        self.allele0.len()
    }

    pub fn allele0(&self) -> &[u8] {
        &self.allele0
    }

    pub fn allele1(&self) -> &[u8] {
        &self.allele1
    }

    pub fn is_ignored(&self) -> bool {
        self.is_ignored
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_snv() {
        let variant = Variant::new_snv(
            0, 1,
            b"A".to_vec(), b"C".to_vec(),
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::Snv);
        assert_eq!(variant.position(), 1);
        assert_eq!(variant.ref_len(), 1);
        assert_eq!(variant.match_allele(b"A"), 0);
        assert_eq!(variant.match_allele(b"C"), 1);
        assert_eq!(variant.match_allele(b"G"), 2);
        assert_eq!(variant.match_allele(b"T"), 2);
    }

    #[test]
    fn test_basic_deletion() {
        // this is the deletion we mostly expect
        let variant = Variant::new_deletion(
            0, 10,
            b"AGT".to_vec(), b"A".to_vec()
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::Deletion);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.ref_len(), 3);
        assert_eq!(variant.match_allele(b"AGT"), 0);
        assert_eq!(variant.match_allele(b"A"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_basic_insertion() {
        let variant = Variant::new_insertion(
            0, 20,
            b"A".to_vec(), b"AGT".to_vec()
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::Insertion);
        assert_eq!(variant.position(), 20);
        assert_eq!(variant.ref_len(), 1);
        assert_eq!(variant.match_allele(b"A"), 0);
        assert_eq!(variant.match_allele(b"AGT"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_basic_indel() {
        // models AG -> A / AGT
        let variant = Variant::new_indel(
            0, 20,
            b"AG".to_vec(), b"AGT".to_vec()
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::Indel);
        assert_eq!(variant.position(), 20);
        assert_eq!(variant.ref_len(), 2);
        assert_eq!(variant.match_allele(b"A"), 2);
        assert_eq!(variant.match_allele(b"AGT"), 1);
        assert_eq!(variant.match_allele(b"AG"), 0);
    }

    #[test]
    fn test_sv_insertion() {
        let variant = Variant::new_sv_insertion(
            0, 20,
            b"A".to_vec(), b"AGT".to_vec()
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::SvInsertion);
        assert_eq!(variant.position(), 20);
        assert_eq!(variant.ref_len(), 1);

        // TODO: replace this with the matching we will do with SVs
        assert_eq!(variant.match_allele(b"A"), 0);
        assert_eq!(variant.match_allele(b"AGT"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_sv_deletion() {
        let variant = Variant::new_sv_deletion(
            0, 10,
            b"AGT".to_vec(), b"A".to_vec()
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::SvDeletion);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.ref_len(), 3);

        // TODO: replace this with the matching we will do with SVs
        assert_eq!(variant.match_allele(b"AGT"), 0);
        assert_eq!(variant.match_allele(b"A"), 1);
        assert_eq!(variant.match_allele(b"AG"), 2);
    }

    #[test]
    fn test_tr_expansion() {
        let variant = Variant::new_tr_expansion(
            0, 10,
            b"AAAC".to_vec(),
            b"AAACAAAC".to_vec()
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::TrExpansion);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.ref_len(), 4);

        assert_eq!(variant.match_allele(b"AAAC"), 0);
        assert_eq!(variant.match_allele(b"AAACAAAC"), 1);
        assert_eq!(variant.match_allele(b"AAACAA"), 2);
    }

    #[test]
    fn test_tr_contraction() {
        let variant = Variant::new_tr_contraction(
            0, 10,
            b"AAACAAAC".to_vec(),
            b"AAAC".to_vec(),
        ).unwrap();
        assert_eq!(variant.variant_type(), VariantType::TrContraction);
        assert_eq!(variant.position(), 10);
        assert_eq!(variant.ref_len(), 8);

        assert_eq!(variant.match_allele(b"AAACAAAC"), 0);
        assert_eq!(variant.match_allele(b"AAAC"), 1);
        assert_eq!(variant.match_allele(b"AAACAA"), 2);
    }
}
