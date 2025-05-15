

/// Simple enum to restrict thing to hap1 or hap2
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Haplotype {
    Hap1,
    Hap2
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub enum Allele {
    /// Indicates an unknown allele, usually '.' in a file
    Unknown,
    /// Typically refers to reference allele, usually 0 in a file; could also mean the first alternate in multi-allelic context
    Reference,
    /// Typically refers to alternate allele, usually 1 in a file; could also mean the second alternate in multi-allelic context
    Alternate
}

impl Allele {
    /// Converts into a basic count representation. Useful for counting ALT alleles in a batch.
    pub fn to_allele_count(&self) -> u8 {
        match self {
            Allele::Unknown |
            Allele::Reference => 0,
            Allele::Alternate => 1,
        }
    }
}

/// Captures the zygosity types we can observe in a diploid organism
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhasedZygosity {
    /// ./., but could just be anything unhandled so far
    Unknown,
    /// 0/0
    HomozygousReference,
    /// 0/1, 1/2, etc.
    UnphasedHeterozygous,
    /// 0|1, 1|2, etc.
    PhasedHet01,
    /// 1|0, 2|1, etc.
    PhasedHet10,
    /// 1/1, 2/2, etc.
    HomozygousAlternate

    // TODO: do we need to handle haploid as well? for now, we can throw it into unknown
}

impl PhasedZygosity {
    /// Returns true if this is a phased heterozygous variant.
    pub fn is_phased(&self) -> bool {
        match self {
            PhasedZygosity::PhasedHet01 |
            PhasedZygosity::PhasedHet10 => true,

            PhasedZygosity::Unknown |
            PhasedZygosity::HomozygousReference |
            PhasedZygosity::UnphasedHeterozygous |
            PhasedZygosity::HomozygousAlternate => false,
        }
    }

    /// Returns true if this is a heterozygous variant
    pub fn is_heterozygous(&self) -> bool {
        match self {
            PhasedZygosity::UnphasedHeterozygous |
            PhasedZygosity::PhasedHet01 |
            PhasedZygosity::PhasedHet10 => true,

            PhasedZygosity::Unknown |
            PhasedZygosity::HomozygousReference |
            PhasedZygosity::HomozygousAlternate => false,
        }
    }

    /// Returns true if this is a homozygous variant
    pub fn is_homozygous(&self) -> bool {
        match self {
            PhasedZygosity::HomozygousReference |
            PhasedZygosity::HomozygousAlternate => true,

            PhasedZygosity::UnphasedHeterozygous |
            PhasedZygosity::PhasedHet01 |
            PhasedZygosity::PhasedHet10 |
            PhasedZygosity::Unknown => false
        }
    }

    /// Splits the zygosity into two Alleles
    pub fn decompose_alleles(&self) -> (Allele, Allele) {
        match self {
            PhasedZygosity::Unknown => (Allele::Unknown, Allele::Unknown),
            PhasedZygosity::HomozygousReference => (Allele::Reference, Allele::Reference),
            PhasedZygosity::UnphasedHeterozygous |
            PhasedZygosity::PhasedHet01 => (Allele::Reference, Allele::Alternate),
            PhasedZygosity::PhasedHet10 => (Allele::Alternate, Allele::Reference),
            PhasedZygosity::HomozygousAlternate => (Allele::Alternate, Allele::Alternate),
        }
    }
}
