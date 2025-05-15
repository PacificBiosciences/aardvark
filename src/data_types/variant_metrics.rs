
use anyhow::bail;

/// Each variant has a source
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum VariantSource {
    Truth,
    Query
}

/// Each variant has a classification
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd, strum_macros::AsRefStr)]
pub enum Classification {
    #[strum(serialize = "UNK")]
    Unknown=0,
    #[strum(serialize = "TP")]
    TruePositive,
    #[strum(serialize = "FN")]
    FalseNegative,
    #[strum(serialize = "FP")]
    FalsePositive
}

/// Collects variant-level assignments based on the GT comparison
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VariantMetrics {
    /// Source of the variant
    source: VariantSource,
    /// Classification, which is based on the entire GT
    classification: Classification,
    /// Expected allele count, range: [0, 2]
    expected_allele_count: u8,
    /// Observed allele count, range: [0, 2]
    observed_allele_count: u8
}

impl VariantMetrics {
    /// Constructor
    /// # Arguments
    /// * `source` - the variant source VCF
    /// * `expected_allele_count` - the number of expected alleles
    /// * `observed_allele_count` - the number of observed alleles; variant classification made by comparing to `expected_allele_count`
    pub fn new(source: VariantSource, expected_allele_count: u8, observed_allele_count: u8) -> anyhow::Result<Self> {
        if expected_allele_count > 2 {
            bail!("Expected allele count must be in range: [0, 2]");
        }
        if observed_allele_count > 2 {
            bail!("Observed allele count must be in range: [0, 2]");
        }

        let classification = match expected_allele_count.cmp(&observed_allele_count) {
            // expected < observed, we have extra that are FP
            std::cmp::Ordering::Less => Classification::FalsePositive,
            // expected == observed, great!
            std::cmp::Ordering::Equal => {
                if expected_allele_count == 0 {
                    bail!("Variant metrics does not support expected and observed allele counts of 0");
                }
                Classification::TruePositive
            },
            // expected > observed, we are missing some that are FN
            std::cmp::Ordering::Greater => Classification::FalseNegative
        };
        
        Ok(Self {
            source,
            classification,
            expected_allele_count,
            observed_allele_count
        })
    }

    /// Given a previous set of VariantMetrics, this will convert the orientation of the comparison.
    /// Source is toggled, classifications are adjusted (FP <=> FN), and expected/observed counts are swapped.
    /// # Arguments
    /// * `original` - the source VariantMetric we are toggling
    pub fn toggle_source(original: &Self) -> Self {
        let source = match original.source() {
            VariantSource::Truth => VariantSource::Query,
            VariantSource::Query => VariantSource::Truth
        };
        let classification = match original.classification() {
            // these stay the same
            Classification::Unknown => Classification::Unknown,
            Classification::TruePositive => Classification::TruePositive,
            // these two switch
            Classification::FalseNegative => Classification::FalsePositive,
            Classification::FalsePositive => Classification::FalseNegative
        };

        // swap these
        let expected_allele_count = original.observed_allele_count();
        let observed_allele_count = original.expected_allele_count();

        Self {
            source,
            classification,
            expected_allele_count,
            observed_allele_count
        }
    }

    // getters
    pub fn source(&self) -> VariantSource {
        self.source
    }

    pub fn classification(&self) -> Classification {
        self.classification
    }

    pub fn expected_allele_count(&self) -> u8 {
        self.expected_allele_count
    }

    pub fn observed_allele_count(&self) -> u8 {
        self.observed_allele_count
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_toggle_source() {
        let original = VariantMetrics::new(VariantSource::Truth, 0, 1).unwrap();
        assert_eq!(original.classification(), Classification::FalsePositive);

        // check the toggled result
        let toggled = VariantMetrics::toggle_source(&original);
        assert_eq!(toggled, VariantMetrics::new(VariantSource::Query, 1, 0).unwrap());
        assert_eq!(toggled.classification(), Classification::FalseNegative);

        // toggling again should get the original back
        let double_toggled = VariantMetrics::toggle_source(&toggled);
        assert_eq!(double_toggled, original);
    }
}