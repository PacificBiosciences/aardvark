
use anyhow::ensure;

use crate::data_types::coordinates::Coordinates;
use crate::data_types::multi_region::MultiRegion;
use crate::data_types::phase_enums::PhasedZygosity;
use crate::data_types::variants::Variant;

/// Structure containing all the necessary information to solve a problem.
/// In theory, almost no file I/O should be necessary after loading a problem.
/// Intended to be a self-contained unit-of-work that can be parallelized out.
#[derive(Debug)]
pub struct CompareRegion {
    /// Unique identifier for the comparison region
    region_id: u64,
    /// The full region we are comparing
    coordinates: Coordinates,
    /// Truth variants
    truth_variants: Vec<Variant>,
    /// Zygosity of the truth variants
    truth_zygosity: Vec<PhasedZygosity>,
    /// Query variants
    query_variants: Vec<Variant>,
    /// Zygosity of the query variants
    query_zygosity: Vec<PhasedZygosity>,
}

impl CompareRegion {
    /// General constructor with checks
    /// # Arguments
    /// * `region_id` - unique ID for the region
    /// * `coordinates` - coordinates for the full region getting compared
    /// * `truth_variants` - variants in the truth set
    /// * `truth_zygosity` - the zygosity of the variants in the truth set
    /// # Errors
    /// * if the number of truth variants and calls are different
    pub fn new(
        region_id: u64, coordinates: Coordinates,
        truth_variants: Vec<Variant>, truth_zygosity: Vec<PhasedZygosity>,
        query_variants: Vec<Variant>, query_zygosity: Vec<PhasedZygosity>
    ) -> anyhow::Result<Self> {
        // verify equal lengths of inputs
        ensure!(truth_variants.len() == truth_zygosity.len(), "Number of truth variants and zygosities must be equal");
        ensure!(query_variants.len() == query_zygosity.len(), "Number of query variants and zygosities must be equal");

        // verify there are no multi-ALT sites provided; i.e., they have been pre-split for us
        ensure!(
            truth_variants.iter().all(|v| v.convert_index(0) == 0),
            "Unsupported multi-alt variant detected in truth_variants"
        );ensure!(
            query_variants.iter().all(|v| v.convert_index(0) == 0),
            "Unsupported multi-alt variant detected in query_variants"
        );
        Ok(Self {
            region_id, coordinates,
            truth_variants, truth_zygosity,
            query_variants, query_zygosity
        })
    }

    /// This will return the coordinates containing all of the truth and query variants.
    /// It should always be a subset of `coordinates` if the inputs are valid.
    pub fn var_coordinates(&self) -> Coordinates {
        let chrom = self.coordinates.chrom().to_string();
        let min_truth = self.truth_variants.first().map(|v| v.position()).unwrap_or(u64::MAX);
        let min_query = self.query_variants.first().map(|v| v.position()).unwrap_or(u64::MAX);
        let start = min_truth.min(min_query);

        let max_truth = self.truth_variants.last().map(|v| v.position() + v.ref_len() as u64).unwrap_or(u64::MIN);
        let max_query = self.query_variants.last().map(|v| v.position() + v.ref_len() as u64).unwrap_or(u64::MIN);
        let end = max_truth.max(max_query);

        Coordinates::new(chrom, start, end)
    }

    // various getters
    pub fn region_id(&self) -> u64 {
        self.region_id
    }

    pub fn coordinates(&self) -> &Coordinates {
        &self.coordinates
    }

    pub fn truth_variants(&self) -> &[Variant] {
        &self.truth_variants
    }

    pub fn truth_zygosity(&self) -> &[PhasedZygosity] {
        &self.truth_zygosity
    }

    pub fn query_variants(&self) -> &[Variant] {
        &self.query_variants
    }

    pub fn query_zygosity(&self) -> &[PhasedZygosity] {
        &self.query_zygosity
    }
}

impl TryFrom<MultiRegion> for CompareRegion {
    type Error = anyhow::Error;

    fn try_from(value: MultiRegion) -> Result<Self, Self::Error> {
        ensure!(value.variants().len() == 2, "Cannot convert MultiRegion with variants != 2");
        ensure!(value.zygosity().len() == 2, "Cannot convert MultiRegion with zygosity != 2");

        Ok(Self {
            region_id: value.region_id(),
            coordinates: value.coordinates().clone(),
            truth_variants: value.variants()[0].clone(),
            truth_zygosity: value.zygosity()[0].clone(),
            query_variants: value.variants()[1].clone(),
            query_zygosity: value.zygosity()[1].clone(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compare_region() {
        // basic single variant SNV region
        let _good_region = CompareRegion::new(0,
            Coordinates::new("chr1".to_string(), 0, 20),
            vec![Variant::new_snv(0, 10, b"A".to_vec(), b"C".to_vec()).unwrap()],
            vec![PhasedZygosity::UnphasedHeterozygous],
            vec![],
            vec![]
        ).unwrap();

        // region with un-equal inputs
        let bad_region_result = CompareRegion::new(1,
            Coordinates::new("chr1".to_string(), 0, 20),
            vec![],
            vec![PhasedZygosity::UnphasedHeterozygous],
            vec![],
            vec![]
        );
        assert!(bad_region_result.is_err());

        // region with un-equal calls
        let bad_region_result = CompareRegion::new(1,
            Coordinates::new("chr1".to_string(), 0, 20),
            vec![],
            vec![],
            vec![],
            vec![PhasedZygosity::UnphasedHeterozygous]
        );
        assert!(bad_region_result.is_err());
    }
}