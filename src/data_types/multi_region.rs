
use anyhow::ensure;

use crate::data_types::coordinates::Coordinates;
use crate::data_types::variants::Variant;
use crate::data_types::phase_enums::PhasedZygosity;

/// Data structure with a generic number of input VCF files that variants can get pulled from
#[derive(Debug)]
pub struct MultiRegion {
    /// Unique identifier for the comparison region
    region_id: u64,
    /// The full region we are comparing
    coordinates: Coordinates,
    /// List of variants, one per input
    variants: Vec<Vec<Variant>>,
    /// List of zygosity, one per input and same length as corresponding entry in `variants`
    zygosity: Vec<Vec<PhasedZygosity>>,
}

impl MultiRegion {
    /// General constructor with checks
    /// # Arguments
    /// * `region_id` - unique ID for the region
    /// * `coordinates` - coordinates for the full region getting compared
    /// * `variants` - variants in each set
    /// * `zygosity` - the zygosity of the variants in each set
    /// # Errors
    /// * if the number of entries in `variants` and `zygosity` is different
    /// * if within an entry, the length of variants and zygosity are different
    /// * if any variants are multi-allelic
    pub fn new(
        region_id: u64, coordinates: Coordinates,
        variants: Vec<Vec<Variant>>, zygosity: Vec<Vec<PhasedZygosity>>
    ) -> anyhow::Result<Self> {
        ensure!(variants.len() == zygosity.len(), "Number of entries in variants and zygosity must be equal");

        // verify equal lengths of inputs
        for (i, (v, z)) in variants.iter().zip(zygosity.iter()).enumerate() {
            ensure!(v.len() == z.len(), "Number of truth variants and zygosities are not equal for input {i}");

            // verify no multi-ALT inputs
            ensure!(
                v.iter().all(|v| v.convert_index(0) == 0),
                "Unsupported multi-alt variant detected in input {i}"
            );
        }

        Ok(Self {
            region_id, coordinates,
            variants, zygosity
        })
    }

    // various getters
    pub fn region_id(&self) -> u64 {
        self.region_id
    }

    pub fn coordinates(&self) -> &Coordinates {
        &self.coordinates
    }

    pub fn variants(&self) -> &[Vec<Variant>] {
        &self.variants
    }

    pub fn zygosity(&self) -> &[Vec<PhasedZygosity>] {
        &self.zygosity
    }
}