<h1 align="center"><img width="300px" src="docs/img/logo_aardvark.svg"/></h1>

<h1 align="center">Aardvark</h1>

<p align="center">A tool for sniffing out the differences in vari-Ants</p>

***

Aardvark was built to quickly compare different variant call sets.
Key features of Aardvark include:

* Ability to [benchmark](./docs/compare.md) variants from a "query" set against those from a "truth" set
  * Constructs haplotype sequences, allowing for basepair-level comparisons that are variant-type agnostic and that enable implicit partial credit for inexact matching variant calls
  * Genotype assessment, allowing for a strict exact-match scoring (similar to traditional metrics)
  * Provides output summary statistics and VCF files with annotations
* Ability to [merge](./docs/merge.md) variants from multiple input sets
  * Variants are compared in a haplotype context (sequence level), harmonizing different variant representations
  * Multiple merge strategies for different scenarios from strict exact-matching to looser majority-vote schemes
* Quality of life additions: Efficient methods with a fast run-time distributed as a single binary

Authors: [Matt Holt](https://github.com/holtjma), [Zev Kronenberg](https://github.com/zeeev)

## Availability
* [Latest release with binary](https://github.com/PacificBiosciences/Aardvark/releases/latest)

## Documentation
* [Installation instructions](./docs/install.md)
* [User guide with quickstart](./docs/user_guide.md)
* [Output files](./docs/user_guide.md#output-files)
* [Recommended settings](./docs/recommended_settings.md)
* [Methods](./docs/methods.md)
* [Performance](./docs/performance.md)

## Citation
If you use Aardvark, please cite our bioRxiv pre-print:

[Holt, James Matthew, et al. "Aardvark: Sifting through differences in a mound of variants." bioRxiv (2025): 2025-10.](https://doi.org/10.1101/2025.10.03.680257)

## Support information
If you notice any missing features, bugs, or need assistance with analyzing the output of Aardvark, 
please [open a GitHub issue](https://github.com/PacificBiosciences/aardvark/issues).

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
