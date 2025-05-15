<h1 align="center"><img width="300px" src="images/logo_aardvark.svg"/></h1>

<h1 align="center">Aardvark</h1>

<p align="center">A tool sniffing out the differences in vari-Ants</p>

***

Aardvark is designed to quickly compare variant calls from different call sets.
Key features of Aardvark include:

* Ability to [benchmark](./docs/compare.md) variants from a "query" set against those from a "truth" set
  * Constructs haplotype sequences, allowing for basepair-level comparisons that are variant-type agnostic and that enable implicit partial credit for inexact matching variant calls
  * Assess variants by genotype, allowing for more traditional exact-matching analyses and scoring
  * Provides output summary statistics and VCF files with annotations
* Ability to [merge](./docs/merge.md) variants from multiple input sets
  * Variants are compared at the haplotype sequence (basepair) level, allowing for merging despite variant representation
  * Multiple merge strategies for different scenarios from strict exact-matching to looser majority-vote schemes
* Efficient methods, a typical aardvark comparison will finish in <5 minutes of CPU time

Authors: [Matt Holt](https://github.com/holtjma), [Zev Kronenberg](https://github.com/zeeev)

## Availability
* [Latest release with binary](https://github.com/PacificBiosciences/Aardvark/releases/latest)

## Documentation
* [Installation instructions](./docs/install.md)
* [User guide with quickstart](./docs/user_guide.md)
* [Output files](./docs/user_guide.md#output-files)
* [Methods](./docs/methods.md)
* [Performance](./docs/performance.md)

## Citation
Aardvark does not currently have a publication associated with it.

## Need help?
If you notice any missing features, bugs, or need assistance with analyzing the output of Aardvark, 
please don't hesitate to open a GitHub issue.

## Support information
Aardvark is a pre-release software intended for research use only and not for use in diagnostic procedures. 
While efforts have been made to ensure that Aardvark lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As Aardvark is not covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any Aardvark release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
