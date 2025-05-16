# User Guide
Table of contents:

* [Installation instruction](./install.md)
* [Main workflows](#main-workflows)

# Main workflows
There are currently two main workflows that are supported by Aardvark:
1. [Compare](./compare.md) - Benchmarks a "query" variant call set against a "truth" variant call set, reporting summary level statistics such as precision, recall, and F1 score. It also generates VCF files with the classifications of individual variants.
2. [Merge](./merge.md) - Given two or more input call sets, this will compare them against each other and "merge" or "collapse" regions with identical basepair-level haplotype sequences. Multiple merge strategies are available to enable some forms of conflict to get merged also.
