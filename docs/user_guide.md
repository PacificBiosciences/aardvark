# User Guide
Table of contents:

* [Installation instruction](./install.md)
* [Main workflows](#main-workflows)
* [Input preparation](#input-preparation)

# Main workflows
There are currently two main workflows that are supported by Aardvark:
1. [Compare](./compare.md) - Benchmarks a "query" variant call set against a "truth" variant call set, reporting summary level statistics such as precision, recall, and F1 score. It also generates VCF files with the classifications of individual variants.
2. [Merge](./merge.md) - Given two or more input call sets, this will compare them against each other and "merge" or "collapse" regions with identical basepair-level haplotype sequences. Multiple merge strategies are available to enable some forms of conflict to get merged also.

# Input preparation
Aardvark does not apply any filtering to input VCF records for common options like the `FILTER` or `QUAL` fields of the VCF.
All variants from the input VCFs that fall within the provided `--regions` are used as-is.
If you want to restrict the analysis to specific variants, pre-filter your VCF files before passing them to Aardvark; for example:

```bash
# keep variants with no filter
bcftools view -f "PASS,." input.vcf.gz -Oz -o filtered.vcf.gz
```

See the [bcftools documentation](https://samtools.github.io/bcftools/bcftools.html#view) for a comprehensive set of common filtration options.
