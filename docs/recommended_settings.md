# Recommended settings
This document contains recommended settings for running Aardvark on different types of variants.

Table of contents:

* [Compare mode](#compare-mode)
  * [Comparing small variants](#comparing-small-variants)
  * [Comparing a mixed callset to small variants](#comparing-a-mixed-callset-to-small-variants)
  * [Comparing structural variants](#comparing-structural-variants)
  * [Comparing tandem repeat variants](#comparing-tandem-repeat-variants)

# Compare mode
## Comparing small variants
The default settings are configured for small-variant only runs.
This includes SNVs and indels < 50bp.

## Comparing a mixed callset to small variants
A query or truth VCF may contain both small variants and longer alleles.
Dropping every allele longer than 50 bp before comparison can hide a case where one 51 bp allele is the same sequence change as two smaller alleles, such as a 25 bp deletion plus a 26 bp deletion.

Set a soft max and a hard max so alleles between them are kept only when they lower the truth-query edit distance:

```
aardvark compare \
    --soft-max-variant-length 50 \
    --hard-max-variant-length 100 \
    ...
```

Alleles at or below 50 bp are always scored.
Alleles from 51 bp through 100 bp are included only when they improve the sequence match.
Alleles longer than 100 bp are ignored.
The same band is applied to both VCFs.
Increase `--min-variant-gap` when the long allele and the shorter alleles would otherwise fall into different sub-regions.

## Comparing structural variants
When adding structural variants, we recommend increasing the window size parameter of Aardvark.
This tends to increase the final scores at the cost of additional compute time:

```
aardvark compare \
    --min-variant-gap 1000 \
    ...
```

## Comparing tandem repeat variants
When adding tandem repeat variants, we recommend increasing the window size parameter of Aardvark.
This tends to increase the final scores at the cost of additional compute time.
Depending on how the variants are represented, the `RECORD_BP` scoring scheme may provide a complementary score for the tandem repeat variants.
See the [record basepair description](./methods.md#record-basepair) for more details:

```
aardvark compare \
    --min-variant-gap 1000 \
    --enable-record-basepair-metrics \
    ...
```
