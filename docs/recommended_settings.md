# Recommended settings
This document contains recommended settings for running Aardvark on different types of variants.

Table of contents:

* [Compare mode](#compare-mode)
  * [Comparing small variants](#comparing-small-variants)
  * [Comparing small and medium structural variants](#comparing-small-and-medium-structural-variants)
  * [Comparing structural variants](#comparing-structural-variants)
  * [Comparing tandem repeat variants](#comparing-tandem-repeat-variants)

# Compare mode
The [compare presets](./compare.md#presets) cover the common cases, which are outlined here.

## Comparing small variants
Use `--preset small` for SNVs and indels.
Alleles at or below 50 bp are always scored.
Alleles from 51 bp through 200 bp are included only when they improve the comparison.
Alleles with event length above 200 bp are ignored.

```
aardvark compare \
    --preset small \
    ...
```

## Comparing small and medium structural variants
Use `--preset medium` when the call sets include small variants and structural variants up to 10 kbp.
The variant gap is 1 kbp, so a longer allele stays in the same sub-region as nearby smaller alleles.
Alleles from 10 kbp through 20 kbp are included only when they improve the comparison, and longer event lengths are ignored.

```
aardvark compare \
    --preset medium \
    ...
```

## Comparing structural variants
Use `--preset large` when structural variants should be compared without a length band from the preset.
The variant gap is 1 kbp.
Any length band passed on the command line is kept.
With no length band, a REF or ALT sequence longer than 1 Mbp is still ignored.

```
aardvark compare \
    --preset large \
    ...
```

## Comparing tandem repeat variants
When adding tandem repeat variants, the [medium preset](#comparing-small-and-medium-structural-variants) is the recommended approach.
Depending on how the variants are represented, the `RECORD_BP` scoring scheme may provide a complementary score for the tandem repeat variants.
See the [record basepair description](./methods.md#record-basepair) for more details:

```
aardvark compare \
    --preset medium \
    --enable-record-basepair-metrics \
    ...
```
