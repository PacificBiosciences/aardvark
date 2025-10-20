# Recommended settings
This document contains recommended settings for running Aardvark on different types of variants.

Table of contents:

* [Compare mode](#compare-mode)
  * [Comparing small variants](#comparing-small-variants)
  * [Comparing structural variants](#comparing-structural-variants)
  * [Comparing tandem repeat variants](#comparing-tandem-repeat-variants)

# Compare mode
## Comparing small variants
The default settings are configured for small-variant only runs.
This includes SNVs and indels < 50bp.

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
    --enable-record-basepair \
    ...
```
