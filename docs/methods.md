# Methods
Table of contents:

* [Overview](#overview)
  * [Sub-region identification](#sub-region-identification)
  * [Compare command](#compare-command)
  * [Merge command](#merge-command)
* [Comparison types](#comparison-types)
  * [Basepair](#basepair)
  * [Haplotype](#haplotype)
  * [Genotype](#genotype)

# Overview
Both `compare` and `merge` commands follow a three step process:

1. Sub-region identification - Load the input files and identify sub-regions of the input that can be analyzed independently.
2. Comparison - For each sub-region, perform the desired comparison. This step is performed in parallel since each sub-region is independent. We include descriptions for each below.
3. Output - Create any summary or joint outputs and write them to file. We refer to our users guide for details on output files.

## Sub-region identification
Aardvark starts by identifying sub-regions of the full input that can be solved independently.
First, all valid BED regions are loaded from the user-provided input (`--confidence-region` or `--regions`).
Aardvark then loads all variants that are **fully-contained** within each region from all provided VCF files.
Lastly, Aardvark splits the BED regions further by looking for gaps between consecutive variants across all files.
Any gaps that larger than a defined window (50 bp) are split into sub-regions.
Effectively, this clusters variants that are within 50 bp of each other into sub-regions, while accounting for gaps in the user-provided BED regions.

## Compare command
In `compare` mode, Aardvark performs a comparison between "truth" and "query" sets.
Each sub-region is compared independently, in parallel, using the following steps:

1. Identify the optimal diplotype sequences - In this step, the algorithm searches for the optimal set of *phased* variant zygosities that minimizes the edit distance between the truth haplotype sequences and the query haplotype sequences. If pre-phased inputs are provided in the truth set, they are retained by Aardvark. In contrast, query phasing is ignored, allowing for phasing errors in the query set. This step outputs the optimal variant phasing and haplotype sequence pairs such that the edit distance between truth and query is minimized.
2. Compare the optimized sequences - The sequences are directly compared to calculate the total number of basepairs that are true positives (TP), false negatives (FN), and false positives (FP). This is calculated by measuring the edit distances between the reference genome, the truth optimal haplotype sequence, and the query optimal haplotype sequence. These values go into a system of equations to solve for TP, FN, and FP. The final results goes into the `BASEPAIR` comparison type.
3. Compare the optimized alleles - Aardvark will also compare the individual alleles on a truth and query haplotype to determine if each is a TP, FN, or FP. Given a set of present alternate alleles on the truth and query haplotypes, it will search for a solution with the fewest number of errors (FN or FP) such that the generated haplotype sequences are identical. Conceptually, this process is very similar to that of [rtg vcfeval](https://github.com/RealTimeGenomics/rtg-tools). The results from each individual haplotype go into the `HAP` and `WEIGHTED_HAP` comparison types. The `HAP` evaluations are combined into a genotype-level score in the `GT` comparison type, where an error in either haplotype counts as an error in the `GT` score.

## Merge command
In `merge` mode, Aardvark performs comparisons between multiple input sets.
Each sub-region is compared independently, in parallel, using the following steps:

1. For each pair of inputs, compare the `BASEPAIR` scores - This process uses the same approach described in the [Compare command](#compare-command). However, Aardvark is only searching for identical regions (e.g., `BASEPAIR` F1 score = 1.0), indicating that the variant sets represent the same collection of variants.
2. The pair-wise `BASEPAIR` scores are used to categorize/label each sub-region. All non-`different` sub-regions are written to the output merged file. We briefly describe each label type, which are presented the priority order used for labeling:
  * `identical` - Indicates that all inputs have identical basepair-level haplotype sequences, even if the representations of those variants and zygosities are not identical. The variants written to output are always from the first input VCF.
  * `no_conflict` - Indicates that all inputs either 1) have identical basepair-level haplotypes or 2) have no variants in the region. This type is intended to identify regions where one or more inputs may have sequencing drop-out, leading to no variant calls. These regions do not "conflict", but are rather caused by technical artifacts. This type of output can be enabled with the `--enable-no-conflict` flag or an appropriate [merge strategy](#merge-strategies). The variants written to output are from the first VCF file with a non-empty set of variants in the region.
  * `majority` - Indicates that the majority (>50%) of the inputs have identical basepair-level haplotypes. Note that the majority may be an empty set of variants if `no_conflict` is not enabled. This type of output can be enabled with the `--enable-voting` flag or an appropriate [merge strategy](#merge-strategies). The variants written to output are from the first VCF file that is part of the majority.
  * `conflict_select` - Indicates that the region would have been labeled as `different`, except a specific input has been selected as priority when conflict occurs (via `--conflict-select`). This overrides the normal `different` output, and will output the variants from the specified input file, disregarding any conflicts from the other inputs.
  * `different` - Indicates that all other labels do not apply to the region. These variants are not included in the output merged VCF.

### Merge strategies
In `merge` mode, a merge strategy (`--merge-strategy`) can be specified that enables different type of labels to pass.
We outline each mode below:

* `exact` - The default mode, which only allows `identical` regions to pass merging.
* `no_conflict` - Enables `no_conflict` regions to pass merging. Same as using `--enable-no-conflict`.
* `majority` - Enables `majority` regions to pass merging. Same as using `--enable-voting`.
* `all` - Enables all optional region types to pass merging. Same as using `--enable-no-conflict --enable-voting`.

# Comparison types
Aardvark was designed to optimize comparisons on a sequence (or haplotype) level.
As such, the `BASEPAIR` comparison type is the primary output of the tool.
However, it also provides haplotype and genotype comparisons that are similar to those of other tools.
We provide details on the trade-offs of each below.

## Basepair
The basepair (`BASEPAIR`) type measures the total number of bases that are accurately changed/inserted/deleted in the query relative to truth.
In contrast to the other compare types, the concept of a "variant" is completely removed, and all comparisons are based on the edit distance between sequences.
Thus, this approach inherently gives higher weight to changes that impact multiple bases (i.e., insertions and deletions).
For example, an SNV will always have a weight of 1 (single basepair change) whereas an insertion/deletion of length `L` will have a weight of `L`.
Homozygous variants also have an implicit higher weight because they impact are scored on each haplotype sequence.
Additionally, since the concept of a "variant" is removed entirely, the number of truth TP bases and query TP bases will always be identical.
This is in constrast to GT and HAP, where the truth and query inputs likely have a different number of variants of each type.

An interesting characteristic of the `BASEPAIR` type is that some bases can be scored as partially correct (e.g., TP=0.5, FN=0.5).
To avoid floating-point calculations, all values reported by aardvark are doubled.
For example, if Aardvark reports TP=5, FP=2, and FN=1; then 2.5 bases were detected as true positive, 1 base was detected as a false positive, and 0.5 bases were detected as false negatives.

## Haplotype
The haplotype (`HAP`) type measures the total number of ALT variant alleles that are accurately detected.
For this type, every ALT allele in the truth is labeled as either TP or FN, and every ALT allele in the query is labeled as either TP or FP.
This labeling occurs on each of the optimized haplotype sequences (i.e., twice), generating a higher implicit weight for homozygous variants relative to heterozygous variants.
Additionally, greater weight is given to variant types with more representation, which may skew the results when indels are split into multiple entries.
Critically, recall is measured solely by the labels of the truth variants and precision solely by the labels of query variants.
This means that truth.TP and query.TP may not match if variants are represented differently or have different zygosities.

## Weighted Haplotype
The weighted haplotype (`WEIGHTED_HAP`) metric is the same as the haplotype (`HAP`) metric except each allele is weighted by the edit distance between the REF and ALT sequences, which effectively provides higher weight to allele changes that impact more basepairs.
For SNVs, this weight is always 1, so the `HAP` and `WEIGHTED_HAP` scores are identical.
For indel variant types, the weight is typically the length change relative to the reference allele.
Thus, a 20 basepair insertion will have 20x the weight of a 1 basepair insertion.
Conceptually, this is quite similar to the `BASEPAIR` metric, but without allowing for partial-credit matches.

## Genotype
The genotype (`GT`) type checks if the total zygosity of ALT variant alleles are accurately detected.
This type is conceptually the same as `HAP`, except homozygous variants are collapsed into a single TP, FN, or FP label.
Similar to `HAP`, this method gives greater weight to variant types with more representation, which may skew the results when indels are split into multiple entries.
We note that labels in truth may not match those in query if the representations or zygosities are different.
For example, if a truth genotype is heterozygous (0/1) and the query genotype is homozygous (1/1), then it will receive both truth.TP and query.FP labels, indicating the the one ALT allele in truth was detected but that the extra ALT allele in query was not matched in truth.

Despite these complexities, this mode does generate a single classification (TP, FP, or FN) for each variant in the input VCF files.
As such, Aardvark uses these `GT` classification to determine which VCF to save each variant in our [debug VCF outputs](./compare.md#debug-vcf-details).
