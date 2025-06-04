# Compare guide
The `aardvark compare` command is used to compare (or benchmark) a "query" call set against a "truth" call set. 
The `compare` command will load the query and truth variant sets, group them into smaller clusters that are interrogated independently, and then perform the comparisons of the sets.
The main output is a summary statistic file reporting  precision, recall, and F1 score measuring how similar and different the "query" and "truth" are.

Table of contents:

* [Quickstart](#quickstart)
* [Output files](#output-files)
  * [Summary statistics file](#summary-statistics-file)
  * [Labeled VCF files](#labeled-vcf-files)
  * [Debug folder](#debug-folder)
* [Complex inputs](#complex-inputs)
* [Frequently asked questions](#frequently-asked-questions)

# Quickstart
```bash
aardvark compare \
    --reference <FASTA> \
    --truth-vcf <VCF> \
    --query-vcf <VCF> \
    --confidence-regions <BED> \
    --output-dir <DIR>
```

Required parameters:
* `--reference <FASTA>` - Input reference genome FASTA file, gzip allowed
* `--truth-vcf <VCF>` - Input truth variant call set, which must be tabix-indexed
* `--query-vcf <VCF>` - Input query variant call set, which must be tabix-indexed
* `--confidence-regions <BED>` - Input regions to perform comparison in. Any variants from either truth or query that are not **fully-contained** within these regions will be excluded from the entire analysis.
* `--output-dir <DIR>` - Primary output folder containing the [summary statistics file](#summary-statistics-file) and [labeled VCF files](#labeled-vcf-files)

## Additional options
* `--threads <THREADS>` - The number of compute threads to use for the comparison step
* `--output-debug <DIR>` - Creates a [debug folder](#debug-folder) and populates it with more detailed statistics from the comparison
* `--truth-sample <SAMPLE>` and `--query-sample <SAMPLE>` - Allows for the specification of the truth or query sample name in the provided VCF. This is usually only needed when they are multi-sample VCFs.
* `--stratifications` - The root TSV file for a stratification. See [stratifications](#stratifications) for details.

# Output files
## Summary statistics file
The summary statistics file (`summary.tsv`) is a TSV output containing high-level summary metrics for the comparison including variant counts and distance metrics.

### Fields
* `compare_label` - A user-provided comparison label that is just passed through to this output. Specified via `--compare-label` option.
* `comparison` - The comparison type, which will be one of `GT` (genotype), `HAP` (haplotype), or `BASEPAIR` (sequence/basepair level). See [methods](./methods.md#comparison-types) for more details on each comparison type.
* `filter` - Indicates if any filter was applied.
* `region_label` - The region label from stratification inputs. By default, only `ALL` is provided which contains all variants analyzed. If [stratifications](#stratifications) are provided, then additional rows for each stratification label will be added and this column will contain the label.
* `variant_type` - The type of variant that the assessment corresponds to.
  * `Snv` - Single Nucleotide Variant, requires both REF and ALT to be exactly 1 bp
  * `Insertion` - Insertion variant, required REF to be exactly 1 bp and ALT to be >1 bp
  * `Deletion` - Deletion variant, requires REF to be >1 bp and Alt to be exactly 1 bp
  * `Indel` - Insertion/deletion variant, requires both REF and ALT to be >1 bp
  * `JointIndel` - Group category that includes `Insertion`, `Deletion`, and `Indel` types
  * `ALL` - Group category that includes all variants types
* `truth_total` - The total number of truth values that were assessed.
* `truth_tp` - The number of truth values that were also identified in the query set, indicating a true positive (TP).
* `truth_fn` - The number of truth values that were not identified in the query set, indicating a false negative (FN).
* `query_total` - The total number of query values that were assessed.
* `query_tp` - The number of query values that were matched to the truth set, indicating TP. For `comparison` mode `BASEPAIR`, this value should be equal to `truth_tp`. Other modes may be different due to variant representation.
* `query_fp` - The number of query values that were not matched to the truth set, indicating a false positive (FP).
* `metric_recall` - The recall metric, which is calculated as `truth_tp / truth_total`.
* `metric_precision` - The precision metric, which is calculated as `query_tp / (query_tp + query_fp)`.
* `metric_f1` - The F1 score metric, which is calculated as the harmonic mean of recall and precision. It is often used as an overall summary metric for comparisons.

### Example
This snippet of the summary file shows the `GT` and `BASEPAIR` comparison types, filtered to `ALL`, `Snv`, and `JointIndel` variant types.

```
compare_label	comparison	region_label	filter	variant_type	truth_total	truth_tp	truth_fn	query_total	query_tp	query_fp	metric_recall	metric_precision	metric_f1
compare	GT	ALL	ALL	ALL	3751311	3737771	13540	3757515	3741335	16180	0.9963905951812579	0.9956939626322183	0.9960421571004356
compare	GT	ALL	ALL	Snv	3256086	3250157	5929	3256468	3252130	4338	0.9981791021490218	0.9986678818892125	0.9984234321984007
compare	GT	ALL	ALL	JointIndel	495225	487614	7611	501047	489205	11842	0.9846312282295926	0.9763654906625526	0.980480939116689
compare	BASEPAIR	ALL	ALL	ALL	12661262	12624128	37134	12657970	12624128	33842	0.9970671170061879	0.997326427539329	0.9971967554150142
compare	BASEPAIR	ALL	ALL	Snv	9087802	9075500	12302	9087548	9080279	7269	0.9986463173383399	0.9992001142662466	0.9989231390468849
compare	BASEPAIR	ALL	ALL	JointIndel	3587078	3561871	25207	3585212	3562195	23017	0.9929728319261527	0.9935800170254925	0.9932763316834902
```

## Labeled VCF files
The output folder will also contain two VCF files with labeled variants, which are based on the genotype (GT) comparison mode.
Index files (.tbi) are automatically generated for each VCF.

* `truth.vcf.gz` - Contains variants from the truth call set that are labeled as true positives (TP) or false negatives (FN)
* `query.vcf.gz`- Contains variants from the query call set that are labeled as true positives (TP) or false positives (FP)

### VCF details
Each VCF file is constructed using either the truth or query VCF as a source.
The header of each VCF is adjusted to add the Aardvark version (`aardvark_version`), Aardvark command that was used (`aardvark_command`), and any new VCF fields.
Variants inside the VCF have been reformatted to match the internal aardvark representation.
Typically, this just means that multi-allelic sites have been split into multiple entries.
Additionally, Aardvark reports the following fields for each variant:

* `GT` - The GenoType from the input VCF file
* `BD` - The Benchmark Decision for the variant, which will be one of TP (true positive), FP (false positive), or FN (false negative)
* `EA` - The Expected Allele count, which is 2 for a homozygous variant and 1 for a heterozygous variant
* `OA` - The Observed Allele count. All entries in a TP file will have `EA==OA`, all entries in a FN file will have `EA > OA`, and all entries in a FP file will have `EA < OA`.
* `RI` - The Region Id the variant was grouped into for calculating all values. These will correspond to the `region_id` of the [region summary file](#region-summary-file).

### Example
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG001
chr1	783006	.	A	G	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:0
chr1	783175	.	T	C	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:1
chr1	784860	.	T	C	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:2
chr1	785417	.	G	A	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:3
chr1	797392	.	G	A	.	.	.	GT:BD:EA:OA:RI	0/1:TP:1:1:4
chr1	798618	.	C	T	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:5
chr1	798662	.	G	A	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:5
chr1	800046	.	G	A	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:6
chr1	801142	.	A	T	.	.	.	GT:BD:EA:OA:RI	1/1:TP:2:2:7
...
```

## Debug folder
The debug folder (`--debug-folder`) contains many files that may be useful for debugging errors from aardvark, or for getting greater details around specific comparisons that were performed.
The following files are currently created when this option is specified:

* `cli_settings.json` - JSON dump of the exact parameters that were used to run Aardvark
* `region_sequences.tsv.gz` - A [region-specific sequence file](#region-sequence-file) containing the Aardvark-constructed haplotype sequences for each sub-problem that was created as part of the Aardvark process. This file is gzip-compressed to reduce storage.
* `region_summary.tsv.gz` - A [region-specific summary file](#region-summary-file), containing summary statistics for each sub-problem that was created as part of the Aardvark process. This file is gzip-compressed to reduce storage.

### Region sequence file
This file captures the constructed haplotype sequences for each region.

Fields:
* `region_id` - A unique region identifier corresponding to the `RI` field of the [debug VCFs](#debug-vcf-details)
* `coordinates` - The coordinates of the corresponding region
* `ref_seq` - The reference sequence within the region, which is copied from the provided reference file
* `truth_seq1` - The first constructed truth haplotype sequence
* `truth_seq2` - The second constructed truth haplotype sequence
* `query_seq1` - The first constructed query haplotype sequence. In an exact match, this is expected to be identical to `truth_seq1`.
* `query_seq2` - The second constructed query haplotype sequence. In an exact match, this is expected to be identical to `truth_seq2`.

Partial example:
```
region_id	coordinates	ref_seq	truth_seq1	truth_seq2	query_seq1	query_seq2
0	chr1:782956-783056	TAATTTTTTATATTGATTGTATACTGCAGTGATAATATTTTGGATGTATCAGGTTAAATAAAATTGACTGATTTCACCTTTTTCCTATTTTAAAAGTGGCT	TAATTTTTTATATTGATTGTATACTGCAGTGATAATATTTTGGATGTATCGGGTTAAATAAAATTGACTGATTTCACCTTTTTCCTATTTTAAAAGTGGCT	TAATTTTTTATATTGATTGTATACTGCAGTGATAATATTTTGGATGTATCGGGTTAAATAAAATTGACTGATTTCACCTTTTTCCTATTTTAAAAGTGGCT	TAATTTTTTATATTGATTGTATACTGCAGTGATAATATTTTGGATGTATCGGGTTAAATAAAATTGACTGATTTCACCTTTTTCCTATTTTAAAAGTGGCT	TAATTTTTTATATTGATTGTATACTGCAGTGATAATATTTTGGATGTATCGGGTTAAATAAAATTGACTGATTTCACCTTTTTCCTATTTTAAAAGTGGCT
1	chr1:783125-783225	GTTGATGAAAAATATTGTTGGTGAGCTCTGCTTAGGTAATATATAGGACATGAGCAGAGAGGAGGCACGTGAACAGTTCTGGCCTGGAGTAGGCTTCATTG	GTTGATGAAAAATATTGTTGGTGAGCTCTGCTTAGGTAATATATAGGACACGAGCAGAGAGGAGGCACGTGAACAGTTCTGGCCTGGAGTAGGCTTCATTG	GTTGATGAAAAATATTGTTGGTGAGCTCTGCTTAGGTAATATATAGGACACGAGCAGAGAGGAGGCACGTGAACAGTTCTGGCCTGGAGTAGGCTTCATTG	GTTGATGAAAAATATTGTTGGTGAGCTCTGCTTAGGTAATATATAGGACACGAGCAGAGAGGAGGCACGTGAACAGTTCTGGCCTGGAGTAGGCTTCATTG	GTTGATGAAAAATATTGTTGGTGAGCTCTGCTTAGGTAATATATAGGACACGAGCAGAGAGGAGGCACGTGAACAGTTCTGGCCTGGAGTAGGCTTCATTG
...
```

### Region summary file
This file captures summary metrics for each discrete region, similar to the overall summary statistics file.

Fields:
* `region_id` - A unique region identifier corresponding to the `RI` field of the [debug VCFs](#debug-vcf-details)
* `coordinates` - The coordinates of the corresponding region
* `comparison` - The comparison type for the metrics on the row, which will be one of `GT`, `HAP`, or `BASEPAIR`. This is a summary for `ALL` variant types.
* `truth_total`, `truth_tp`, `truth_fn`, `query_tp`, `query_fp`, `metric_recall`, `metric_precision`, `metric_f1` - See definitions for the [summary statistics file](#summary-statistics-file)

Partial example:
```
region_id	coordinates	comparison	truth_total	truth_tp	truth_fn	query_tp	query_fp	metric_recall	metric_precision	metric_f1
0	chr1:782956-783056	GT	1	1	0	1	0	1.0	1.0	1.0
0	chr1:782956-783056	HAP	2	2	0	2	0	1.0	1.0	1.0
0	chr1:782956-783056	BASEPAIR	4	4	0	4	0	1.0	1.0	1.0
1	chr1:783125-783225	GT	1	1	0	1	0	1.0	1.0	1.0
1	chr1:783125-783225	HAP	2	2	0	2	0	1.0	1.0	1.0
1	chr1:783125-783225	BASEPAIR	4	4	0	4	0	1.0	1.0	1.0
...
```

# Complex inputs
## Stratifications
Sometimes it is useful to stratify the results into known regions for further analysis.
Aardvark accepts the stratification format [supported by Hap.py](https://github.com/Illumina/hap.py/blob/master/doc/happy.md#stratification-via-bed-regions) and distributed by [Genome in a Bottle](https://github.com/genome-in-a-bottle/genome-stratifications).
In short, the `--stratifications` parameter accepts the "root" TSV from the stratification, where each row contains a label and corresponding BED file.
This file-of-files may contain paths relative to the root TSV, and those BED files may be gzip-compressed.
A short example is in our [test files](../test_data/example_stratification/strat.tsv), and a real example is available through [GIAB](https://github.com/genome-in-a-bottle/genome-stratifications/blob/master/GRCh38/v3.1-GRCh38-all-stratifications.tsv).

Aardvark takes a haplotype-centric approach to all analysis, including stratification.
This means Aardvark compares the entire region against the stratifications when determining labels.
Most variants in a typical benchmark are isolated, so for those variants Aardvark will just check if the reference bases are fully-contained within the stratification regions.
If two or more variants are located in the same Aardvark region, it instead checks if the bases from the start of the first variant through the end of the last variant are fully-contained.
As a result, stratification regions that are small and co-located with increased variation may be under-reported in the stratifications.

# Frequently asked questions
## Which comparison mode is most like existing benchmark tools?
Other tools, like [hap.py](https://github.com/Illumina/hap.py) or [rtg vcfeval](https://github.com/RealTimeGenomics/rtg-tools), typically use the genotype accuracy for benchmarking.
Aardvark's GenoType (`GT`) comparison type is most similar to these other approaches, requiring that the alternate sequences exactly match while not allowing for partial credit.
