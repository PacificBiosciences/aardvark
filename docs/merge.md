# Merge guide
The `aardvark merge` command provides a method to merge variant calls.
This functionality is helpful when multiple variant callers or technologies are used for the same sample.
Merging variant sets can be complex, and Aardvark simplifies this process by modeling the variants as haplotypes.
This allows for different merge strategies for resolving conflicting variant sets. 

Table of contents:

* [Quickstart](#quickstart)
* [Output files](#output-files)
  * [Summary statistics file](#summary-statistics-file)
  * [Debug folder](#debug-folder)
* [Frequently asked questions](#frequently-asked-questions)

# Quickstart
```bash
aardvark merge \
    --reference <FASTA> \
    --input-vcf <VCF> \
    --input-vcf <VCF> \
    --regions <BED> \
    --output-vcfs <DIR>
```

Required parameters:
* `--reference <FASTA>` - Input reference genome FASTA file, gzip allowed
* `--input-vcf <VCF>` - Input variant call set, which must be tabix-indexed. The `merge` command expects to see this flag multiple times, once for each VCF file to be merged.
* `--regions <BED>` - Input regions to perform merge in. Any variants from any input VCF that are not **fully-contained** within these regions will be excluded from the entire analysis.
* `--output-vcfs <DIR>` - [Output merged VCF files](#output-vcf-folder) and supporting indices

## Additional options
* `--output-summary <TSV>` - Creates a [merge summary file](#merge-summary-file) containing summary metrics on how many variants were in agreement from each source file.
* `--merge-strategy <STRAT>` - Sets the merge strategy for the variants, which must be one of `exact`, `no_conflict`, `majority`, or `all`. See [Methods](./methods.md#merge-strategies) for details on each merge strategy.
* `--conflict-select <INDEX>` - By default, any conflicting regions are not written to the output VCF and are marked as failed regions. This option will enable an input VCF to be selected instead of generating a failed region. The expected input is a 0-based index (e.g., "0" selects the first VCF input, "1" the second, etc.).
* `--threads <THREADS>` - The number of compute threads to use for the comparison step
* `--output-debug <DIR>` - Creates a [debug folder](#debug-folder) and populates it with additional debug outputs
* `--vcf-sample <SAMPLE>` - Allows for the specification of the sample name in the provided VCF, which is usually only required for multi-sample VCFs. These must be specified in the same order as the `--input-vcf` option.
* `--vcf-tag <SAMPLE>` - A optional tag that can be specified for variants that successfully merge. These must be specified in the same order as `--input-vcf`. Default is "vcf_#".

# Output files
## Output VCF folder
The output VCF folder is the core output of the `aardvark merge` process.
It currently contains the following files:

* `passing.vcf.gz` - The merged VCF file containing any variants that pass the provided merge strategy. Details on this file can be found [below](#passing-vcf-file).
* `regions.bed.gz` - The passing regions for the merged variants, which are labeled with `{merge_reason}_{region_id}`. See [Methods](./methods.md#merge-command) for details on each passing reason. The `merge_reason` and `region_id` values correspond to the `MR` and `RI` tags, respectively, of the [passing VCF file](#passing-vcf-file).
* `failed_regions.bed.gz` - The failed regions for the merged variants, which are labeled with `different_{region_id}`.

### Passing VCF file
The passing VCF file is constructed using variants from all input VCFs.
These inputs may have different representations for a collection of variants in a region.
When a region is determined to be passing, the highest-priority (earliest) input VCF is that is part of the passing set is selected to have its variants copied.
For example, on `exact` mode, all variants must match across all inputs to pass.
In this situation, all output variants will *always* be from the first provided VCF.
In contrast, `no_conflict` mode may allow variants from a subset of the inputs. 
If a passing set is found that includes the second VCF file but not the first, then the second VCF file variants will be written to the output.

The header of this VCF is based off the first provided VCF file, but it is adjusted to add the Aardvark version (`aardvark_version`), Aardvark command that was used (`aardvark_command`), and any new VCF fields.
Variants inside the VCF have been reformatted to match the internal aardvark representation.
Typically, this just means that multi-allelic sites have been split into multiple rows.
Additionally, Aardvark reports the following fields for each variant:
* `INPUT/SOURCES` - The sources that contained the matching set of variants for this region. This can be modified using the `--vcf-tag` parameter.
* `INPUT/MR` - The Merge Reason indicating why this variant (and the region) was allowed into the passing VCF. This will be one of `identical`, `no_conflict`, or `majority`. See [Methods](./methods.md#merge-command) for details on each merge reason.
* `FORMAT/GT` - The GenoType from the input VCF file
* `FORMAT/RI` - The Region Id the variant was grouped into for calculating all values. The `MR` and `RI` tags will corresponding to the labels in the [passing regions BED file](#output-vcf-folder).

Partial example:
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
chr1	38232	.	A	G	.	.	SOURCES=pb,ont;MR=no_conflict	GT:RI	1/1:0
chr1	38907	.	C	T	.	.	SOURCES=pb,ilmn,ont;MR=identical	GT:RI	1/1:1
chr1	40244	.	C	T	.	.	SOURCES=pb,ont;MR=no_conflict	GT:RI	1/1:2
chr1	105279	.	G	C	.	.	SOURCES=ilmn;MR=no_conflict	GT:RI	1/1:3
chr1	106534	.	A	G	.	.	SOURCES=pb,ilmn,ont;MR=identical	GT:RI	1/1:4
chr1	106544	.	C	G	.	.	SOURCES=pb,ilmn,ont;MR=identical	GT:RI	1/1:4
chr1	116549	.	C	T	.	.	SOURCES=pb,ont;MR=majority	GT:RI	1/1:5
chr1	120458	.	T	C	.	.	SOURCES=pb,ilmn,ont;MR=identical	GT:RI	1/1:6
chr1	120705	.	T	C	.	.	SOURCES=pb,ilmn,ont;MR=identical	GT:RI	1/1:7
...
```

## Merge summary file
The merge summary file contains statistics on how many variants from each VCF file were included or excluded.
The outputs are grouped by the merge reason, variant type, and VCF.

### Fields
* `merge_reason` - The merge reason assigned to the collection of variants. The label here will match the `MR` tag from the VCF files, but may have additional index information provided for clarity. For example, the `majority_0_2` label indicates variants the variants from input 0 and input 2 were in agreement and are thus "passing", but that input 1 had a different set of variants that were instead "failed".
* `variant_type` - The variant type counted on the row
* `vcf_index` - The 0-based index of the VCF counted on the row
* `vcf_label` - The corresponding label (`--vcf-tag`) for the VCF
* `pass_variants` - The total number of passing variants with the given `merge_reason`, `variant_type`, and `vcf_index`. "Passing" indicates that the variant is part of a region that is basepair-equivalent to the region in the merged VCF. It **does not** mean that the exact representation is in the merged VCF. By definition, the `different` merge reason will never have passing variants.
* `fail_variants` - The total number of failing variants with the given `merge_reason`, `variant_type`, and `vcf_index`. "Failing" indicates that the variant is part of a region that is *not* basepair-equivalent to the region in the merged VCF. By definition, the `identical` and `no_conflict` merge reasons will never have failing variants.

### Example
This partial example only contains the "Snv" type for brevity:
```
merge_reason	variant_type	vcf_index	vcf_label	pass_variants	fail_variants
different	Snv	0	pb	0	88072
different	Snv	1	ilmn	0	51250
different	Snv	2	ont	0	75279
no_conflict_0	Snv	0	pb	26207	0
no_conflict_0_1	Snv	0	pb	144889	0
no_conflict_0_1	Snv	1	ilmn	144424	0
no_conflict_0_2	Snv	0	pb	131610	0
no_conflict_0_2	Snv	2	ont	131575	0
no_conflict_1	Snv	1	ilmn	19896	0
no_conflict_1_2	Snv	1	ilmn	13295	0
no_conflict_1_2	Snv	2	ont	13290	0
no_conflict_2	Snv	2	ont	44727	0
majority_0_1	Snv	0	pb	153882	0
majority_0_1	Snv	1	ilmn	150677	0
majority_0_1	Snv	2	ont	0	111297
majority_0_2	Snv	0	pb	40561	0
majority_0_2	Snv	1	ilmn	0	37810
majority_0_2	Snv	2	ont	40529	0
majority_1_2	Snv	0	pb	0	16904
majority_1_2	Snv	1	ilmn	17900	0
majority_1_2	Snv	2	ont	18133	0
identical	Snv	0	pb	3135228	0
identical	Snv	1	ilmn	3132738	0
identical	Snv	2	ont	3134498	0
```

## Debug folder
The debug folder (`--debug-folder`) contains many files that may be useful for debugging errors from aardvark, or for getting greater details around specific comparisons that were performed.
The following files are currently created when this option is specified:

* `cli_settings.json` - JSON dump of the exact parameters that were used to run Aardvark

# Frequently asked questions
Will populate as needed.
