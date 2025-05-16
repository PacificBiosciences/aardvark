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

## Debug folder
The debug folder (`--debug-folder`) contains many files that may be useful for debugging errors from aardvark, or for getting greater details around specific comparisons that were performed.
The following files are currently created when this option is specified:

* `cli_settings.json` - JSON dump of the exact parameters that were used to run Aardvark

# Frequently asked questions
Will populate as needed.
