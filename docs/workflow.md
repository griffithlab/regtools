# RegTools example workflow

This is an example workflow for running RegTools on a cohort of samples. This analysis requires that there be a VCF and RNA bam file for each sample. The workflow described below was used to run our own analysis on TCGA data.

By the end of the analysis, the directory structure should look like the example below. The `*` in the example below refers to the tag/parameter used to run `regtools cis-splice-effects identify` with. At the bottom of this page, we provide a description of each of these files.

```bash
- Project/
  - all_splicing_variants*.bed
  - dir_names.tsv
  - variants_all_sorted.vcf.gz
  - variants_all_sorted.vcf.gz.tbi
  - samples/
    - Sample_1/
      - tumor_rna_alignments.bam
      - tumor_rna_alignments.bam.bai
      - variants.vcf.gz
      - variants.vcf.gz.tbi
      - logs/
      - output/
        - cse_identify_filtered_*
        - cse_identify_filtered_compare_*
        - variants*.bed
    - Sample_2/
      - tumor_rna_alignments.bam
      - tumor_rna_alignments.bam.bai
      - variants.vcf.gz
      - variants.vcf.gz.tbi
      - logs/
      - output/
        - cse_identify_filtered_*
        - cse_identify_filtered_compare_*
        - variants*.bed
  - compare_junctions/
    - hist/
      - junction_pvalues_*.tsv
```

## RegTools commands

**Set tag and parameter shell arguments**

```bash
tag=<tag>
param=<run option>
# (e.g. tag=default param=""; tag=E param="-E"; tag=i20e5 param="-i 20 -e 5")
```

**Run `regtools cis-splice-effects identify` with desired options for selecting variant and window size**

```bash
for i in samples/*/; do regtools cis-splice-effects identify $param -o ${i}/output/cse_identify_filtered_$tag.tsv -j ${i}/output/cse_identify_filtered_$tag.bed -v ${i}/output/cse_identify_filtered_$tag.vcf ${i}/variants.per_gene.vep.vcf.gz ${i}/tumor_rna_alignments.bam /reference.fa reference.gtf; done
```

**Make `variant.bed` for each sample**

```bash
for i in samples/*/; do bash variants.sh ${i}/output/cse_identify_filtered_$tag.tsv ${i}/output/variants_$tag.bed; done
```

**Combine each sample's `variant.bed` file per tag to get all variants that were deemed significant to splicing across all samples for a given tag**

```bash
echo -e 'chrom\tstart\tend\tsamples' > all_splicing_variants_$tag.bed
for i in samples/*/; do j=${i##samples/}; uniq ${i}output/variants_$tag.bed | awk -v var=${j%%/} '{print $0 "\t" var}' >> all_splicing_variants_$tag.bed; done
```

**Make vcf of all variants across all samples (from each sample's variants.vcf). Then, compress it and index it.**

```bash
vcf-concat samples/*/variants.vcf.gz | vcf-sort > all_variants_sorted.vcf

bgzip all_variants_sorted.vcf

tabix all_variants_sorted.vcf.gz
```

**Run `regtools cis-splice effects identify` on all samples with all variants (with `$tag` options as example)**

```bash
for i in samples/*/; do bsub -oo $i/logs/regtools_compare_$tag.lsf regtools cis-splice-effects identify $param -o ${i}/output/cse_identify_filtered_compare_$tag.tsv -j ${i}/output/cse_identify_filtered_compare_$tag.bed -v ${i}/output/cse_identify_filtered_compare_$tag.vcf all_variants_sorted.vcf.gz ${i}/tumor_rna_alignments.bam reference.fa reference.gtf; done
```

## Statistical analysis

**Make directory to store comparison results**

```bash
mkdir -p compare_junctions/hist
```

**Run `compare_junctions_hist.py` on sample data**

This script compares RegTools associations across all samples in a defined cohort to determine significance. By default, this script divides the input data into chunks of 50,000 variants in order to not max out memory for large cohorts. We provide three options for how to group case vs control groups with respect to a particular event (i.e. which samples do you consider to have the variant of interest vs. those that do not). The options are as follows:

1) only samples with exactly the same variant will be included in the case group, all others samples (even potentially those with very similar variants) will be included with the control group (`--variant-grouping=strict`)
2) similar variants will be excluded from both the case and control groups (`--variant-grouping=exclude`)
3) similar variants will be included in the case group (`--variant-grouping=include`)

*NOTE: "similar" is defined as associated with the same junction*

We recommend using `--variant-grouping=exclude` since determining whether nearby variants have a similar effect in terms of splicing consequence is often difficult.

```bash
python3 compare_junctions_hist.py --tag <tag> --variant-grouping <parameter> 
```

**Run `filter_and_BH.R` to adjust p values and filter results**

This script filters significant results on junction support and junction type and then performs the Benjamini-Hochberg Procedure to apply multiple test correction to the data.

```bash
Rscript --vanilla filter_and_BH.R <tag>
```

## File description

* **`all_splicing_variants*.bed`** - a file containing all of the variants that regtools identified as being associated with a junction for the given parameters used to run `cis-splice-effects identify`.
* **`dir_names.tsv`** - a file containing a list of each of the sample directories with each directory on a new line. This can be obtained by using `ls samples/ > dir_names.tsv`. For this example, it would look like:

  ```bash
  Sample_1
  Sample_2
  ```

* **`variants_all_sorted.vcf.gz`** - a compressed vcf file containing all variants from all samples.
* **`variants_all_sorted.vcf.gz.tbi`** - an index file for the vcf file mentioned above.
* **`samples/`** - a directory containing each of the samples to be analyzed alongside each other.
* **`Sample_1/`** - a sample directory. This will contain input data files as well as output files from RegTools.
* **`tumor_rna_alignments.bam`** - file containing aligned RNA-seq reads for the given sample.
* **`tumor_rna_alignments.bam.bai`** - index file for the above RNA-seq alignment file.
* **`variants.vcf.gz`** - a compressed vcf file containing all variants from a given samples.
* **`variants.vcf.gz.tbi`** - an index file for the vcf file mentioned above.
* **`logs/`** - directory containing log or error files for a given sample.
* **`output/`** - directory containing RegTools output files for a given sample.
* **`cse_identify_filtered_*`** - RegTools output files from the initial RegTools run for a given sample. This will contain results for this sample's variants only.
* **`cse_identify_filtered_compare_*`** - RegTools output files from the second RegTools run for a given sample. This will contain results for all samples' variants.
* **`variants*.bed`** - a bedfile containing the variants considered to be splicing relevant for a given RegTools parameter. This is used later to make `all_splicing_variants*.bed`.
* **`compare_junctions/hist/`** - directory containing output from the statistics script analyze all variants across all samples.
* **`junction_pvalues_*.tsv`** - a file containing the output from the statistic analysis script.
