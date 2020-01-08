# RegTools example workflow

This is an example workflow for running RegTools on a cohort of samples. This analysis requires that there be a vcf and RNA bam/cram file for each samples. The outline described below was used to run our own analysis on TCGA data.

By the end of the analysis, the directory structure should look like the example below. The `*` in the example below refers to the tag/parameter used to run `regtools cis-splice-effects identify` with.

```bash
- Project/
  - all_splicing_variants*.bed
  - paths.tsv
  - make_vcfs.sh
  - dir_names.tsv
  - variants_all_sorted.vcf.gz
  - variants_all_sorted.vcf.gz.tbi
  - samples/
    - Sample_1/
      - tumor_rna_alignments.bam
      - tumor_rna_alignments.bam.bai
      - variants.per_gene.vep.vcf.gz
      - variants.per_gene.vep.vcf.gz.tbi
      - variants.ensembl
      - logs/
      - output/
        - cse_identify_filtered_*
        - cse_identify_filtered_compare_*
        - variants*.bed
    - Sample_2/
      - tumor_rna_alignments.bam
      - tumor_rna_alignments.bam.bai
      - variants.per_gene.vep.vcf.gz
      - variants.per_gene.vep.vcf.gz.tbi
      - variants.ensembl
      - logs/
      - output/
        - cse_identify_filtered_*
        - cse_identify_filtered_compare_*
        - variants*.bed
  - compare_junctions/
    - hist/
      - junction_pvalues_*.tsv
```

### Set tag and parameter shell arguments

```bash
tag=<tag>
param=<run option>
# (e.g. tag=default param=""; tag=E param="-E"; tag=i20e5 param="-i 20 -e 5")
```

### Run `regtools cis-splice-effects identify` with desired options for selecting variant and window size

```bash
for i in samples/*/; do regtools cis-splice-effects identify $param -o ${i}/output/cse_identify_filtered_$tag.tsv -j ${i}/output/cse_identify_filtered_$tag.bed -v ${i}/output/cse_identify_filtered_$tag.vcf ${i}/variants.per_gene.vep.vcf.gz ${i}/tumor_rna_alignments.bam /reference.fa reference.gtf; done
```

### Make `variant.bed` for each sample

```bash
for i in samples/*/; do bash variants.sh ${i}/output/cse_identify_filtered_$tag.tsv ${i}/output/variants_$tag.bed; done
```

### Combine each sample's `variant.bed` file per tag to get all variants that were deemed significant to splicing across all samples for a given tag

```bash
echo -e 'chrom\tstart\tend\tsamples' > all_splicing_variants_$tag.bed
for i in samples/*/; do j=${i##samples/}; uniq ${i}output/variants_$tag.bed | awk -v var=${j%%/} '{print $0 "\t" var}' >> all_splicing_variants_$tag.bed; done
```

### Make vcf of all variants across all samples (from each sample's variants.vcf). Then, compress it and index it

```bash
vcf-concat samples/*/variants.vcf.gz | vcf-sort > all_variants_sorted.vcf

###### optional ######
bgzip all_variants_sorted.vcf

tabix all_variants_sorted.vcf.gz
```

### Run `regtools cis-splice effects identify` on all samples with all variants (with `$tag` options as example)

```bash
for i in samples/*/; do bsub -oo $i/logs/regtools_compare_$tag.lsf regtools cis-splice-effects identify $param -o ${i}/output/cse_identify_filtered_compare_$tag.tsv -j ${i}/output/cse_identify_filtered_compare_$tag.bed -v ${i}/output/cse_identify_filtered_compare_$tag.vcf all_variants_sorted.vcf.gz ${i}/tumor_rna_alignments.bam reference.fa reference.gtf; done
```

## Beginning of compare junctions analysis

### Make directory to store comparison results

```bash
mkdir -p compare_junctions/hist
```

### Run `compare_junctions_hist.R` on sample data

```bash
Rscript --vanilla compare_junctions_hist.R <tag>
```

### Run `filter_and_BH.R` to adjust p values and filter results

```bash
Rscript --vanilla filter_and_BH.R <tag>
```
