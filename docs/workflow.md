# RegTools example workflow

This is an example workflow for running RegTools on a cohort of samples. This analysis requires that there be a vcf and RNA bam/cram file for each samples. The outline described below was used to run our own analysis on TCGA data.

By the end of the analysis, the directory structure should look like so:

- Project/ (SCLC/)
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
  - compare_junctions/
    - hist/
      - junction_pvalues_*.tsv

## Set tag and parameter shell arguments

tag=<tag> param=<run option> (e.g. tag=default param=""; tag=E param="-E"; tag=i20e5 param="-i 20 -e 5")

# run regtools cse identify with desired options for selecting variant and window size
for i in samples/*/; do bsub -oo $i/logs/regtools_actual_$tag.lsf regtools cis-splice-effects identify $param -o ${i}/output/cse_identify_filtered_$tag.tsv -j ${i}/output/cse_identify_filtered_$tag.bed -v ${i}/output/cse_identify_filtered_$tag.vcf ${i}/variants.per_gene.vep.vcf.gz ${i}/tumor_rna_alignments.bam /gscmnt/gc2602/griffithlab/regtools/yafeng/GRCh37.fa /gscmnt/gc2602/griffithlab/regtools/GRCh37.87.exons.sorted.gtf; done

# make variant.bed
for i in samples/*/; do bsub -oo $i/logs/make_variant_bed_$tag.lsf bash /gscmnt/gc2602/griffithlab/regtools/yafeng/scripts/variants.sh ${i}/output/cse_identify_filtered_$tag.tsv ${i}/output/variants_$tag.bed; done 

# make bed (really just tsv with columns: chrom start end samples) with all variants that were deemed significant to splicing across all samples
echo -e 'chrom\tstart\tend\tsamples' > all_splicing_variants_$tag.bed
for i in samples/*/; do j=${i##samples/}; uniq ${i}output/variants_$tag.bed | awk -v var=${j%%/} '{print $0 "\t" var}' >> all_splicing_variants_$tag.bed; done

# make vcf of all variants across all samples (from variants.vcf)
vcf-concat samples/*/variants.per_gene.vep.vcf.gz | vcf-sort > all_variants_sorted.vcf
bgzip all_variants_sorted.vcf 
tabix all_variants_sorted.vcf.gz 

# run cis-splice effects identify on all samples with all variants (with $tag options as example)
for i in samples/*/; do bsub -oo $i/logs/regtools_compare_$tag.lsf regtools cis-splice-effects identify $param -o ${i}/output/cse_identify_filtered_compare_$tag.tsv -j ${i}/output/cse_identify_filtered_compare_$tag.bed -v ${i}/output/cse_identify_filtered_compare_$tag.vcf variants_all_sorted.vcf.gz ${i}/tumor_rna_alignments.bam /gscmnt/gc2602/griffithlab/regtools/yafeng/GRCh37.fa /gscmnt/gc2602/griffithlab/regtools/GRCh37.87.exons.sorted.gtf; done

<===== JUNCTION COMPARISON ANALYSIS =====>

# make directories
mkdir compare_junctions/
mkdir compare_junctions/outlier
mkdir compare_junctions/ratio

# outlier analysis
bsub -M 650000000 -R 'select[mem>65000] span[hosts=1] rusage[mem=65000]' /gscuser/zskidmor/R-3.3.0/bin/Rscript --vanilla /gscmnt/gc2602/griffithlab/regtools/yafeng/scripts/compare_junctions_outlier.R $tag

# ratio analysis
bsub -M 650000000 -R 'select[mem>65000] span[hosts=1] rusage[mem=65000]' /gscuser/zskidmor/R-3.3.0/bin/Rscript --vanilla /gscmnt/gc2602/griffithlab/regtools/yafeng/scripts/compare_junctions_ratio.R $tag

# vep comparison (outlier)
bsub /gscuser/zskidmor/R-3.3.0/bin/Rscript --vanilla /gscmnt/gc2602/griffithlab/regtools/yafeng/scripts/vep_compare.R outlier $tag

# vep comparison (ratio)
bsub /gscuser/zskidmor/R-3.3.0/bin/Rscript --vanilla /gscmnt/gc2602/griffithlab/regtools/yafeng/scripts/vep_compare.R ratio $tag

NOTE: in the vep comparison, since regtools doesn't really have a concept of "variants" per se but rather "variant positions" (doesn't care about the actual alleles), we are really counting (positions of variants found splicing-significant by vep AND found significant by regtools) / (positions of variants found significant by found significant by regtools)

To count:

echo -e "anchor\tvep_significant_vars\ttotal_significant_vars\tpercent_vep_significant" > comparison_counts.tsv
for i in default i20e5 i50e5 E I; do 
	echo -n $i >> comparison_counts.tsv; 
	echo -en "\t" >> comparison_counts.tsv; 
	vep=$(cut -f 1,2,7 vep_comparison_$i.tsv | sort | uniq | grep "TRUE" | wc -l)
	echo -n $vep >> comparison_counts.tsv; 
	echo -en "\t" >> comparison_counts.tsv;
	total=$(($(cut -f 1,2,7 vep_comparison_$i.tsv | sort | uniq | wc -l)-1)) 
	echo -n $total >> comparison_counts.tsv;
	echo -en "\t" >> comparison_counts.tsv;
	echo $(bc -l <<< "$vep/$total") >> comparison_counts.tsv;
done
