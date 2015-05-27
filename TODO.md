###Task 1
- Implement a tools that starts with an RNA-seq BAM and produces an annotated junctions file
- Start with an RNA-seq BAM and produce a junctions.bed file (use code from TopHat)
- Then reimplement the functionality of this script: genome/lib/perl/Genome/Model/ClinSeq/OriginalScripts/rnaseq/annotateObservedJunctions.pl
- Note that this script was reimplemented in the GMS here:
  genome/lib/perl/Genome/Model/Tools/Transcriptome/SpliceJunctionSummary.pm

- Example input: /gscmnt/gc13001/info/model_data/a0517e48eafd4c189fc0d4e105b32fe2/buildd6edd280e32b45dc80ff904053a22be1/AML103/rnaseq/tumor/tophat_junctions_absolute/AlignmentResult_69b9e9c7887f40db825b64b99a3b7e03_junctions.bed
- Example output: /gscmnt/gc13001/info/model_data/a0517e48eafd4c189fc0d4e105b32fe2/buildd6edd280e32b45dc80ff904053a22be1/AML103/rnaseq/tumor/tophat_junctions_absolute/observed.junctions.anno.NCBI-human.ensembl-67_37l_v2.tsv
- As a starting point, lets create a tool: 'regtools annotate-junctions' that creates an file with no diffs to that output file
- Other input files will be the reference genome fasta, and the ensembl GTF file (same version)

- Test case to be reproduced (from within a TGI system where the GMS is installed):

```bash  
git clone git@github.com:griffithlab/regtools.git
cd regtools/tests
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/brain_vs_uhr_w_ercc/downsampled_5pc_chr22/chr22.fa.gz
wget https://xfer.genome.wustl.edu/gxfer1/project/gms/testdata/bams/brain_vs_uhr_w_ercc/downsampled_5pc_chr22/genes_chr22.gtf.gz
gunzip *.gz
rm -fr /tmp/junction_summary/
mkdir /tmp/junction_summary/
gmt transcriptome splice-junction-summary --output-directory='/tmp/junction_summary/' --observed-junctions-bed12-file='junctions.chr22.bed' --reference-fasta-file='chr22.fa' --annotation-gtf-file='genes_chr22.gtf' --annotation-name='Ensembl'
cd /tmp/junction_summary/
rm -fr SpliceJunctionSummary.R.stderr summary SpliceJunctionSummary.R.stdout Ensembl.Junction.TranscriptExpression.top1percent.tsv Ensembl.Junction.GeneExpression.top1percent.tsv
```

