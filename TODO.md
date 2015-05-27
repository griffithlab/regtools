###Task 1
- Implement a tool that starts with an RNA-seq BAM and produces an annotated junctions file
 
- Start with an RNA-seq BAM and produce a junctions.bed file (use code from TopHat)

- Then reimplement the functionality of this GMT in the GMS:

genome/lib/perl/Genome/Model/Tools/Transcriptome/SpliceJunctionSummary.pm

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

