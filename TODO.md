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

- Test case to be reproduced:
gmt transcriptome splice-junction-summary --output-directory='/gscuser/mgriffit/tmp/' --observed-junctions-bed12-file='/gscmnt/gc12001/info/model_data/2887640332/build127105782/alignments/junctions.bed' --reference-fasta-file='/gscmnt/sata420/info/model_data/2857786885/build102671028/all_sequences.fa' --annotation-gtf-file='/gscmnt/gc12001/info/model_data/2772828715/buildd00a39c84382427fa0efdec3229e8f5f/annotation_data/rna_annotation/106942997-all_sequences.gtf' --annotation-name='Ensembl'

