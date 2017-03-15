###Synopsis
The `cis-ase identify` command is used to identify allele-specific expression events. This command takes in a list of germline variants and somatic variants in the VCF format. The module also needs RNAseq alignments produced with a splice-aware aligner in the BAM format and an alignment of the DNA reads in the BAM format. The tool then proceeds to identify polymorphisms that show allele specific expression near the somatic variant sites.

###Usage
`regtools cis-ase identify [options] somatic_variants.vcf polymorphism.vcf dna_alignments.bam rna_alignments.bam ref.fa annotations.gtf`

###Input
| Input                  | Description |
| ------                 | ----------- |
| somatic-variants.vcf   | Somatic variant calls in VCF format. The tool looks for allele specific expression at polymorphic loci near the somatic variants|
| polymorphisms.vcf   | List of polymorphic loci in the VCF format. RNA expression is checked at these sites to identify evidence of allele speciific expression|
| dna-alignments.bam | Aligned DNA reads in the BAM format that has been indexed for example with `samtools index`. We have tested this command with alignments from BWA.|
| dna-alignments.bam | Aligned RNAseq BAM produced with a splice aware aligner, that has been indexed for example with `samtools index`. We have tested this command with alignments from TopHat.|
| ref.fa          | The reference FASTA file. The donor and acceptor sequences used in the "splice-site" column of the annotated junctions are extracted from the FASTA file. |
| annotations.gtf | The GTF file specifies the transcriptome that is used to annotate the junctions and variants. For examples, the Ensembl GTFs for release78 are [here](ftp://ftp.ensembl.org/pub/release-78/gtf/).|

**Note** - Please make sure that the version of the annotation GTF that you use corresponds with the version of the assembly build (ref.fa) and that the co-ordinates in the VCF file are also from the same build.

###Options
| Option  | Description |
| ------  | ----------- |
| -o      | Output file containing the variants that show evidence for allele specific expression. [STDOUT] |
| -d      | INT, Minimum number of reads at the somatic variant and the polymorphic loci to be considered. [10] |
| -w      | INT, Window around the somatic variant to look for transcripts containing polymorphisms exhibiting ASE. [1000] |
| -E      | Flag to look at all neighboring polymorphisms for ASE, not just the exonic polymorphisms. |
| -B      | Flag to use the binomial model and not the default beta model. This feature is under test. |

###Output
The output is in the VCF format and contains a list of polymorphic sites that show evidence for allele specific expression.
TODO - add details about the model parameters here.
