[csei]: ../images/csei_examples.png

###Synopsis
The `cis-splice-effects identify` command is used to identify splicing misregulation events. This command takes in a list of variants in the VCF format and RNAseq alignments produced with a splice-aware aligner in the BAM format. The tool then proceeds to identify non-canonical splicing junctions near the variant sites.

###Usage
`regtools cis-splice-effects identify [options] variants.vcf alignments.bam ref.fa annotations.gtf`

###Input
| Input                  | Description |
| ------                 | ----------- |
| variants.vcf | Variant call in VCF format from which to look for cis-splice-effects.|
| alignments.bam | Aligned RNAseq BAM produced with a splice aware aligner, that has been indexed for example with `samtools index`. We have tested this command with alignments from TopHat.|
| ref.fa          | The reference FASTA file. The donor and acceptor sequences used in the "splice-site" column of the annotated junctions are extracted from the FASTA file. |
| annotations.gtf | The GTF file specifies the transcriptome that is used to annotate the junctions and variants. For examples, the Ensembl GTFs for release78 are [here](ftp://ftp.ensembl.org/pub/release-78/gtf/).|

**Note** - Please make sure that the version of the annotation GTF that you use corresponds with the version of the assembly build (ref.fa) and that the co-ordinates in the VCF file are also from the same build.

###Options
| Option  | Description |
| ------  | ----------- |
| -o      | Output file containing the aberrant splice junctions. [STDOUT] |
| -v      | Output file containing variants annotated as splice relevant (VCF format). |
| -w      | Window around the variant file (in basepairs) to identify splicing events in. If specified the tool looks at +/- n b.p around the variant start position. For example -w 500 will look at a 1kb window around the variant. If this option is not specified, the default option is to look at a window that ranges from the start co-ordinate of the previous exon and ends at the end co-ordinate of the next exon i.e by treating the current exon as a cassette exon. |
| -j      | Optional file containing the aberrant junctions in BED12 format. |
| -e      | Maximum distance from the start/end of an exon to annotate a variant as relevant to splicing, the variant is in exonic space, i.e a coding variant. [default = 3] |
| -i      | Maximum distance from the start/end of an exon to annotate a variant as relevant to splicing, the variant is in intronic space. [default = 2] |
| -I      | Annotate variants in intronic space within a transcript (not to be used with -i).
| -E      | Annotate variants in exonic space within a transcript (not to be used with -e).
| -S      | Dont skip single exon transcripts. The default is to skip the single exon transcripts. |

###Output
For an explanation of the annotated junctions that are identified by this command please refer to the output of the `junctions annotate` command [here](junctions-annotate.md#output)
For an explanation of the annotated variants that are identified by this command when using the -v option, please refer to the output of the `variants annotate` command [here](variants-annotate.md#output)

###Examples
![cis-splice-effects identify example][csei]
