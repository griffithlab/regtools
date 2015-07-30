# Examples
[junction_annotation]: images/junction_annotation_examples.png
[anchor_annotation]: images/anchor_examples.png

##Junctions

The transcriptome structure is often summarized from a RNAseq experiment with a BED file.
This BED file contains the exon-exon boundary co-ordinates which are referred to as junctions.
For example, TopHat outputs a file called 'junctions.bed' which contains this information.
This file is very useful if you are interested in studying which exons/transcripts are expressed i.e
splicing effects and so on. The series of sub-sections below explain how regtools works with these junction
files. On the command line these commands can be accessed using the `regtools junctions` command.

###Annotate
The `regtools junctions annotate` command is a tool to annotate the observed junctions with respect to a known
transcript structure. The known transcript structure is in the form of a GTF file obtained from one of the standard
Gene Annotation databases such as Ensembl/RefSeq/UCSC etc. The goal of the annotation step is to help identify novel/unusual junctions.

####Annotating observed junctions with known donor/acceptor/junction information
It is useful to annotate the ends of junction with respect to known acceptors,
donors and junctions in the transcriptome. The known acceptor, donor and junction
information is computed from the GTF file and this information is then used to annotate the observed
junctions.

The junctions are annotated using the following nomenclature(and as shown in the figure below.)

1. DA - This exon-exon junction is present in the transcriptome provided by the user(GTF)
The ends of this junction are hence known donor and known acceptor sites according to the GTF file.

2. NDA - This exon-exon junction is not present(novel) in the transcriptome provided by the user(GTF)
The ends of this junction are known donor and known acceptor sites according to the GTF file.

3. D - This exon-exon junction is not present(novel) in the transcriptome provided by the user(GTF)
The donor of this junction is a known donor but the acceptor is novel.

4. A - This exon-exon junction is not present(novel) in the transcriptome provided by the user(GTF)
The acceptor of this junction is a known acceptor but the donor is novel.

5. N - This exon-exon junction is not present(novel) in the transcriptome provided by the user(GTF)
The ends of this junction are hence not known donor/acceptor sites according to the GTF file.


![Anchor-annotation example][anchor_annotation]

####Annotating a junction


![Junction-annotation example][junction_annotation]

