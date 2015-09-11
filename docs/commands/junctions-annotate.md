# Examples
[junction_annotation]: images/junction_annotation_examples.png
[anchor_annotation]: images/anchor_examples.png

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

####Annotating a junction with number of donors/acceptors/exons skipped
Exon skipping is a form of RNA splicing that can be identified using RNAseq data. It is hence useful
to compute for every observed putative exon-exon junction, the number of exons skipped, the number of
known donor sites skipped and the number of known acceptor sites skipped. The known exons, donors and
acceptors are calculated from the user supplied GTF file.

In the example shown below, the observed junction has skipped 2 known donor sites and 2 known acceptor sites.
For the number of exons skipped we consider two situations. The second exon in all the three transcripts overlap,
if these overlapping exons are merged together, the number of exons skipped is just 1. If these three exons are
considered to be different the number of exons skipped is 3. We try and provide both these annotations.

![Junction-annotation example][junction_annotation]

If any of the examples are not clear or if you would like more information please feel free to open an issue on GitHub [here](https://github.com/griffithlab/regtools)
or post on the discussion page [here.](https://groups.google.com/d/forum/regtools)
