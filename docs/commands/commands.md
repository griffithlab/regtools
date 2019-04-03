##Synopsis
To get a list of all available regtools commands `regtools --help`

The main regtools commands are

- [cis-splice-effects](#cis-splice-effects)
- [junctions](#junctions)
- [variants](#variants)

##cis-splice-effects
This set of tools helps identify and work with aberrant splicing events near variants, these could be somatic variants or germline polymorphisms/mutations. These variants are hypothesized to act in cis and affect how the gene is transcribed.

Below are links to detailed explanations of the `cis-splice-effects` sub-commands:

- [identify](cis-splice-effects-identify.md)
- [associate](cis-splice-effects-associate.md)

##cis-ase
This set of tools helps identify and work with allele-specific-expression near variants, these could be somatic variants or germline polymorphisms/mutations. These variants are hypothesized to act in cis and affect how the gene is transcribed.

Below are links to detailed explanations of the `cis-ase` sub-commands:

- [identify](cis-ase-identify.md)

##junctions
The transcriptome structure is often summarized from a RNAseq experiment with a BED file. This BED file contains the exon-exon boundary co-ordinates which are referred to as junctions. For example, TopHat outputs a file called 'junctions.bed' which contains this information. This file is very useful if you are interested in studying which exons/transcripts are expressed, splicing effects etc. On the command line these commands can be accessed using the `regtools junctions` command.

Listed below are links to detailed explanations of the `junctions` sub-commands:

- [extract](junctions-extract.md)
- [annotate](junctions-annotate.md)

##variants
The variants sub-command contains a list of tools that deal with variants that are potentially regulatory in nature. Variants are generally accepted in the standard VCF format unless specified otherwise.

Below are links to detailed explanations of the `variants` sub-commands:

- [annotate](variants-annotate.md)
