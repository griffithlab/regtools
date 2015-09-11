#Regtools commands

This page summarizes the regtools commands.

To get a list of all available regtools, you can run
```
    regtools --help
```

##Junctions

The transcriptome structure is often summarized from a RNAseq experiment with a BED file. This BED file contains the exon-exon boundary co-ordinates which are referred to as junctions. For example, TopHat outputs a file called 'junctions.bed' which contains this information. This file is very useful if you are interested in studying which exons/transcripts are expressed i.e splicing effects and so on. The series of sub-sections below explain how regtools works with these junction files. On the command line these commands can be accessed using the `regtools junctions` command.

[junctions-annotate](commands/junctions-annotate.md)
