# regtools

regtools is a set of tools that integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context.

##Features

- Extract exon-exon junctions from a RNAseq BAM file.
- Annotate exon-exon junctions with information from a known transcriptome.
- Annotate variants with splice-region(the definition of this region is configurable) annotations.

##Installation

Clone and install regtools by running:
```
    git clone https://github.com/griffithlab/regtools
    cd regtools/
    mkdir build
    cd build/
    cmake ..
    make
```

##Usage

```
    regtools --help
```
For information about the individual regtools commands, please see [the Commands page](commands/commands.md)

##Contribute

- [Issue Tracker](https://github.com/griffithlab/regtools/issues)
- [Source Code](https://github.com/griffithlab/regtools)

##Support

If you have issues using the project, please let us know.
We have a mailing list located at: [regtools@googlegroups.com](mailto:regtools@googlegroups.com) and the forum is [here](https://groups.google.com/forum/#!forum/regtools).
Github issues are another option to contact the project about potential bugs.

##Acknowledgements

regtools uses several open-source libraries. We would like to thank the
developers of htslib and bedtools. We would also like to thank Travis Abbott for
useful comments and code.

##License

The project is licensed under the MIT license.
