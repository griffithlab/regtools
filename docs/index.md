# RegTools

RegTools is a set of tools that integrate DNA-seq and RNA-seq data to help interpret mutations in a regulatory and splicing context. You can find the source code at our [GitHub repository](https://github.com/griffithlab/regtools) or just use our [Docker image](https://hub.docker.com/r/griffithlab/regtools/) without need for installation.

## Features

- Extract exon-exon junctions from a RNAseq BAM/CRAM file.
- Annotate exon-exon junctions with information from a known transcriptome.
- Annotate variants with splice-region(the definition of this region is configurable) annotations.

## Installation

Clone and install regtools by running:

```sh
git clone https://github.com/griffithlab/regtools
cd regtools/
mkdir build
cd build/
cmake ..
make
```

## Usage

```sh
regtools --help
```

For information about the individual RegTools commands, please see [the Commands page](commands/commands.md)

## Contribute

- [Issue Tracker](https://github.com/griffithlab/regtools/issues)
- [Source Code](https://github.com/griffithlab/regtools)

## Citation

You can find a bioRxiv preprint describing our intial validation of RegTools [here](https://www.biorxiv.org/content/10.1101/436634v2)

## Support

If you have issues using the project, please let us know.
We have a mailing list located at: [regtools@googlegroups.com](mailto:regtools@googlegroups.com) and the forum is [here](https://groups.google.com/forum/#!forum/regtools). Another options is to use our [GitHub issues](https://github.com/griffithlab/regtools/issues) page to contact the project about potential bugs.

## Acknowledgements

RegTools uses several open-source libraries. We would like to thank the
developers of htslib and bedtools. We would also like to thank Travis Abbott for useful comments and code.

## License

The project is licensed under the MIT license.
