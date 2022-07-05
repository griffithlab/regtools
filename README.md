![Build Status](https://github.com/griffithlab/regtools/actions/workflows/cmake.yml/badge.svg?branch=master)
[![Documentation Status](https://readthedocs.org/projects/regtools/badge/?version=latest)](https://readthedocs.org/projects/regtools/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/griffithlab/regtools/badge.svg?branch=master&service=github)](https://coveralls.io/github/griffithlab/regtools?branch=master)

# RegTools

Tools that integrate DNA-seq and RNA-seq data to help interpret mutations
in a regulatory and splicing context.

## Features

- Identify evidence for aberrant splicing in RNA reads near a list of variants.
- Extract exon-exon junctions from a RNAseq BAM file.
- Annotate exon-exon junctions with information from a known transcriptome.
- Annotate variants with splice-region(the definition of this region is configurable) annotations.

## Hardware requirements
RegTools  requires only a standard computer with enough RAM to support the in-memory operations.

## Software requirements
OS Requirements
This package is supported for macOS and Linux. The package has been tested on the following systems:

macOS: macOS 10.12 (Sierra), macOS 10.13 (High Sierra), macOS 10.14 (Mojave), macOS 10.15 (Catalina), macOS 11 (Big Sur), macOS 12 (Monterey)

Linux: Ubuntu 16.04, Ubuntu 18.04, Ubuntu 20.04

## Installation

Clone and install regtools by running the following:
```
    git clone https://github.com/griffithlab/regtools
    cd regtools/
    mkdir build
    cd build/
    cmake ..
    make
```

For convienience we also maintain a docker image available at [https://hub.docker.com/r/griffithlab/regtools/](https://hub.docker.com/r/griffithlab/regtools/)

## Usage:

```
    regtools --help
```

## Contribute

- Issue Tracker: github.com/griffithlab/regtools/issues
- Source Code: github.com/griffithlab/regtools

## Support

If you have issues using the project, please let us know.
We have a mailing list located at: regtools@googlegroups.com and the
forum is here - https://groups.google.com/forum/#!forum/regtools.
Github issues are another option to contact the project about
potential bugs.

## Documentation

The documentation for the project is hosted on
[Read the Docs.](https://regtools.readthedocs.org/en/latest/)

If you would like to build the documentation locally, please install
[mkdocs](http://www.mkdocs.org/), `pip install mkdocs --user` should
work on most machines. Then run `mkdocs serve` from within the `regtools`
base directory.


## Acknowledgements

Regtools uses several open-source libraries. We would like to thank the
developers of htslib and bedtools. We would also like to thank Travis Abbott for
useful comments and code.

## License

The project is licensed under the [MIT license](https://opensource.org/licenses/MIT).

## Stable release with DOI

[![DOI](https://zenodo.org/badge/35841695.svg)](https://zenodo.org/badge/latestdoi/35841695)


