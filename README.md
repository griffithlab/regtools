[![Build Status](https://travis-ci.org/griffithlab/regtools.svg?branch=master)](https://travis-ci.org/griffithlab/regtools)
[![Documentation Status](https://readthedocs.org/projects/regtools/badge/?version=latest)](https://readthedocs.org/projects/regtools/?badge=latest)
[![Coverage Status](https://coveralls.io/repos/griffithlab/regtools/badge.svg?branch=master&service=github)](https://coveralls.io/github/griffithlab/regtools?branch=master)

# regtools

Tools that integrate DNA-seq and RNA-seq data to help interpret mutations
in a regulatory and splicing context.

## Features

- Identify evidence for aberrant splicing in RNA reads near a list of variants.
- Extract exon-exon junctions from a RNAseq BAM file.
- Annotate exon-exon junctions with information from a known transcriptome.
- Annotate variants with splice-region(the definition of this region is configurable) annotations.

## Installation

Clone and install regtools by running:
```
    git clone https://github.com/griffithlab/regtools
    cd regtools/
    mkdir build
    cd build/
    cmake ..
    make
```

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
