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

Installation should take 1-5 minutes.

For convienience we also maintain a docker image available at [https://hub.docker.com/r/griffithlab/regtools/](https://hub.docker.com/r/griffithlab/regtools/)

## Usage:

```
    regtools --help
```

If one wishes to test their installation, we include test data under `test_data`. 

Here's an example command using that data along with the example output. This should run in under a minute.

```sh
regtools cis-splice-effects identify -s RF -e 10 -i 10 test_data/HCC1395_chr22.vcf.gz test_data/HCC1395_tumor.bam test_data/chr22_with_ERCC92.fa test_data/chr22_with_ERCC92.gtf

Variant 22	42129188	42129189	-1
Variant region is 22:42128784-42130813

chrom	start	end	name	score	strand	splice_site	acceptors_skipped	exons_skipped	donors_skipped	anchor	known_donor	known_acceptor	known_junction	gene_names	gene_ids	transcripts	variant_info
position = 22:42125408-42125409
position = 22:42130565-42130566
22	42125407	42130567	JUNC00000001	4	+	GT-AG	0	0	0	D	1	0	0	NDUFA6-AS1	ENSG00000237037	ENST00000439129	22:42129188-42129189
position = 22:42128881-42128882
position = 22:42129670-42129671
22	42128880	42129672	JUNC00000002	3	+	GT-AG	0	0	0	N	0	0	0	NA	NA	NA	22:42129188-42129189
position = 22:42128944-42128945
position = 22:42129031-42129032
22	42128943	42129033	JUNC00000003	4	-	GT-GG	1	0	0	D	1	0	0	CYP2D6	ENSG00000100197	ENST00000360608,ENST00000389970,ENST00000488442	22:42129188-42129189
position = 22:42129783-42129784
position = 22:42143453-42143454
22	42129782	42143455	JUNC00000004	2	+	GC-AG	9	8	9	N	0	0	0	NA	NA	NA	22:42129188-42129189
position = 22:42130224-42130225
position = 22:42130565-42130566
22	42130223	42130567	JUNC00000005	2	+	GT-AG	0	0	0	N	0	0	0	NA	NA	NA	22:42129188-42129189
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


