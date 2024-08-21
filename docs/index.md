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

If one wishes to test their installation, we include test data under `tests/integration-test/data`.

Here's an example command using that data along with the example output.

```sh
regtools cis-splice-effects identify -s RF tests/integration-test/data/vcf/test1.vcf tests/integration-test/data/bam/test_hcc1395.2.bam tests/integration-test/data/fa/test_chr22.fa tests/integration-test/data/gtf/test_ensemble_chr22.2.gtf

chrom	start	end	name	score	strand	splice_site	acceptors_skipped	exons_skipped	donors_skipped	anchor	known_donor	known_acceptor	known_junction	gene_names	gene_ids	transcripts	variant_info
22	93668	97252	JUNC00000001	5	+	GT-AG	1	1	1	NDA	1	1	0	EP300	ENSG00000100393	ENST00000263253	22:94626-94627
```

For information about the individual RegTools commands, please see [the Commands page](commands/commands.md)

## Contribute

- [Issue Tracker](https://github.com/griffithlab/regtools/issues)
- [Source Code](https://github.com/griffithlab/regtools)

## Citation

Cotto KC*, Feng YY*, Ramu A, Richters M, Freshour SL, Skidmore ZL, Xia H, McMichael JF, Kunisaki J, Campbell KM, Chen TH, Rozycki EB, Adkins D, Devarakonda S, Sankararaman S, Lin Y, Chapman WC, Maher CA, Arora V, Dunn GP, Uppaluri R, Govindan R, Griffith OL, Griffith M. **Integrated analysis of genomic and transcriptomic data for the discovery of splice-associated variants in cancer.** Nature Communications. 2023 Mar.
pmid: [36949070](https://pubmed.ncbi.nlm.nih.gov/36949070/); doi: [10.1038/s41467-023-37266-6](https://www.nature.com/articles/s41467-023-37266-6)

## Support

If you have issues using the project, please let us know.
We have a mailing list located at: [regtools@googlegroups.com](mailto:regtools@googlegroups.com) and the forum is [here](https://groups.google.com/forum/#!forum/regtools). Another options is to use our [GitHub issues](https://github.com/griffithlab/regtools/issues) page to contact the project about potential bugs.

## Acknowledgements

RegTools uses several open-source libraries. We would like to thank the developers of htslib and bedtools. We would also like to thank Travis Abbott for useful comments and code.

## License

The project is licensed under the MIT license.
