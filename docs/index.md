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

If one wishes to test their installation, we include test data under `test_data`. 

Here's an example command using that data along with the example output.

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
