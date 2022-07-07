# Overview of `junctions extract` command

The `junctions extract` command can be used to extract exon-exon junctions from an RNAseq BAM file. The output is a BED file in the BED12 format. We have tested this command with alignments from HISAT2, TopHat2, STAR, kallisto, and minimap2 and by comparing the exon-exon junctions with the `junctions.bed` file produced from TopHat.


## Usage

`regtools junctions extract [options] indexed_alignments.bam`

## Input

| Input                  | Description |
| ------                 | ----------- |
| indexed_alignments.bam | Aligned RNAseq BAM/CRAM which has been indexed for example with `samtools index`. We have tested this command with alignments from TopHat.|

## Options

| Option  | Description |
| ------  | ----------- |
| -a      | Minimum anchor length. 8bp by default. Junctions having a minimum overlap of this much on both ends are reported. Note - the required overlap can be observed amongst separate reads, for example one read might have sufficient left overlap and another read might have sufficient right overlap, this is sufficient for the junction to be reported. No mismatches are allowed in the anchor regions.|
| -m      | Minimum intron size. 70bp by default. The intron size is the same as junction.end - junction.start. (Not to be confused with chromStart and chromEnd below, the required blockSizes need to be added/subtracted.)|
| -M      | Maximum intron size. 500,000bp by default. The intron size the same as junction.end - junction.start. (Not to be confused with chromStart and chromEnd below, the required blockSizes need to be added/subtracted.)|
| -o      | File to write output to. STDOUT by default.|
| -r      | Region to extract junctions in. This is specified in the format "chr:start-end". If not specified, junctions are extracted from the entire BAM file.|
| -h      | Display help message for this command.|
| -s      | Strand specificity of RNA library preparation, where the options XS, use XS tags provided by aligner; RF, first-strand; FR, second-strand. This option is required. If your alignments contain XS tags, these will be used in the "unstranded" mode. If you are unsure, we have created this [table](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/) to help.

## Output

The output is in the BED12 format which is described in detail [here.](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) Each line is an exon-exon junction as explained below.

| Column-name       | Description |
| -----------       | ----------- |
| chrom | The name of the chromosome.
| chromStart | The starting position of the junction-anchor. This includes the maximum overhang for the junction on the left. For the exact junction start add blockSizes[0].
| chromEnd | The ending position of the junction-anchor. This includes the maximum overhang for the juncion on the left. For the exact junction end subtract blockSizes[1].
| name | The name of the junctions, the junctions are just numbered JUNC1 to JUNCn.
| score | The number of reads supporting the junction.
| strand | Defines the strand - either '+' or '-'. This is calculated using the XS tag in the BAM file.
| thickStart | Same as chromStart.
| thickEnd | Same as chromEnd.
| itemRgb | RGB value - "255,0,0" by default.
| blockCount | The number of blocks, 2 by default.
| blockSizes | A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
| blockStarts | A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
