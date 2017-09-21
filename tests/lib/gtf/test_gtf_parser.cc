/*  test_gtf_parser.cc -- Unit-tests for the GtfParser class

    Copyright (c) 2015, The Griffith Lab

    Author: Avinash Ramu <aramu@genome.wustl.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <array>
#include <string>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include "gtf_parser.h"

class GtfParserTest : public ::testing::Test {
    public:
        GtfParser gp1;
};

//Check if the gtf filename is properly parsed
TEST_F(GtfParserTest, FileTest) {
    string gtf_ = "test.gtf";
    gp1.set_gtffile(gtf_);
    EXPECT_EQ(gtf_, gp1.gtffile());
}

//Parse an exon line into a gtf struct
TEST_F(GtfParserTest, ParseExonLineTest) {
    string line1 = "22\tprotein_coding\texon\t12791\t14103\t.\t+\t."
                   "\tccds_id \"CCDS14010\"; exon_id \"ENSE00001343011"
                   "\"; exon_number \"1\"; gene_biotype \"protein_coding\";"
                   " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                   "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                   "\"CCDS\"; transcript_id \"ENST00000263253\"; transcript_name "
                   "\"EP300-001\"; transcript_source \"ensembl_havana\"; "
                   "tss_id \"TSS138009\"";
    Gtf expected_gtf1;
    expected_gtf1.seqname = "22";
    expected_gtf1.source = "protein_coding";
    expected_gtf1.feature = "exon";
    expected_gtf1.start = 12791;
    expected_gtf1.end = 14103;
    expected_gtf1.score = ".";
    expected_gtf1.strand = "+";
    expected_gtf1.frame = '.';
    expected_gtf1.attributes = "ccds_id \"CCDS14010\"; exon_id \"ENSE00001343011"
                               "\"; exon_number \"1\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    expected_gtf1.is_exon = true;
    EXPECT_EQ(expected_gtf1, gp1.parse_exon_line(line1));
}

//Parse the required field from attributes column
TEST_F(GtfParserTest, ParseAttributeTest) {
    vector<string> attributes1;
    attributes1.push_back("ccds_id \"CCDS14010\"");
    attributes1.push_back("gene_source \"ensembl_havana\"");
    attributes1.push_back("tss_id \"TSS138009\"");
    EXPECT_EQ("TSS138009", gp1.parse_attribute(attributes1, "tss_id"));
    EXPECT_EQ("CCDS14010", gp1.parse_attribute(attributes1, "ccds_id"));
    EXPECT_EQ("ensembl_havana", gp1.parse_attribute(attributes1, "gene_source"));
    EXPECT_EQ("NA", gp1.parse_attribute(attributes1, "fake_attr"));
}

//Add an exon to the transcript map
TEST_F(GtfParserTest, AddExonToTranscriptTest) {
    Gtf test_gtf1;
    test_gtf1.seqname = "22";
    test_gtf1.source = "protein_coding";
    test_gtf1.feature = "exon";
    test_gtf1.start = 12791;
    test_gtf1.end = 14103;
    test_gtf1.score = ".";
    test_gtf1.strand = "+";
    test_gtf1.frame = '.';
    test_gtf1.attributes = "ccds_id \"CCDS14010\"; exon_id \"ENSE00001343011"
                               "\"; exon_number \"1\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    test_gtf1.is_exon = true;
    gp1.add_exon_to_transcript_map(test_gtf1);
    EXPECT_EQ("EP300",
              gp1.get_gene_from_transcript("ENST00000263253")[0]);
    EXPECT_EQ("ENSG00000100393",
              gp1.get_gene_from_transcript("ENST00000263253")[1]);
    EXPECT_EQ("NA, NA",
              gp1.get_gene_from_transcript("ENSTfake")[0]);
    EXPECT_EQ("NA, NA",
              gp1.get_gene_from_transcript("ENSTfake")[1]);
    gp1.annotate_transcript_with_bins();
    EXPECT_EQ(37359u,
              gp1.bin_from_transcript("ENST00000263253"));
    std::vector<string> expected_transcript;
    expected_transcript.push_back("ENST00000263253");
    EXPECT_EQ(expected_transcript, gp1.transcripts_from_bin("22", 37359));
}

//Test sorting of exons within a positive-strand transcript
TEST_F(GtfParserTest, SortExonTranscriptPsTest) {
    Gtf gtf1;
    gtf1.seqname = "22";
    gtf1.source = "protein_coding";
    gtf1.feature = "exon";
    gtf1.start = 10100;
    gtf1.end = 10200;
    gtf1.score = ".";
    gtf1.strand = "+";
    gtf1.frame = '.';
    gtf1.attributes = "ccds_id \"CCDS14010\"; exon_id \"ENSE00001343011"
                               "\"; exon_number \"3\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    gtf1.is_exon = true;
    Gtf gtf2;
    gtf2.seqname = "22";
    gtf2.source = "protein_coding";
    gtf2.feature = "exon";
    gtf2.start = 9900;
    gtf2.end = 10000;
    gtf2.score = ".";
    gtf2.strand = "+";
    gtf2.frame = '.';
    gtf2.attributes = "ccds_id \"CCDS14011\"; exon_id \"ENSE00001343012"
                               "\"; exon_number \"2\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    gtf2.is_exon = true;
    Gtf gtf3;
    gtf3.seqname = "22";
    gtf3.source = "protein_coding";
    gtf3.feature = "exon";
    gtf3.start = 9700;
    gtf3.end = 9800;
    gtf3.score = ".";
    gtf3.strand = "+";
    gtf3.frame = '.';
    gtf3.attributes = "ccds_id \"CCDS14012\"; exon_id \"ENSE00001343013"
                               "\"; exon_number \"1\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    gtf3.is_exon = true;
    vector<BED> expected_exons1, expected_exons2;
    BED exon1 = BED(gtf1.seqname, gtf1.start,
                    gtf1.end, gtf1.feature,
                    gtf1.score, gtf1.strand);
    BED exon2 = BED(gtf2.seqname, gtf2.start,
                    gtf2.end, gtf2.feature,
                    gtf2.score, gtf2.strand);
    BED exon3 = BED(gtf3.seqname, gtf3.start,
                    gtf3.end, gtf3.feature,
                    gtf3.score, gtf3.strand);
    //Not sorted
    expected_exons1.push_back(exon1);
    expected_exons1.push_back(exon2);
    expected_exons1.push_back(exon3);
    //Sorted
    expected_exons2.push_back(exon3);
    expected_exons2.push_back(exon2);
    expected_exons2.push_back(exon1);
    //Add to parser
    gp1.add_exon_to_transcript_map(gtf1);
    gp1.add_exon_to_transcript_map(gtf2);
    gp1.add_exon_to_transcript_map(gtf3);
    EXPECT_EQ(expected_exons1, gp1.get_exons_from_transcript("ENST00000263253"));
    gp1.sort_exons_within_transcripts();
    EXPECT_EQ(expected_exons2, gp1.get_exons_from_transcript("ENST00000263253"));
}

//Test sorting of exons within a negative-strand transcript
TEST_F(GtfParserTest, SortExonTranscriptNsTest) {
    Gtf gtf1;
    gtf1.seqname = "22";
    gtf1.source = "protein_coding";
    gtf1.feature = "exon";
    gtf1.start = 10100;
    gtf1.end = 10200;
    gtf1.score = ".";
    gtf1.strand = "-";
    gtf1.frame = '.';
    gtf1.attributes = "ccds_id \"CCDS14010\"; exon_id \"ENSE00001343011"
                               "\"; exon_number \"3\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    gtf1.is_exon = true;
    Gtf gtf2;
    gtf2.seqname = "22";
    gtf2.source = "protein_coding";
    gtf2.feature = "exon";
    gtf2.start = 9900;
    gtf2.end = 10000;
    gtf2.score = ".";
    gtf2.strand = "-";
    gtf2.frame = '.';
    gtf2.attributes = "ccds_id \"CCDS14011\"; exon_id \"ENSE00001343012"
                               "\"; exon_number \"2\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    gtf2.is_exon = true;
    Gtf gtf3;
    gtf3.seqname = "22";
    gtf3.source = "protein_coding";
    gtf3.feature = "exon";
    gtf3.start = 9700;
    gtf3.end = 9800;
    gtf3.score = ".";
    gtf3.strand = "-";
    gtf3.frame = '.';
    gtf3.attributes = "ccds_id \"CCDS14012\"; exon_id \"ENSE00001343013"
                               "\"; exon_number \"1\"; gene_biotype \"protein_coding\";"
                               " gene_id \"ENSG00000100393\"; gene_name \"EP300\"; "
                               "gene_source \"ensembl_havana\"; p_id \"P5137\"; tag "
                               "\"CCDS\"; transcript_id \"ENST00000263253\"; "
                               "transcript_name \"EP300-001\"; transcript_source "
                               "\"ensembl_havana\"; tss_id \"TSS138009\"";
    gtf3.is_exon = true;
    vector<BED> expected_exons1, expected_exons2;
    BED exon1 = BED(gtf1.seqname, gtf1.start,
                    gtf1.end, gtf1.feature,
                    gtf1.score, gtf1.strand);
    BED exon2 = BED(gtf2.seqname, gtf2.start,
                    gtf2.end, gtf2.feature,
                    gtf2.score, gtf2.strand);
    BED exon3 = BED(gtf3.seqname, gtf3.start,
                    gtf3.end, gtf3.feature,
                    gtf3.score, gtf3.strand);
    //Not sorted
    expected_exons1.push_back(exon2);
    expected_exons1.push_back(exon3);
    expected_exons1.push_back(exon1);
    //Sorted
    expected_exons2.push_back(exon1);
    expected_exons2.push_back(exon2);
    expected_exons2.push_back(exon3);
    //Add to parser
    gp1.add_exon_to_transcript_map(gtf2);
    gp1.add_exon_to_transcript_map(gtf3);
    gp1.add_exon_to_transcript_map(gtf1);
    EXPECT_EQ(expected_exons1, gp1.get_exons_from_transcript("ENST00000263253"));
    gp1.sort_exons_within_transcripts();
    EXPECT_EQ(expected_exons2, gp1.get_exons_from_transcript("ENST00000263253"));
}
