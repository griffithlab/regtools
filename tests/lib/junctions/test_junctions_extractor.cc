/*  test_junctions_extractor.cc -- Unit-tests for the JunctionsExtractor class

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

#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include "junctions_extractor.h"

class JunctionsExtractTest : public ::testing::Test {
    public:
        JunctionsExtractor jc1;
};

TEST_F(JunctionsExtractTest, ParseInput) {
    int argc = 4;
    char * argv[] = {"extract", "-s", "0", "test_input.bam"};
    int ret = jc1.parse_options(argc, argv);
    string expected_bam("test_input.bam");
    ASSERT_EQ(expected_bam, jc1.get_bam());
    ASSERT_EQ(0, ret);
}

TEST_F(JunctionsExtractTest, ParseNoInput) {
    int argc = 1;
    char * argv[] = {"extract"};
    ASSERT_THROW(jc1.parse_options(argc, argv), std::runtime_error);
}

TEST_F(JunctionsExtractTest, ParseIncorrectOption) {
    int argc = 2;
    char * argv[] = {"extract", "-k", "24", "test_input.bam"};
    ASSERT_THROW(jc1.parse_options(argc, argv), std::runtime_error);
}

// TEST_F(JunctionsExtractTest, Usage) {
//     ostringstream out, out2;
//     out << "Usage:" 
//         << "\t\t" << "regtools junctions extract [options] indexed_alignments.bam" << endl;
//     out << "Options:" << endl;
//     out << "\t\t" << "-a INT\tMinimum anchor length. Junctions which satisfy a minimum \n"
//         << "\t\t\t " << "anchor length on both sides are reported. [8]" << endl;
//     out << "\t\t" << "-i INT\tMinimum intron length. [70]" << endl;
//     out << "\t\t" << "-I INT\tMaximum intron length. [500000]" << endl;
//     out << "\t\t" << "-o FILE\tThe file to write output to. [STDOUT]" << endl;
//     out << "\t\t" << "-r STR\tThe region to identify junctions \n"
//         << "\t\t\t " << "in \"chr:start-end\" format. Entire BAM by default." << endl;
//     out << "\t\t" << "-s INT\tStrand specificity of RNA library preparation \n"
//         << "\t\t\t " << "(0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR). [1]" << endl;
//     out << endl;
//     jc1.usage(out2);
//     ASSERT_EQ(out.str(), out2.str()) << "Error parsing as expected";
// }

TEST_F(JunctionsExtractTest, JunctionName) {
    string j1_name = jc1.get_new_junction_name();
    ASSERT_EQ(j1_name, string("JUNC00000001"));
    Junction j1("chr1", 10000, 10200,
            9500, 10700, "+");
    jc1.add_junction(j1);
    string j2_name = jc1.get_new_junction_name();
    ASSERT_EQ(string("JUNC00000002"), j2_name);
}

TEST_F(JunctionsExtractTest, PrintJunction) {
    stringstream ss1, expected;
    Junction j1("chr1", 10000, 10200,
            9500, 10700, "+");

    expected << "chr1" <<
        "\t" << 9500 << "\t" << 10700 <<
        "\t" << "NA" << "\t" << 0 << "\t" << "+" <<
        "\t" << 9500 << "\t" << 10700 <<
        "\t" << "255,0,0" << "\t" << 2 <<
        "\t" << 10000 - 9500 << "," << 10700 - 10200 <<
        "\t" << "0," << 10200 - 9500  << endl;

    j1.print(ss1);
    ASSERT_EQ(expected.str(), ss1.str());
}

TEST_F(JunctionsExtractTest, AddJunction) {
    stringstream ss1, expected;
    //Add one junction with differing thick start/ends
    jc1.add_junction(Junction("chr1", 10000, 10200,
            9900, 10300, "+"));
    jc1.add_junction(Junction("chr1", 10000, 10200,
            9500, 10200, "+"));
    jc1.add_junction(Junction("chr1", 10000, 10200,
            9950, 10700, "+"));
    //Add second junction
    jc1.add_junction(Junction("chr1", 8000, 8500,
            7000, 10000, "+"));
    //Add second junction with different strand
    jc1.add_junction(Junction("chr1", 8000, 8500,
            7000, 10000, "-"));
    expected << "chr1" <<
        "\t" << 7000 << "\t" << 10000 <<
        "\t" << "JUNC00000002" << "\t" << 1 << "\t" << "+" <<
        "\t" << 7000 << "\t" << 10000 <<
        "\t" << "255,0,0" << "\t" << 2 <<
        "\t" << 8000 - 7000 << "," << 10000 - 8500 <<
        "\t" << "0," << 8500 - 7000  << endl;
    expected << "chr1" <<
        "\t" << 7000 << "\t" << 10000 <<
        "\t" << "JUNC00000003" << "\t" << 1 << "\t" << "-" <<
        "\t" << 7000 << "\t" << 10000 <<
        "\t" << "255,0,0" << "\t" << 2 <<
        "\t" << 8000 - 7000 << "," << 10000 - 8500 <<
        "\t" << "0," << 8500 - 7000  << endl;
    expected << "chr1" <<
        "\t" << 9500 << "\t" << 10700 <<
        "\t" << "JUNC00000001" << "\t" << 3 << "\t" << "+" <<
        "\t" << 9500 << "\t" << 10700 <<
        "\t" << "255,0,0" << "\t" << 2 <<
        "\t" << 10000 - 9500 << "," << 10700 - 10200 <<
        "\t" << "0," << 10200 - 9500  << endl;
    //This will sort and then print
    jc1.print_all_junctions(ss1);
    ASSERT_EQ(expected.str(), ss1.str());
}
