/*  test_junctions_creator.cc -- Unit-tests for the JunctionsCreator class

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
#include "junctions_creator.h"

class JunctionsCreateTest : public ::testing::Test {
    public:
        JunctionsCreator jc1;
};

TEST_F(JunctionsCreateTest, ParseInput) {
    int argc = 2;
    char * argv[] = {"create", "test_input.bam"};
    int ret = jc1.parse_options(argc, argv);
    string expected_bam("test_input.bam");
    ASSERT_EQ(expected_bam, jc1.get_bam());
    ASSERT_EQ(0, ret);
}

TEST_F(JunctionsCreateTest, ParseNoInput) {
    int argc = 1;
    char * argv[] = {"create"};
    ASSERT_THROW(jc1.parse_options(argc, argv), std::runtime_error);
}

TEST_F(JunctionsCreateTest, ParseIncorrectOption) {
    int argc = 2;
    char * argv[] = {"create", "-k", "24", "test_input.bam"};
    ASSERT_THROW(jc1.parse_options(argc, argv), std::runtime_error);
}

TEST_F(JunctionsCreateTest, Usage) {
    ostringstream out, out2;
    out << "\nUsage:\t\t" << "regtools junctions create [options] indexed_alignments.bam";
    out << "\nOptions:";
    out << "\t" << "-a INT\tMinimum anchor length. Junctions which satisfy a minimum "
                     "anchor length on both sides are reported. [8]";
    out << "\n\t\t" << "-i INT\tMinimum intron length. [70]";
    out << "\n\t\t" << "-I INT\tMaximum intron length. [500000]";
    out << "\n\t\t" << "-o FILE\tThe file to write output to. [STDOUT]";
    out << "\n\t\t" << "-r STR\tThe region to identify junctions "
                     "in \"chr:start-end\" format. Entire BAM by default.";
    out << "\n";
    jc1.usage(out2);
    ASSERT_EQ(out.str(), out2.str()) << "Error parsing as expected";
}
