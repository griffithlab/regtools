/*  test_cis_splice_effects_identifier.cc -- Unit tests for the CisSpliceEffectsIdentifier class

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
#include <stdexcept>
#include "cis_splice_effects_identifier.h"
#include "common.h"

class CisSpliceEffectsIdentifierTest : public ::testing::Test {
    public:
        CisSpliceEffectsIdentifier csei1;
};

//Check for missing input
TEST_F(CisSpliceEffectsIdentifierTest, ParseInput) {
    int argc = 4;
    char * argv[] = {"identify",
                     "test.vcf",
                     "test.fa",
                     "test.gtf"};
    EXPECT_THROW(csei1.parse_options(argc, argv), std::runtime_error);
}

//Check for unwanted input
TEST_F(CisSpliceEffectsIdentifierTest, ParseInput2) {
    int argc = 5;
    char * argv[] = {"identify",
                     "test.vcf",
                     "test.bam",
                     "test.fa",
                     "test.gtf",
                     "random_input"};
    EXPECT_THROW(csei1.parse_options(argc, argv), std::runtime_error);
}

//Check the -h option
TEST_F(CisSpliceEffectsIdentifierTest, ParseInput3) {
    int argc = 5;
    char * argv[] = {"identify",
                     "-h"};
    EXPECT_THROW(csei1.parse_options(argc, argv), common::cmdline_help_exception);
}

//Test if constructor works as expected.
TEST_F(CisSpliceEffectsIdentifierTest, ConstructorTest) {
    EXPECT_EQ(csei1.vcf(), "NA");
    EXPECT_EQ(csei1.output_file(), "NA");
    EXPECT_EQ(csei1.annotated_variant_file(), "NA");
    EXPECT_EQ(csei1.window_size(), 0);
}
