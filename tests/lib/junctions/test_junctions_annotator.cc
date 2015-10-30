/*  test_junctions_annotator.cc -- Unit-tests for the JunctionsAnnotator class

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
#include "junctions_annotator.h"

class JunctionsAnnotatorTest : public ::testing::Test {
    public:
        JunctionsAnnotator ja1;
};

TEST_F(JunctionsAnnotatorTest, ParseInput) {
    int argc = 4;
    char * argv[] = {"annotate",
                     "test.bed",
                     "test.fa",
                     "test.gtf"};
    int ret = ja1.parse_options(argc, argv);
    string expected_gtf("test.gtf");
    ASSERT_EQ(expected_gtf, ja1.gtf_file());
    ASSERT_EQ(0, ret);
}
