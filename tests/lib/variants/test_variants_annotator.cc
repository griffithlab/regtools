/*  test_variants_annotator.cc -- Unit-tests for the VariantsAnnotator class

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
#include "variants_annotator.h"

class VariantsAnnotatorTest : public ::testing::Test {
    public:
        AnnotatedVariant av1, av2, av3;
};

TEST_F(VariantsAnnotatorTest, LessThanOperator) {
    av1 = AnnotatedVariant(std::string("chr1"), 100, 101, "NA");
    av2 = AnnotatedVariant(std::string("chr1"), 100, 101, "NA");
    av3 = AnnotatedVariant(std::string("chr1"), 200, 300, "NA");
    ASSERT_EQ(av1, av2);
    ASSERT_LT(av1, av3);
}
