/*  test_beta_model.cc -- Unit tests for the Beta model

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
#include "cis_ase_identifier.h"
#include "beta_model.h"

class BetaModelTest : public ::testing::Test {
    public:
        BetaModel bm1, bm2, bm3, bm4;
        void SetUp() {
            bm2 = BetaModel(50, 50);
            bm3 = BetaModel(75, 25);
            bm4 = BetaModel(100, 3);
        }
};

//Uninitialized case
TEST_F(BetaModelTest, Construct1) {
    genotype geno1;
    bm1.calculate_beta_phet(geno1);
    EXPECT_EQ(-1, geno1.p_het);
}

//CHECK NO ASE
TEST_F(BetaModelTest, NoAse) {
    genotype geno1;
    bm2.calculate_beta_phet(geno1);
    EXPECT_NEAR(1, geno1.p_het, 0.001);
    EXPECT_EQ("NOASE", geno1.het_type);
}

//CHECK MOD ASE
TEST_F(BetaModelTest, ModAse) {
    genotype geno1;
    bm3.calculate_beta_phet(geno1);
    EXPECT_NEAR(0, geno1.p_het, 0.001);
    EXPECT_EQ("MODASE", geno1.het_type);
}

//CHECK STRONG ASE
TEST_F(BetaModelTest, StrongAse) {
    genotype geno1;
    bm4.calculate_beta_phet(geno1);
    EXPECT_NEAR(0, geno1.p_het, 0.001);
    EXPECT_EQ("STRONGASE", geno1.het_type);
}
