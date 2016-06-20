/*  test_binomial_model.cc -- Unit tests for the binomial model

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
#include "binomial_model.h"

//Check hom-ref
TEST(BinomialGermlineModel, HomRef) {
    bcf_call_t bc1;
    int n_samples = 1;
    memset(&bc1, 0, sizeof(bcf_call_t));
    bc1.PL = (int32_t *) malloc(15 * n_samples * sizeof(bc1.PL));
    bc1.n_alleles = 2;
    bc1.PL[0] = 0;
    bc1.PL[1] = 255;
    bc1.PL[2] = 255;
    genotype geno1;
    calculate_binomial_germline_phet(bc1, geno1);
    EXPECT_NEAR(0, geno1.p_het, 0.001);
    free(bc1.PL);
}

//Check het
TEST(BinomialGermlineModel, Het) {
    bcf_call_t bc2;
    int n_samples = 1;
    memset(&bc2, 0, sizeof(bcf_call_t));
    bc2.PL = (int32_t *) malloc(15 * n_samples * sizeof(bc2.PL));
    bc2.n_alleles = 2;
    bc2.PL[0] = 255;
    bc2.PL[1] = 0;
    bc2.PL[2] = 255;
    genotype geno1;
    calculate_binomial_germline_phet(bc2, geno1);
    EXPECT_NEAR(1, geno1.p_het, 0.001);
    free(bc2.PL);
}

//Check hom-alt
TEST(BinomialGermlineModel, HomAlt) {
    bcf_call_t bc3;
    int n_samples = 1;
    memset(&bc3, 0, sizeof(bcf_call_t));
    bc3.PL = (int32_t *) malloc(15 * n_samples * sizeof(bc3.PL));
    bc3.n_alleles = 2;
    bc3.PL[0] = 255;
    bc3.PL[1] = 255;
    bc3.PL[2] = 0;
    genotype geno1;
    calculate_binomial_germline_phet(bc3, geno1);
    EXPECT_NEAR(0, geno1.p_het, 0.001);
    free(bc3.PL);
}
