/*  binomial_model.h -- Use the binomial likelihoods to calculate prob_het

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

#ifndef BINOMIAL_MODEL_H
#define BINOMIAL_MODEL_H

#include <cmath>
#include "cis_ase_identifier.h"
#include "Rmath/Rmath.h"

//Calculate binomial p_het for germline variants
inline void calculate_binomial_germline_phet(const bcf_call_t& bc, genotype& geno) {
    uint32_t ref_count_ = bc.anno[0] + bc.anno[1];
    uint32_t alt_count_ = bc.anno[2] + bc.anno[3];
    // Assume a uniform beta prior of beta(1, 9)
    double alpha = 1 + alt_count_;
    double beta = 1 + ref_count_;
    //Last two arguments of pbeta specifies if lower.tail, returned prob is log
    double p_het = pbeta(0.6, alpha, beta, true, false) - pbeta(0.4, alpha, beta, true, false);
    geno.p_het = (double) p_het;
    //0.0 - 0.1
    double p_homref = pbeta(0.1, alpha, beta, true, false);
    //0.9 - 1.0, note: upper tail used here for pbeta
    double p_homalt = pbeta(0.9, alpha, beta, false, false);
    //geno.p_het = (double) p_het / (double) (p_het + p_homref + p_homalt);
    cerr << "inside beta " << ref_count_ << "\t" << alt_count_ << "\t" << p_het << "\t" << p_homref <<
            "\t" << p_homalt << "\t" << geno.p_het << endl;
}

//Calculate binomial p_het for somatic variants - be more permissive with AF
inline void calculate_binomial_somatic_phet(const bcf_call_t& bc, genotype& geno) {
    uint32_t ref_count_ = bc.anno[0] + bc.anno[1];
    uint32_t alt_count_ = bc.anno[2] + bc.anno[3];
    // Assume a uniform beta prior of beta(1, 9)
    double alpha = 1 + alt_count_;
    double beta = 1 + ref_count_;
    //Last two arguments of pbeta specifies if lower.tail, returned prob is log
    double p_het = pbeta(0.8, alpha, beta, true, false) - pbeta(0.2, alpha, beta, true, false);
    geno.p_het = (double) p_het;
    //0.0 - 0.1
    double p_homref = pbeta(0.25, alpha, beta, true, false);
    //0.9 - 1.0, note: upper tail used here for pbeta
    double p_homalt = pbeta(0.75, alpha, beta, false, false);
    //geno.p_het = (double) p_het / (double) (p_het + p_homref + p_homalt);
    cerr << "inside beta " << ref_count_ << "\t" << alt_count_ << "\t" << p_het << "\t" << p_homref <<
            "\t" << p_homalt << "\t" << geno.p_het << endl;
}

#endif
