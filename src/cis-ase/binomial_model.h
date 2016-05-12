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

#include "cis_ase_identifier.h"

//Calculate binomial p_het
inline void calculate_binomial_phet(const bcf_call_t& bc, genotype& geno) {
    //RR, RA1, A1A1, RA2, A1A2, A2A2, RA3, A1A3, A2A3, A3A3, RA4 ..
    bool gt_het[15] = {0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0};
    double sum_lik = 0, max_het_lik = 0;
    int n_gt = bc.n_alleles * (bc.n_alleles + 1) / 2;
    for (int i=0; i < n_gt; i++) {
        //convert back from phred
        double lik = pow(10.0, (-1.0 / 10.0 * bc.PL[i]));
        //printf(" PL %d lik %f", bc.PL[i], lik);
        sum_lik += lik;
        //True if GT is het
        if(gt_het[i]) {
            //perhaps switch from max to sum
            if (lik > max_het_lik) {
                max_het_lik = lik;
            }
        }
    }
    geno.p_het = max_het_lik/sum_lik;
}

#endif
