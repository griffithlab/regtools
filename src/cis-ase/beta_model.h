/*  beta_model.h -- implement the model in http://dx.doi.org/10.1101/007211

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
#include "Rmath/Rmath.h"

//parameters for no ASE model
N_alpha = 2000;
N_beta = 2000;
//parameters for the moderate ASE model
M_alpha = 36;
M_beta = 12;
//parameters for the strong ASE model
S_alpha = 80;
S_beta = 1;


class BetaModel {
    private:
        int ref_count_;
        int alt_count_;
        //likelihood under no ASE
        float lik_N;
        //likelihood under moderate ASE
        float lik_M;
        //likelihood under strong ASE
        float lik_S;
    public:
        BetaModel() {
            ref_count_ = -1;
            alt_count_ = -1;
        }
        BetaModel(int ref_count, int alt_count) {
            ref_count_ = ref_count;
            alt_count_ = alt_count;
        }
        BetaModel(const bcf_call_t& bc) {
            ref_count_ = bc.anno[0] + bc.anno[1];
            alt_count_ = bc.anno[2] + bc.anno[3];
        }
        calculate_beta_phet(genotype geno) {
            if(ref_count_ + alt_count_ == 0) {
                geno.p_het = -1;
            }
            geno.p_het = max_het_lik/sum_lik;
        }
        calc_S_lik() {
            float AF = alt_count_/(ref_count_ + alt_count_);
            lik_S = 0.5 * (pbeta(
        }
}

#endif
