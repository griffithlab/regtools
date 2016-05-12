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

#ifndef BETA_MODEL_H
#define BETA_MODEL_H

#include "cis_ase_identifier.h"
#include "Rmath/Rmath.h"
#include <algorithm>

//parameters for no ASE model
int N_alpha = 2000;
int N_beta = 2000;
//parameters for the moderate ASE model
int M_alpha = 36;
int M_beta = 12;
//parameters for the strong ASE model
int S_alpha = 80;
int S_beta = 1;


class BetaModel {
    private:
        int ref_count_;
        int alt_count_;
        //likelihood under no ASE
        float lik_N_;
        //likelihood under moderate ASE
        float lik_M_;
        //likelihood under strong ASE
        float lik_S_;
        //Posterior prob for the S model
        float pp_S_;
        //Posterior prob for the N model
        float pp_N_;
        //Posterior prob for the M model
        float pp_M_;
    public:
        BetaModel() {
            ref_count_ = -1;
            alt_count_ = -1;
            lik_N_ = lik_M_ = lik_S_ = 0;
            pp_N_ = pp_M_ = pp_S_ = 0;
        }
        BetaModel(int ref_count, int alt_count) {
            ref_count_ = ref_count;
            alt_count_ = alt_count;
            lik_N_ = lik_M_ = lik_S_ = 0;
            pp_N_ = pp_M_ = pp_S_ = 0;
        }
        BetaModel(const bcf_call_t& bc) {
            ref_count_ = bc.anno[0] + bc.anno[1];
            alt_count_ = bc.anno[2] + bc.anno[3];
            lik_N_ = lik_M_ = lik_S_ = 0;
            pp_N_ = pp_M_ = pp_S_ = 0;
        }
        //Calculate prob of het
        void calculate_beta_phet(genotype& geno) {
            if(ref_count_ + alt_count_ <= 0) {
                geno.p_het = -1;
                return;
            }
            calc_S_lik();
            calc_M_lik();
            calc_N_lik();
            calculate_posteriors();
            cout << "pp_S_ " << pp_S_;
            cout << " pp_M_ " << pp_M_;
            cout << " pp_N_ " << pp_N_;
            cout << endl;
            if(pp_M_ > 0.5) {
                geno.het_type = "MODASE";
            } else if(pp_S_ > 0.5) {
                geno.het_type = "STRONGASE";
            } else if(pp_N_ > 0.5) {
                geno.het_type = "NOASE";
            }
            //Assign p_het to max of MOD/STRONG
            geno.p_het = std::max(pp_M_, pp_S_);
        }
        //Calculate likelihood under the S model
        void calc_S_lik() {
            float AF = alt_count_/(ref_count_ + alt_count_);
            bool log_density = false;
            lik_S_ = 0.5 * (dbeta(AF, S_alpha, S_beta, log_density) +
                           dbeta(AF, S_beta, S_alpha, log_density));
        }
        //Calculate likelihood under the M model
        void calc_M_lik() {
            float AF = alt_count_/(ref_count_ + alt_count_);
            bool log_density = false;
            lik_M_ = 0.5 * (dbeta(AF, M_alpha, M_beta, log_density) +
                           dbeta(AF, M_beta, M_alpha, log_density));
        }
        //Calculate likelihood under the N model
        void calc_N_lik() {
            float AF = alt_count_/(ref_count_ + alt_count_);
            bool log_density = false;
            lik_N_ = 0.5 * (dbeta(AF, N_alpha, N_beta, log_density) +
                           dbeta(AF, N_beta, N_alpha, log_density));
        }
        //Calculate pp_S, pp_M, pp_N
        void calculate_posteriors() {
            if(lik_S_ + lik_M_ + lik_N_ == 0) {
                throw runtime_error("All likelihoods zero, unable to calculate "
                                    "posterior for beta model\n");
            }
            float total_lik = lik_M_ + lik_N_ + lik_S_;
            pp_M_ = lik_M_/total_lik;
            pp_N_ = lik_N_/total_lik;
            pp_S_ = lik_S_/total_lik;
        }
};

#endif
