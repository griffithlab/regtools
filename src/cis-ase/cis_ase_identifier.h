/*  cis_ase_identifier.h -- 'cis-ase identify'

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

#ifndef CIS_ASE_IDENTIFIER_
#define CIS_ASE_IDENTIFIER_

#include <iostream>
#include <fstream>
#include <map>
#include "htslib/sam.h"
#include "bam2bcf.h"
#include "bam_plcmd.h"
#include "hts.h"
#include "vcf.h"
#include "sample.h"

using namespace std;

extern const double MIN_HET_PROB;

extern "C" {
    void *bed_read(const char *fn);
    void bed_destroy(void *_h);
    int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len);
    void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
            int n, char *const*fn, int *n_plp,
            const bam_pileup1_t **plp, int ignore_rg);
    int bed_overlap(const void *_h, const char *chr, int beg, int end);
    int mpileup_with_likelihoods(mplp_conf_t *conf, int n, char **fn, mplp_aux_t **data, bcf_callaux_t *bca, bcf_callret1_t *bcr, bcf_call_t *bc, mplp_pileup_t *gplp, htsFile *bcf_fp, bcf_hdr_t *bcf_hdr, bam_sample_t *sm, bam_hdr_t **h, mplp_ref_t *mp_ref, int* beg0, int* end0);
}

//Results of genotype call
struct genotype {
    bool is_het;
    //prob of het
    double p_het;
    genotype() {
        is_het = false;
        p_het = -1.0;
    }
    void determine_het() {
        if(p_het == -1.0) {
            throw runtime_error("Uninitialized genotype!");
        }
        if(p_het > MIN_HET_PROB) {
            is_het = true;
        }
    }
};

//results of mpileup for a variant
struct locus_info {
    //Flag for hom variants in RNA.
    bool is_hom_rna;
    //probability of hom genotype in RNA
    double p_hom_rna;
    //Flag for het variants in DNA.
    bool is_het_dna;
    //probability of het variant in DNA.
    double p_het_dna;
};

//Workhorse for "cis-ase identify"
class CisAseIdentifier {
    private:
        //VCF file with somatic variants
        string somatic_vcf_;
        //VCF file with polymorphisms
        string poly_vcf_;
        //Tumor RNAseq alignments
        string tumor_rna_;
        //Tumor DNA alignments
        string tumor_dna_;
        //Reference sequence FASTA
        string ref_;
        //GTF file with annotations
        string gtf_;
        //File to write output to
        string output_file_;
        //output stream to output annotated junctions file
        ofstream ofs_;
        //Somatic VCF file handle
        htsFile *somatic_vcf_fh_;
        //Somatic VCF Header
        bcf_hdr_t *somatic_vcf_header_;
        //Reference FASTA object
        faidx_t *ref_fai_;
        //Somatic VCF record
        bcf1_t *somatic_vcf_record_;
        //Configuration for somatic mpileup
        mplp_conf_t somatic_conf_;
        //Configuration for germline mpileup
        mplp_conf_t germline_conf_;
        //Get info about a variant - key is chr:start
        //Bi-allelic assumption
        map<string, locus_info> germline_variants_;
    public:
        //Constructor
        CisAseIdentifier() : somatic_vcf_("NA"),
                             tumor_rna_("NA"),
                             tumor_dna_("NA"), ref_("NA"), gtf_("NA"),
                             output_file_("NA"),
                             somatic_vcf_fh_(NULL),
                             somatic_vcf_header_(NULL),
                             somatic_vcf_record_(NULL) {}
        //Destructor
        ~CisAseIdentifier() {
            if(ofs_.is_open()) {
                ofs_.close();
            }
        }
        //Parse command line arguments
        void parse_options(int argc, char* argv[]);
        //Usage for this tool
        void usage(ostream &out);
        //The workhorse
        void run();
        //Open somatic VCF file
        void open_somatic_vcf();
        //Read in next record
        bool read_somatic_record();
        //Run mpileup and get the genotype likelihoods
        bool run_mpileup(string bam, mplp_conf_t *conf, bool (CisAseIdentifier::*f)(bcf_hdr_t*, int, int, const bcf_call_t&, bcf1_t*));
        //Call genotypes using the posterior prob
        genotype call_geno(const bcf_call_t& bc);
        //Get the SNPs within relevant window
        void process_snps_in_window(string window);
        //Process somatic variants
        bool process_somatic_het(bcf_hdr_t* bcf_hdr, int tid, int pos, const bcf_call_t& bc, bcf1_t* bcf_rec);
        bool process_germline_het(bcf_hdr_t* bcf_hdr, int tid, int pos, const bcf_call_t& bc, bcf1_t* bcf_rec);
        //Set the region as the somatic-vcf
        void set_mpileup_conf_somatic_vcf(mplp_conf_t &mplp_conf);
        //Set the region as the region-string
        void set_mpileup_conf_region(mplp_conf_t & mplp_conf, string region);
        //Free relevant pointers
        void cleanup();
        //Load pointer to reference
        void load_reference() {
            ref_fai_ = fai_load(ref_.c_str());
            if (ref_fai_ == NULL) throw runtime_error("Unable to open reference FASTA");
        }
        //Destroy pointer to reference
        void destroy_reference() {
            if(ref_fai_)
                fai_destroy(ref_fai_);
        }
        //Set the configuration for mpileup
        mplp_conf_t get_default_mpileup_conf() {
            mplp_conf_t mplp_conf;
            memset(&mplp_conf, 0, sizeof(mplp_conf_t));
            mplp_conf.min_baseQ = 13;
            mplp_conf.capQ_thres = 0;
            mplp_conf.max_depth = 250;
            mplp_conf.max_indel_depth = 250;
            mplp_conf.openQ = 40;
            mplp_conf.extQ = 20; mplp_conf.tandemQ = 100;
            mplp_conf.min_frac = 0.002; mplp_conf.min_support = 1;
            mplp_conf.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
            //uncompressed VCF
            mplp_conf.flag |= MPLP_BCF | MPLP_VCF | MPLP_NO_COMP;
            mplp_conf.fai = ref_fai_;
            mplp_conf.fai_fname = (char*)ref_.c_str();
            return mplp_conf;
        }
        void free_mpileup_conf(mplp_conf_t conf) {
            if (conf.bed)
                bed_destroy(conf.bed);
            if (conf.reg)
                free(conf.reg);
        }
};

#endif //CIS_ASE_IDENTIFIER_
