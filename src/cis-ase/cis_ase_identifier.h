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
#include "htslib/sam.h"
#include "bam2bcf.h"
#include "bam_plcmd.h"
#include "hts.h"
#include "vcf.h"
#include "sample.h"

using namespace std;

extern "C" {
    void *bed_read(const char *fn);
    void bed_destroy(void *_h);
    int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len);
    void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
            int n, char *const*fn, int *n_plp,
            const bam_pileup1_t **plp, int ignore_rg);
    int bed_overlap(const void *_h, const char *chr, int beg, int end);
    int mpileup_with_likelihoods(mplp_conf_t *conf, int n, char **fn, mplp_aux_t **data, bcf_callaux_t *bca, bcf_callret1_t *bcr, bcf_call_t *bc, mplp_pileup_t *gplp, htsFile *bcf_fp, bcf_hdr_t *bcf_hdr, bam_sample_t *sm, bam_hdr_t **h, mplp_ref_t *mp_ref);
}

//Workhorse for "cis-ase identify"
class CisAseIdentifier {
    private:
        //VCF file with somatic variants
        string somatic_vcf_;
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
        //Somatic VCF record
        bcf1_t *somatic_vcf_record_;
        //Configuration for mpileup
        mplp_conf_t mplp_conf_;
    public:
        //Constructor
        CisAseIdentifier() : somatic_vcf_("NA"),
                             tumor_rna_("NA"),
                             tumor_dna_("NA"), ref_("NA"), gtf_("NA"),
                             output_file_("NA"),
                             somatic_vcf_fh_(NULL),
                             somatic_vcf_header_(NULL),
                             somatic_vcf_record_(NULL) {
            somatic_vcf_record_ = bcf_init();
        }
        //Destructor
        ~CisAseIdentifier() {
            if(ofs_.is_open()) {
                ofs_.close();
            }
            bcf_destroy(somatic_vcf_record_);
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
        void run_mpileup();
        //Set the configuration for mpileup
        void set_mpileup_conf();
        //Free relevant pointers
        void cleanup();
};

#endif //CIS_ASE_IDENTIFIER_
