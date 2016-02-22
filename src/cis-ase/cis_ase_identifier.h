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
#include "hts.h"
#include "vcf.h"

using namespace std;

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
        //Print somatic VCF record
        void print_somatic_record();
        //Free relevant pointers
        void cleanup();
};

#endif //CIS_ASE_IDENTIFIER_
