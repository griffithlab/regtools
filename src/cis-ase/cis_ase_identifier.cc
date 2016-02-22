/*  cis_ase_identifier.cc -- 'cis-ase identify' methods

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

#include <stdexcept>
#include "cis_ase_identifier.h"
#include "hts.h"

using namespace std;

//Usage for this tool
void CisAseIdentifier::usage(ostream& out) {
    out << "\nUsage:\t\t"
        << "regtools cis-splice-effects identify [options] somatic_variants.vcf"
        << " tumor_dna_alignments.bam tumor_rna_alignments.bam ref.fa annotations.gtf";
    out << "\nOptions:";
    out << "\t"   << "-o STR Output file containing the aberrant splice junctions with annotations. [STDOUT]";
    out << "\n\t\t" << "-v STR Output file containing variants annotated as splice relevant (VCF format).";
    out << "\n\t\t" << "-w INT\tWindow size in b.p to identify splicing events in. "
        << "\n\t\t\t" << "The tool identifies events in variant.start +/- w basepairs."
        << "\n\t\t\t" << "Default behaviour is to look at the window between previous and next exons.";
    out << "\n\t\t" << "-j STR Output file containing the aberrant junctions in BED12 format.";
    out << "\n";
}

//Parse command line options
void CisAseIdentifier::parse_options(int argc, char* argv[]) {
    optind = 1; //Reset before parsing again.
    char c;
    while((c = getopt(argc, argv, "o:w:v:j:h")) != -1) {
        switch(c) {
            case 'o':
                output_file_ = string(optarg);
                break;
            default:
                usage(std::cerr);
                throw runtime_error("\nError parsing inputs!(1)");
        }
    }
    if(argc - optind >= 5) {
        somatic_vcf_ = string(argv[optind++]);
        tumor_dna_ = string(argv[optind++]);
        tumor_rna_ = string(argv[optind++]);
        ref_ = string(argv[optind++]);
        gtf_ = string(argv[optind++]);
    }
    if(optind < argc ||
       somatic_vcf_ == "NA" ||
       tumor_dna_ == "NA" ||
       tumor_rna_ == "NA" ||
       ref_ == "NA" ||
       gtf_ == "NA"){
        usage(std::cerr);
        throw runtime_error("\nError parsing inputs!(2)\n");
    }
    cerr << "\nSomatic variants: " << somatic_vcf_;
    cerr << "\nTumor DNA: " << tumor_dna_;
    cerr << "\nTumor RNA: " << tumor_rna_;
    cerr << "\nReference fasta file: " << ref_;
    cerr << "\nAnnotation file: " << gtf_;
    cerr << endl;
}

//Open somatic VCF file
void CisAseIdentifier::open_somatic_vcf() {
    somatic_vcf_fh_ = bcf_open(somatic_vcf_.c_str(), "r");
    if(somatic_vcf_fh_ == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    somatic_vcf_header_ = bcf_hdr_read(somatic_vcf_fh_);
    if(somatic_vcf_header_ == NULL) {
        throw std::runtime_error("Unable to read header.");
    }
}

//Read in next record
bool CisAseIdentifier::read_somatic_record() {
    return (bcf_read(somatic_vcf_fh_, somatic_vcf_header_, somatic_vcf_record_) == 0);
}

//Print somatic VCF record
void CisAseIdentifier::print_somatic_record() {
    if(somatic_vcf_record_) {
        string chr = std::string(bcf_hdr_id2name(somatic_vcf_header_, somatic_vcf_record_->rid));
        cout << chr << "\t" << somatic_vcf_record_->pos << "\n";
    }
}

//Free relevant pointers
void CisAseIdentifier::cleanup() {
    if(somatic_vcf_header_)
        bcf_hdr_destroy(somatic_vcf_header_);
    if(somatic_vcf_fh_)
        bcf_close(somatic_vcf_fh_);
    if(somatic_vcf_record_)
        bcf_destroy(somatic_vcf_record_);
}

//The workhorse
void CisAseIdentifier::run() {
    cout << 1 << endl;
    open_somatic_vcf();
    cout << 2 << endl;
    while(read_somatic_record()) {
        print_somatic_record();
    }
    cout << 3 << endl;
    cleanup();
}

/*    cerr << endl << "In Run";
    for each_dna_somatic_variant {
        if variant_is_het_in_tumor_dna {
            for each_nearby_germline_polymorphism {
                if poly_is_het_tumor_dna and poly_is_not_het_in_tumor_rna {
                    print dna_somatic_variant "\t" dna_poly "\t" rna_poly "\t" evidence_level"
                    add_poly_to_list_of_ase_polys();
                }
            }
        }
    }
*/
