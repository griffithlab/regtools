/*  variants_annotator.h -- Declarations for `variants annotate`

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

#ifndef VARIANTS_ANNOTATOR_H_
#define VARIANTS_ANNOTATOR_H_

#include <iostream>
#include <stdint.h>
#include "bedFile.h"
#include "gtf_parser.h"
#include "hts.h"
#include "junctions_annotator.h"
#include "vcf.h"

using namespace std;

//Hold annotations
struct AnnotatedVariant : public BED {
    string overlapping_genes;
    string overlapping_transcripts;
    string overlapping_distances;
    string annotation;
    AnnotatedVariant() : overlapping_genes("NA"),
                         overlapping_transcripts("NA"),
                         overlapping_distances("NA"){}
    AnnotatedVariant(string chr1, CHRPOS start1, CHRPOS end1):
                         BED(chr1, start1, end1),
                         overlapping_genes("NA"),
                         overlapping_transcripts("NA"),
                         overlapping_distances("NA") {}
};

//The class that does all the annotation
class VariantsAnnotator {
    private:
        //Variant file
        string vcf_;
        //Gene annotations file
        string gtffile_;
        //GTF parser object - holds the GTF in memory
        GtfParser gtf_;
        //Output VCF file
        string vcf_out_;
        //Minimum distance of a variant from
        //edge of an exon(intronic) to be considered
        //a splicing variant
        uint32_t intronic_min_distance_;
        //Minimum distance of a variant from
        //edge of an exon(exonic) to be considered
        //a splicing variant
        uint32_t exonic_min_distance_;
        //Option to skip single exon genes
        bool skip_single_exon_genes_;
        //VCF file handle
        htsFile *vcf_fh_in_;
        //Header of VCF file
        bcf_hdr_t *vcf_header_in_;
        //Output VCF file handle
        htsFile *vcf_fh_out_;
        //Header of output VCF file
        bcf_hdr_t *vcf_header_out_;
        //Each VCF record
        bcf1_t *vcf_record_;
    public:
        //Default constructor
        VariantsAnnotator() : vcf_("NA"), gtffile_("NA"),
                              vcf_out_("NA"),
                              intronic_min_distance_(2),
                              exonic_min_distance_(3),
                              skip_single_exon_genes_(true),
                              vcf_fh_in_(NULL), vcf_header_in_(NULL),
                              vcf_fh_out_(NULL), vcf_header_out_(NULL),
                              vcf_record_(NULL) {
            vcf_record_ = bcf_init();
        }
        //constructor
        VariantsAnnotator(string vcf_f, string gtf_f, string vcf_out) : vcf_(vcf_f),
                              gtffile_(gtf_f),
                              vcf_out_(vcf_out),
                              intronic_min_distance_(2),
                              exonic_min_distance_(3),
                              skip_single_exon_genes_(true),
                              vcf_fh_in_(NULL), vcf_header_in_(NULL),
                              vcf_fh_out_(NULL), vcf_header_out_(NULL),
                              vcf_record_(NULL) {
            vcf_record_ = bcf_init();
            gtf_.set_gtffile(gtffile_);
        }
        //Destructor
        ~VariantsAnnotator() {
            cleanup_vcf();
        }
        //Parse command-line options for this tool
        int parse_options(int argc, char *argv[]);
        //Usage statement for this tool
        int usage(ostream& out);
        //Annotate VCF file
        void annotate_vcf();
        //Read in GTF file
        void load_gtf();
        //Open input VCF file
        void open_vcf_in();
        //Open output VCF file
        void open_vcf_out();
        //Cleanup VCF file data structures
        void cleanup_vcf();
        //Set the GTF parser
        void set_gtf_parser(GtfParser gp1) {
            gtf_ = gp1;
        }
        //Return GTF parser
        GtfParser gtf() {
            return gtf_;
        }
        //Annotate one line of a VCF
        AnnotatedVariant annotate_record_with_transcripts();
        //Given a transcript ID and variant position,
        //check if the variant is in a splice relevant region
        //relevance depends on the user params
        //intronic_min_distance_ and exonic_min_distance_
        //stores result in the AnnotatedVariant object
        void get_variant_overlaps_spliceregion(const vector<BED> &exons,
                                           AnnotatedVariant  &variant);
        //Read next record of VCF.
        bool read_next_record();
        //Write annotation output
        void write_annotation_output(const AnnotatedVariant &v1);
};

#endif
