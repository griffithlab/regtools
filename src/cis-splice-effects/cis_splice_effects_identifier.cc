/*  cis_splice_effects_identifier.cc -- 'cis-splice-effects identify' methods

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
#include <set>
#include "common.h"
#include "cis_splice_effects_identifier.h"
#include "junctions_annotator.h"
#include "junctions_extractor.h"
#include "variants_annotator.h"

//Usage for this tool
void CisSpliceEffectsIdentifier::usage(ostream& out) {
    out << "\nUsage:\t\t"
        << "regtools cis-splice-effects identify [options] variants.vcf"
        << " alignments.bam ref.fa annotations.gtf";
    out << "\nOptions:\t";
    out << "\t" << "-o Output file [STDOUT]";
    out << "\n";
}

//Do QC on files
void CisSpliceEffectsIdentifier::file_qc() {
    if(vcf_ == "NA" || bam_ == "NA" ||
       ref_ == "NA" || gtf_ == "NA") {
        usage(std::cout);
        throw runtime_error("\nError parsing inputs!(2)\n");
    }
    if(!common::file_exists(vcf_) || !common::file_exists(bam_) ||
       !common::file_exists(ref_) || !common::file_exists(gtf_)) {
        throw runtime_error("\nPlease make sure input files exist.\n");
    }
}

//Parse command line options
void CisSpliceEffectsIdentifier::parse_options(int argc, char* argv[]) {
    optind = 1; //Reset before parsing again.
    stringstream help_ss;
    char c;
    while((c = getopt(argc, argv, "o:w:h")) != -1) {
        switch(c) {
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'w':
                window_size_ = atoi(optarg);
                break;
            case 'h':
                usage(help_ss);
                throw common::cmdline_help_exception(help_ss.str());
            default:
                usage(std::cerr);
                throw runtime_error("\nError parsing inputs!(1)");
        }
    }
    if(argc - optind >= 2) {
        vcf_ = string(argv[optind++]);
        bam_ = string(argv[optind++]);
        ref_ = string(argv[optind++]);
        gtf_ = string(argv[optind++]);
    }
    if(optind < argc ||
       vcf_ == "NA" ||
       bam_ == "NA" ||
       ref_ == "NA" ||
       gtf_ == "NA"){
        usage(std::cout);
        throw runtime_error("\nError parsing inputs!(2)\n");
    }
    file_qc();
    cerr << "\nVariant file: " << vcf_;
    cerr << "\nAlignment file: " << bam_;
    cerr << "\nReference fasta file: " << ref_;
    cerr << "\nAnnotation file: " << gtf_;
    cerr << "\nWindow size: " << window_size_;
    cerr << endl;
}

void CisSpliceEffectsIdentifier::identify() {
    GtfParser gp1(gtf_);
    gp1.load();
    VariantsAnnotator va(vcf_, gtf_);
    va.open_vcf_in();
    va.set_gtf_parser(gp1);
    JunctionsAnnotator ja1(ref_, va.gtf());
    ja1.set_gtf_parser(gp1);
    set<Junction> unique_junctions;
    while(va.read_next_record()) {
        AnnotatedVariant v1 = va.annotate_record_with_transcripts(false);
        if(v1.annotation != non_splice_region_annotation) {
            string variant_region = v1.chrom + ":" +
                                    num_to_str(v1.start - window_size_) +
                                    "-" + num_to_str(v1.end + window_size_);
            std::cerr << variant_region << endl;
            JunctionsExtractor je1(bam_, variant_region);
            je1.identify_junctions_from_BAM();
            vector<Junction> junctions = je1.get_all_junctions();
            for (size_t i = 0; i < junctions.size(); i++) {
                unique_junctions.insert(junctions[i]);
            }
        }
    }
    for (set<Junction>::iterator j1 = unique_junctions.begin(); j1 != unique_junctions.end(); j1++) {
        AnnotatedJunction line(*j1);
        ja1.get_splice_site(line);
        ja1.annotate_junction_with_gtf(line);
        if(line.anchor != "DA") {
            line.print();
        }
    }
}
