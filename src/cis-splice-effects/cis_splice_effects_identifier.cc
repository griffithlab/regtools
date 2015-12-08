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
    out << "\t" << "-v Annotated variants(optional, this is in the VCF format)";
    out << "\n";
}

//Return stream to write output to
void CisSpliceEffectsIdentifier::close_ostream() {
    if(ofs_.is_open())
        ofs_.close();
}

//Return stream to write output to
//If output file is not empty, attempt to open
//If output file is empty, set to cout
void CisSpliceEffectsIdentifier::set_ostream() {
    if(output_file_ == "NA") {
        common::copy_stream(cout, ofs_);
        return;
    }
    ofs_.open(output_file_.c_str());
    if(!ofs_.is_open())
        throw runtime_error("Unable to open " +
                            output_file_);
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
    while((c = getopt(argc, argv, "o:v:w:h")) != -1) {
        switch(c) {
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'w':
                window_size_ = atoi(optarg);
                break;
            case 'v':
                annotated_variant_file_ = string(optarg);
                break;
            case 'h':
                usage(help_ss);
                throw common::cmdline_help_exception(help_ss.str());
            default:
                usage(std::cerr);
                throw runtime_error("\nError parsing inputs!(1)");
        }
    }
    if(argc - optind >= 4) {
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
    if(window_size_ != 0) {
        cerr << "\nWindow size: " << window_size_;
    }
    if(output_file_ != "NA")
        cerr << "\nOutput file: " << output_file_;
    if(annotated_variant_file_ != "NA") {
        cerr << "\nAnnotated variants file: " << annotated_variant_file_;
        write_annotated_variants_ = true;
    }
    cerr << endl;
}

//The workhorse
void CisSpliceEffectsIdentifier::identify() {
    //GTF parser object
    GtfParser gp1(gtf_);
    gp1.load();
    //variant annotator
    VariantsAnnotator va(vcf_, gtf_, annotated_variant_file_);
    va.open_vcf_in();
    if(write_annotated_variants_)
        va.open_vcf_out();
    va.set_gtf_parser(gp1);
    JunctionsAnnotator ja1(ref_, va.gtf());
    ja1.set_gtf_parser(gp1);
    //Unique set of junctions near splicing variants
    set<Junction> unique_junctions;
    //Annotate each variant and pay attention to splicing related ones
    while(va.read_next_record()) {
        AnnotatedVariant v1 = va.annotate_record_with_transcripts();
        if(v1.annotation != non_splice_region_annotation_string) {
            string region_start = window_size_ ? common::num_to_str(v1.start - window_size_) :
                                           common::num_to_str(v1.cis_effect_start);
            string region_end = window_size_ ? common::num_to_str(v1.end + window_size_) :
                                           common::num_to_str(v1.cis_effect_end);
            string variant_region = v1.chrom + ":" + region_start + "-" + region_end;
            cerr << "\n\nVariant " << v1;
            cerr << "Variant region is " << variant_region;
            if(write_annotated_variants_)
                va.write_annotation_output(v1);
            //Extract junctions near this variant
            JunctionsExtractor je1(bam_, variant_region);
            je1.identify_junctions_from_BAM();
            vector<Junction> junctions = je1.get_all_junctions();
            //Add all the junctions to the unique set
            for (size_t i = 0; i < junctions.size(); i++) {
                if(window_size_ == 0) {
                    if(junctions[i].start >= v1.cis_effect_start &&
                       junctions[i].end <= v1.cis_effect_end) {
                       unique_junctions.insert(junctions[i]);
                    }
                    continue;
                }
                if(common::coordinate_diff(junctions[i].start, v1.start) < window_size_ &&
                   common::coordinate_diff(junctions[i].end, v1.start) <= window_size_) {
                       unique_junctions.insert(junctions[i]);
                }
            }
        }
    }
    set_ostream();
    //Annotate the junctions in the set and write to file
    AnnotatedJunction::print_header(ofs_);
    //This is ugly, waiting to start using C++11/14
    for (set<Junction>::iterator j1 = unique_junctions.begin(); j1 != unique_junctions.end(); j1++) {
        AnnotatedJunction line(*j1);
        ja1.get_splice_site(line);
        ja1.annotate_junction_with_gtf(line);
        if(line.anchor != "DA") {
            line.print(ofs_);
        }
    }
    close_ostream();
}
