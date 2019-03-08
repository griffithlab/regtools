/*  cis_splice_effects_associator.cc -- 'cis-splice-effects associate' methods

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
#include "cis_splice_effects_associator.h"
#include "junctions_annotator.h"
#include "junctions_extractor.h"
#include "variants_annotator.h"

//Usage for this tool
void CisSpliceEffectsAssociator::usage(ostream& out) {
    out << "Usage:" 
        << "\t\t" << "regtools cis-splice-effects associate [options] variants.vcf"
        << "\t\t " << "junctions.bed ref.fa annotations.gtf" << endl;
    out << "Options:" << endl;
    out << "\t\t" << "-o STR\tOutput file containing the aberrant splice junctions with annotations. [STDOUT]" << endl;
    out << "\t\t" << "-v STR\tOutput file containing variants annotated as splice relevant (VCF format)." << endl;
    out << "\t\t" << "-j STR\tOutput file containing the aberrant junctions in BED12 format." << endl;
    out << "\t\t" << "-s INT\tStrand specificity of RNA library preparation \n"
        << "\t\t\t " << "(0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR). [1]" << endl;
    out << "\t\t" << "-a INT\tMinimum anchor length. Junctions which satisfy a minimum \n"
        << "\t\t\t " << "anchor length on both sides are reported. [8]" << endl;
    out << "\t\t" << "-m INT\tMinimum intron length. [70]" << endl;
    out << "\t\t" << "-M INT\tMaximum intron length. [500000]" << endl;
    out << "\t\t" << "-w INT\tWindow size in b.p to identify splicing events in.\n" 
        << "\t\t\t " << "The tool identifies events in variant.start +/- w basepairs.\n"
        << "\t\t\t " << "Default behaviour is to look at the window between previous and next exons." << endl;
    out << "\t\t" << "-e INT\tMaximum distance from the start/end of an exon \n"
        << "\t\t\t " << "to annotate a variant as relevant to splicing, the variant \n"
        << "\t\t\t " << "is in exonic space, i.e a coding variant. [3]" << endl;
    out << "\t\t" << "-i INT\tMaximum distance from the start/end of an exon \n"
        << "\t\t\t " << "to annotate a variant as relevant to splicing, the variant \n"
        << "\t\t\t " << "is in intronic space. [2]" << endl;
    out << "\t\t" << "-I\tAnnotate variants in intronic space within a transcript(not to be used with -i)." << endl;
    out << "\t\t" << "-E\tAnnotate variants in exonic space within a transcript(not to be used with -e)." << endl;
    out << "\t\t" << "-S\tDon't skip single exon transcripts." << endl;
    out << endl;
}

//Return stream to write output to
void CisSpliceEffectsAssociator::close_ostream() {
    if(ofs_.is_open())
        ofs_.close();
}

//Return stream to write output to
//If output file is not empty, attempt to open
//If output file is empty, set to cout
void CisSpliceEffectsAssociator::set_ostream() {
    if(output_file_ == "NA") {
        common::copy_stream(cout, ofs_);
    } else {
        ofs_.open(output_file_.c_str());
        if(!ofs_.is_open())
            throw runtime_error("Unable to open " +
                                output_file_);
    }
    if(output_junctions_bed_ != "NA") {
        ofs_junctions_bed_.open(output_junctions_bed_.c_str());
        if(!ofs_junctions_bed_.is_open())
            throw runtime_error("Unable to open " +
                                output_junctions_bed_);
    }
}

//Do QC on files
void CisSpliceEffectsAssociator::file_qc() {
    if(vcf_ == "NA" || bed_ == "NA" ||
       ref_ == "NA" || gtf_ == "NA") {
        usage(std::cout);
        throw runtime_error("Error parsing inputs!(2)\n\n");
    }
    if(!common::file_exists(vcf_) || !common::file_exists(bed_) ||
       !common::file_exists(ref_) || !common::file_exists(gtf_)) {
        throw runtime_error("Please make sure input files exist.\n\n");
    }
}

//Parse command line options
void CisSpliceEffectsAssociator::parse_options(int argc, char* argv[]) {
    optind = 1; //Reset before parsing again.
    stringstream help_ss;
    char c;
    while((c = getopt(argc, argv, "o:w:v:j:e:Ei:ISha:m:M:")) != -1) {
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
            case 'j':
                output_junctions_bed_ = string(optarg);
                break;
            case 'i':
                intronic_min_distance_ = atoi(optarg);
                break;
            case 'e':
                exonic_min_distance_ = atoi(optarg);
                break;
            case 'I':
                all_intronic_space_ = true;
                break;
            case 'E':
                all_exonic_space_ = true;
                break;
            case 'S':
                skip_single_exon_genes_ = false;
                break;
            case 'h':
                usage(help_ss);
                throw common::cmdline_help_exception(help_ss.str());
            default:
                usage(std::cerr);
                throw runtime_error("Error parsing inputs!(1)\n\n");
        }
    }
    if(argc - optind >= 4) {
        vcf_ = string(argv[optind++]);
        bed_ = string(argv[optind++]);
        ref_ = string(argv[optind++]);
        gtf_ = string(argv[optind++]);
    }
    if(optind < argc ||
       vcf_ == "NA" ||
       bed_ == "NA" ||
       ref_ == "NA" ||
       gtf_ == "NA"){
        usage(std::cerr);
        throw runtime_error("Error parsing inputs!(2)\n\n");
    }
    file_qc();
    cerr << "Variant file: " << vcf_ << endl;
    cerr << "Junctions BED file: " << bed_ << endl;
    cerr << "Reference fasta file: " << ref_ << endl;
    cerr << "Annotation file: " << gtf_ << endl;
    if(window_size_ != 0) {
        cerr << "Window size: " << window_size_ << endl;
    }
    if(output_file_ != "NA")
        cerr << "Output file: " << output_file_ << endl;
    if(output_junctions_bed_ != "NA")
        cerr << "Output junctions BED file: " << output_junctions_bed_ << endl; 
    if(annotated_variant_file_ != "NA") {
        cerr << "Annotated variants file: " << annotated_variant_file_ << endl;
        write_annotated_variants_ = true;
    }
    cerr << endl;
}

//Call the junctions annotator; duplication of identify
void CisSpliceEffectsAssociator::annotate_junctions(const GtfParser& gp1) { 
    JunctionsAnnotator ja1(ref_, gp1);
    ja1.set_gtf_parser(gp1);
    set_ostream();
    //Annotate the junctions in the set and write to file
    AnnotatedJunction::print_header(ofs_, true);
    int i = 0;
    //This is ugly, waiting to start using C++11/14
    for (set<Junction>::iterator j1 = unique_junctions_.begin(); j1 != unique_junctions_.end(); j1++) {
        Junction j = *j1;
        AnnotatedJunction line(j);
        ja1.get_splice_site(line);
        ja1.annotate_junction_with_gtf(line);
        line.name = j.name = get_junction_name(++i);
        if(output_junctions_bed_ != "NA") {
            j.print(ofs_junctions_bed_);
        }
        line.variant_info = variant_set_to_string(junction_to_variant_[j]);
        line.print(ofs_, true);
    }
    close_ostream();
}

//get name for the junction; duplication of identify
string CisSpliceEffectsAssociator::get_junction_name(int i) {
    stringstream name_ss;
    name_ss << "JUNC" << setfill('0') << setw(8) << i;
    return name_ss.str();
}

//parse BED into vector of junctions
vector<Junction> CisSpliceEffectsAssociator::parse_BED_to_junctions(){
    vector<Junction> junctions;
    JunctionsAnnotator anno(bed_);
    AnnotatedJunction line;
    Junction junc;
    line.reset();
    try {
        anno.open_junctions();
        while(anno.get_single_junction(line)) {
            junc.chrom = line.fields[0];
            junc.start = line.fields[1];
            junc.end = line.fields[2];
            junc.name = line.fields[3];
            junc.read_count = line.fields[4];
            junc.strand = line.fields[5];
            junc.thick_start = line.fields[6];
            junc.thick_end = line.fields[7];
            junc.color = line.fields[8];
            junc.nblocks = line.fields[9];
            junc.has_left_min_anchor = true;
            junc.has_right_min_anchor = true;
            junctions.push_back(junc); // are we making a copy or just adding a pointer?
            line.reset();
        }
        anno.close_junctions();
    } catch(const common::cmdline_help_exception& e) {
        cerr << e.what() << endl;
        return 0;
    } catch(const runtime_error& e) {
        cerr << e.what() << endl;
        return 1;
    }
    return junctions;
}

//The workhorse
void CisSpliceEffectsAssociator::associate() {
    //GTF parser object
    GtfParser gp1(gtf_);
    gp1.load();
    //variant annotator
    VariantsAnnotator va(vcf_, gp1, annotated_variant_file_, intronic_min_distance_, exonic_min_distance_, all_intronic_space_, all_exonic_space_, skip_single_exon_genes_);
    va.open_vcf_in();
    if(write_annotated_variants_)
        va.open_vcf_out();
    cerr << endl;
    //Get junctions from BED
    vector<Junction> junctions = parse_BED_to_junctions();
    //Annotate each variant and pay attention to splicing related ones
    while(va.read_next_record()) {
        AnnotatedVariant v1 = va.annotate_record_with_transcripts();
        if(v1.annotation != non_splice_region_annotation_string) {
            string region_start = window_size_ ? common::num_to_str(v1.start - window_size_) :
                                           common::num_to_str(v1.cis_effect_start);
            string region_end = window_size_ ? common::num_to_str(v1.end + window_size_) :
                                           common::num_to_str(v1.cis_effect_end);
            string variant_region = v1.chrom + ":" + region_start + "-" + region_end;
            cerr << "Variant " << v1;
            cerr << "Variant region is " << variant_region << endl;
            cerr << endl;
            if(write_annotated_variants_)
                va.write_annotation_output(v1);
            //Add all the junctions to the unique set
            for (size_t i = 0; i < junctions.size(); i++) {
                //Allow partial overlap - either junction start or end is within window
                if((junctions[i].start >= v1.cis_effect_start && junctions[i].start <= v1.cis_effect_end) ||
                   (junctions[i].end <= v1.cis_effect_end && junctions[i].end >= v1.cis_effect_start)) {
                   unique_junctions_.insert(junctions[i]);
                   //add to the map of junctions to variants
                   junction_to_variant_[junctions[i]].insert(v1);
                }
            }
        }
    }
    annotate_junctions(gp1);
}
