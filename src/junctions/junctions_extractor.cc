/*  junctions_extractor.h -- Declarations for `junctions extract` command

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

#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "common.h"
#include "junctions_extractor.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"

using namespace std;

//Parse the options passed to this tool
int JunctionsExtractor::parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    while((c = getopt(argc, argv, "ha:m:M:o:r:t:s:b:")) != -1) {
        switch(c) {
            case 'h':
                usage(help_ss);
                throw common::cmdline_help_exception(help_ss.str());
            case 'a':
                min_anchor_length_ = atoi(optarg);
                break;
            case 'm':
                min_intron_length_ = atoi(optarg);
                break;
            case 'M':
                max_intron_length_ = atoi(optarg);
                break;
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'r':
                region_ = string(optarg);
                break;
            case 't':
                strand_tag_ = string(optarg);
                break;
            case 's':
                strandness_ = atoi(optarg);
                break;
            case 'b':
                output_barcodes_file_ = string(optarg);
                break;
            case '?':
            default:
                usage();
                throw runtime_error("Error parsing inputs!(1)\n\n");
        }
    }
    if(argc - optind >= 1) {
        bam_ = string(argv[optind++]);
    }
    if(optind < argc || bam_ == "NA") {
        usage();
        throw runtime_error("Error parsing inputs!(2)\n\n");
    }
    if(strandness_ == -1){
        usage();
        throw runtime_error("Please supply strand specificity with '-s' option!\n\n");
    }
    cerr << "Minimum junction anchor length: " << min_anchor_length_ << endl;
    cerr << "Minimum intron length: " << min_intron_length_ << endl;
    cerr << "Maximum intron length: " << max_intron_length_ << endl;
    cerr << "Alignment: " << bam_ << endl;
    cerr << "Output file: " << output_file_ << endl;
    cerr << endl;
    return 0;
}

//Usage statement for this tool
int JunctionsExtractor::usage(ostream& out) {
    out << "Usage:" 
        << "\t\t" << "regtools junctions extract [options] indexed_alignments.bam" << endl;
    out << "Options:" << endl;
    out << "\t\t" << "-a INT\tMinimum anchor length. Junctions which satisfy a minimum \n"
        << "\t\t\t " << "anchor length on both sides are reported. [8]" << endl;
    out << "\t\t" << "-m INT\tMinimum intron length. [70]" << endl;
    out << "\t\t" << "-M INT\tMaximum intron length. [500000]" << endl;
    out << "\t\t" << "-o FILE\tThe file to write output to. [STDOUT]" << endl;
    out << "\t\t" << "-r STR\tThe region to identify junctions \n"
        << "\t\t\t " << "in \"chr:start-end\" format. Entire BAM by default." << endl;
    out << "\t\t" << "-s INT\tStrand specificity of RNA library preparation \n"
        << "\t\t\t " << "(0 = unstranded, 1 = first-strand/RF, 2, = second-strand/FR). REQUIRED" << endl;
    out << "\t\t" << "-t STR\tTag used in bam to label strand. [XS]" << endl;
        
    out << endl;
    return 0;
}

//Get the BAM filename
string JunctionsExtractor::get_bam() {
    return bam_;
}

//Name the junction based on the number of junctions
// in the map.
string JunctionsExtractor::get_new_junction_name() {
    int index = junctions_.size() + 1;
    stringstream name_ss;
    name_ss << "JUNC" << setfill('0') << setw(8) << index;
    return name_ss.str();
}

//Do some basic qc on the junction
bool JunctionsExtractor::junction_qc(Junction &j1) {
    if(j1.end - j1.start < min_intron_length_ ||
       j1.end - j1.start > max_intron_length_) {
        return false;
    }
    if(j1.start - j1.thick_start >= min_anchor_length_)
        j1.has_left_min_anchor = true;
    if(j1.thick_end - j1.end >= min_anchor_length_)
        j1.has_right_min_anchor = true;
    return true;
}

//Add a junction to the junctions map
//The read_count field is the number of reads supporting the junction.
int JunctionsExtractor::add_junction(Junction j1) {
    //Check junction_qc
    if(!junction_qc(j1)) {
        return 0;
    }

    //Construct key chr:start-end:strand
    stringstream s1;
    string start, end;
    s1 << j1.start; start = s1.str();
    s1 << j1.end; end = s1.str();
    //since ?,+,- sort differently on different systems
    string strand_proxy;
    if(j1.strand == "+") {
        strand_proxy = "0";
    } else if(j1.strand == "-") {
        strand_proxy = "1";
    } else {
        strand_proxy = "2";
    }
    string key = j1.chrom + string(":") + start + "-" + end + ":" + strand_proxy;

    //Check if new junction
    if(!junctions_.count(key)) {
        j1.name = get_new_junction_name();
        j1.read_count = 1;
        j1.score = common::num_to_str(j1.read_count);
    } else { //existing junction
        Junction j0 = junctions_[key];
        
        if (output_barcodes_file_ != "NA"){
            map<string, int>::const_iterator it = j0.barcodes.find(j1.barcodes.begin()->first);
            if (it != j0.barcodes.end()) {// barcode exists already
                j1.barcodes = j0.barcodes;
                j1.barcodes[it->first]++;
            } else {
                pair<string, int> tmp_barcode = *j1.barcodes.begin();
                j1.barcodes = j0.barcodes;
                j1.barcodes.insert(tmp_barcode);
            }
        }
        //increment read count
        j1.read_count = j0.read_count + 1;
        j1.score = common::num_to_str(j1.read_count);
        //Keep the same name
        j1.name = j0.name;
        //Check if thick starts are any better
        if(j0.thick_start < j1.thick_start)
            j1.thick_start = j0.thick_start;
        if(j0.thick_end > j1.thick_end)
            j1.thick_end = j0.thick_end;
        //preserve min anchor information
        j1.has_left_min_anchor = j1.has_left_min_anchor || j0.has_left_min_anchor;
        j1.has_right_min_anchor = j1.has_right_min_anchor || j0.has_right_min_anchor;
    }
    //Add junction and check anchor while printing.
    junctions_[key] = j1;
    return 0;
}

//Print all the junctions - this function needs work
vector<Junction> JunctionsExtractor::get_all_junctions() {
    //Sort junctions by position
    if(!junctions_sorted_) {
        create_junctions_vector();
        sort_junctions(junctions_vector_);
        junctions_sorted_ = true;
    }
    return junctions_vector_;
}

//Print all the junctions - this function needs work
void JunctionsExtractor::print_all_junctions(ostream& out) {
    ofstream fout;
    ofstream fout_barcodes;
    if(output_file_ != string("NA")) {
        fout.open(output_file_.c_str());
    }
    if(output_barcodes_file_!= string("NA")) {
        fout.open(output_barcodes_file_.c_str());
    }
    //Sort junctions by position
    if(!junctions_sorted_) {
        create_junctions_vector();
        sort_junctions(junctions_vector_);
        junctions_sorted_ = true;
    }
    for(vector<Junction> :: iterator it = junctions_vector_.begin();
        it != junctions_vector_.end(); it++) {
        Junction j1 = *it;
        if(j1.has_left_min_anchor && j1.has_right_min_anchor) {
            if(fout.is_open())
                j1.print(fout);
            else
                j1.print(out);
            if(fout_barcodes.is_open())
                j1.print_barcodes(fout_barcodes);
        }
    }
    if(fout.is_open())
        fout.close();
    if(fout_barcodes.is_open())
        fout_barcodes.close();
}

//Get the strand from the XS aux tag
void JunctionsExtractor::set_junction_strand_XS(bam1_t *aln, Junction& j1) {
    uint8_t *p = bam_aux_get(aln, strand_tag_.c_str());
    if(p != NULL) {
        char strand = bam_aux2A(p);
        strand ? j1.strand = string(1, strand) : j1.strand = string(1, '?');
        //cerr <<"XS strand is " << strand << endl;
    } else {
        //cerr <<"XS strand is NULL" << endl;
        j1.strand = string(1, '?');
        return;
    }
}

//Get the strand from the bitwise flag
void JunctionsExtractor::set_junction_strand_flag(bam1_t *aln, Junction& j1) {
    uint32_t flag = (aln->core).flag;
    int reversed = (flag >> 4) % 2;
    int mate_reversed = (flag >> 5) % 2;
    int first_in_pair = (flag >> 6) % 2;
    int second_in_pair = (flag >> 7) % 2;
    // strandness_ is 0 for unstranded, 1 for RF, and 2 for FR
    int bool_strandness = strandness_ - 1;
    int first_strand = !bool_strandness ^ first_in_pair ^ reversed;
    int second_strand = !bool_strandness ^ second_in_pair ^ mate_reversed;
    char strand;
    if (first_strand){
        strand = '+';
    } else {
        strand = '-';
    }
    //cerr << "flag is " << flag << endl;
    // if strand inferences from first and second in pair don't agree, we've got a problem
    if (first_strand == second_strand){
        j1.strand = string(1, strand);
    } else {
        j1.strand = string(1, '?');
    }
    //cerr <<"flag strand is " << j1.strand << endl;
    return;
}

//Get the strand
void JunctionsExtractor::set_junction_strand(bam1_t *aln, Junction& j1) {
    // if unstranded data
    if (strandness_ > 0){
        return set_junction_strand_flag(aln, j1);
    } else {
        return set_junction_strand_XS(aln, j1);
    }
}

//Get the the barcode
void JunctionsExtractor::set_junction_barcode(bam1_t *aln, Junction& j1) {
    uint8_t *p = bam_aux_get(aln, barcode_tag_.c_str());
    if(p != NULL) {
        char barcode = bam_aux2A(p);
        j1.barcodes.insert(pair<string, int>(string(1, barcode),1));
    } else {
        j1.barcodes.insert(pair<string, int>(string(1, '?'),1));
        cerr << 'WARNING: No CB tag found for alignment (id = ' << to_string(aln->id) << ')';
        return;
    }
}

//Parse junctions from the read and store in junction map
int JunctionsExtractor::parse_alignment_into_junctions(bam_hdr_t *header, bam1_t *aln) {
    int n_cigar = aln->core.n_cigar;
    if (n_cigar <= 1) // max one cigar operation exists(likely all matches)
        return 0;

    int chr_id = aln->core.tid;
    int read_pos = aln->core.pos;
    string chr(header->target_name[chr_id]);
    uint32_t *cigar = bam_get_cigar(aln);

    Junction j1;
    j1.chrom = chr;
    j1.start = read_pos; //maintain start pos of junction
    j1.thick_start = read_pos;
    set_junction_strand(aln, j1);
    if (output_barcodes_file_ != "NA"){
        set_junction_barcode(aln, j1);
    }
    bool started_junction = false;
    for (int i = 0; i < n_cigar; ++i) {
        char op =
               bam_cigar_opchr(cigar[i]);
        int len =
               bam_cigar_oplen(cigar[i]);
        switch(op) {
            case 'N':
                if(!started_junction) {
                    j1.end = j1.start + len;
                    j1.thick_end = j1.end;
                    //Start the first one and remains started
                    started_junction = true;
                } else {
                    //Add the previous junction
                    try {
                        add_junction(j1);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    j1.thick_start = j1.end;
                    j1.start = j1.thick_end;
                    j1.end = j1.start + len;
                    j1.thick_end = j1.end;
                    //For clarity - the next junction is now open
                    started_junction = true;
                }
                break;
            case '=':
            case 'M':
                if(!started_junction)
                    j1.start += len;
                else
                    j1.thick_end += len;
                break;
            //No mismatches allowed in anchor
            case 'D':
            case 'X':
                if(!started_junction) {
                    j1.start += len;
                    j1.thick_start = j1.start;
                } else {
                    try {
                        add_junction(j1);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    //Don't include these in the next anchor
                    j1.start = j1.thick_end + len;
                    j1.thick_start = j1.start;
                }
                started_junction = false;
                break;
            case 'I':
            case 'S':
                if(!started_junction)
                    j1.thick_start = j1.start;
                else {
                    try {
                        add_junction(j1);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    //Don't include these in the next anchor
                    j1.start = j1.thick_end;
                    j1.thick_start = j1.start;
                }
                started_junction = false;
                break;
            case 'H':
                break;
            default:
                cerr << "Unknown cigar " << op;
                break;
        }
    }
    if(started_junction) {
        try {
            add_junction(j1);
        } catch (const std::logic_error& e) {
            cout << e.what() << '\n';
        }
    }
    return 0;
}

//The workhorse - identifies junctions from BAM
int JunctionsExtractor::identify_junctions_from_BAM() {
    if(!bam_.empty()) {
        //open BAM for reading
        samFile *in = sam_open(bam_.c_str(), "r");
        if(in == NULL) {
            throw runtime_error("Unable to open BAM/SAM file.\n\n");
        }
        //Load the index
        hts_idx_t *idx = sam_index_load(in, bam_.c_str());
        if(idx == NULL) {
            throw runtime_error("Unable to open BAM/SAM index."
                                " Make sure alignments are indexed\n\n");
        }
        //Get the header
        bam_hdr_t *header = sam_hdr_read(in);
        //Initialize iterator
        hts_itr_t *iter = NULL;
        //Move the iterator to the region we are interested in
        iter  = sam_itr_querys(idx, header, region_.c_str());
        if(header == NULL || iter == NULL) {
            sam_close(in);
            throw runtime_error("Unable to iterate to region within BAM.\n\n");
        }
        //Initiate the alignment record
        bam1_t *aln = bam_init1();
        while(sam_itr_next(in, iter, aln) >= 0) {
            parse_alignment_into_junctions(header, aln);
        }
        hts_itr_destroy(iter);
        hts_idx_destroy(idx);
        bam_destroy1(aln);
        bam_hdr_destroy(header);
        sam_close(in);
    }
    return 0;
}

//Create the junctions vector from the map
void JunctionsExtractor::create_junctions_vector() {
    for(map<string, Junction> :: iterator it = junctions_.begin();
        it != junctions_.end(); it++) {
        Junction j1 = it->second;
        junctions_vector_.push_back(j1);
    }
}
