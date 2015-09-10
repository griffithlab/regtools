/*  junctions_creator.h -- Declarations for `junctions create` command

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
#include "junctions_creator.h"
#include "sam.h"
#include "hts.h"
#include "faidx.h"
#include "kstring.h"

using namespace std;

//Parse the options passed to this tool
int JunctionsCreator::parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    while((c = getopt(argc, argv, "ha:o:r:")) != -1) {
        switch(c) {
            case 'a':
                min_anchor_length = atoi(optarg);
                break;
            case 'o':
                output_file = string(optarg);
                break;
            case 'r':
                region = string(optarg);
                break;
            case '?':
            case 'h':
            default:
                return 1;
        }
    }
    if(argc == optind) {
        cerr << endl << "Error parsing inputs!" << endl;
        return 1;
    }
    bam_ = string(argv[optind]);
    cerr << endl << "Minimum junction anchor length: " << min_anchor_length;
    cerr << endl << "BAM file: " << bam_;
    return 0;
}

//Usage statement for this tool
int JunctionsCreator::usage(ostream& out) {
    out << "\nUsage:\t\t" << "regtools junctions create [options] indexed_alignments.bam";
    out << "\nOptions:\t" << "-a INT\tMinimum anchor length on any one side of the junction. [8]";
    out << "\n\t" << "-o FILE\tThe file to write output to.";
    out << "\n\t" << "-r STR\tThe region to identify junctions in \"chr:start-end\" format.";
    out << "\n";
    return 0;
}

//Get the BAM filename
string JunctionsCreator::get_bam() {
    return bam_;
}

//Name the junction based on the number of junctions
// in the map.
string JunctionsCreator::get_new_junction_name() {
    int index = junctions.size() + 1;
    stringstream name_ss;
    name_ss << "JUNC" << setfill('0') << setw(8) << index;
    return name_ss.str();
}

//Do some basic qc on the junction
bool JunctionsCreator::junction_qc(Junction &j1) {
    if(j1.end - j1.start < min_intron_length ||
       j1.end - j1.start > max_intron_length) {
        cerr << "Failing qc length: " << j1.end - j1.start << endl;
        return false;
    }
    if(j1.start - j1.thick_start >= min_anchor_length)
        j1.has_left_min_anchor = true;
    if(j1.thick_end - j1.end >= min_anchor_length)
        j1.has_right_min_anchor = true;
    return true;
}

//Add a junction to the junctions map
//The read_count field is the number of reads supporting the junction.
int JunctionsCreator::add_junction(Junction j1) {
    //Check junction_qc
    if(!junction_qc(j1)) {
        cerr << endl << "Failed qc";
        return 0;
    }

    //Construct key chr:start-end:strand
    stringstream s1;
    string start, end;
    s1 << j1.start; start = s1.str();
    s1 << j1.end; end = s1.str();
    string key = j1.chrom + string(":") + start + "-" + end + ":" + j1.strand;

    //Check if new junction
    if(!junctions.count(key)) {
        j1.name = get_new_junction_name();
        j1.read_count = 1;
    } else { //existing junction
        Junction j0 = junctions[key];
        //increment read count
        j1.read_count = j0.read_count + 1;
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
    junctions[key] = j1;
    cerr << "Inside add_junction\n";
    print_one_junction(j1, cerr);
    return 0;
}

//Print one junction
void JunctionsCreator::print_one_junction(const Junction j1, ostream& out) {
    out << j1.chrom <<
        "\t" << j1.thick_start << "\t" << j1.thick_end <<
        "\t" << j1.name << "\t" << j1.read_count << "\t" << j1.strand <<
        "\t" << j1.thick_start << "\t" << j1.thick_end <<
        "\t" << j1.color << "\t" << j1.nblocks <<
        "\t" << j1.start - j1.thick_start << "," << j1.thick_end - j1.end <<
        "\t" << "0," << j1.end - j1.thick_start << endl;
}

//Print all the junctions - this function needs work
void JunctionsCreator::print_all_junctions(ostream& out) {
    ofstream fout;
    if(!output_file.empty())
        fout.open(output_file.c_str());
    for(map<string, Junction> :: iterator it = junctions.begin();
        it != junctions.end(); it++) {
        Junction j1 = it->second;
        if(j1.has_left_min_anchor && j1.has_right_min_anchor) {
            if(fout.is_open())
                print_one_junction(j1, fout);
            else
                print_one_junction(j1, out);
        }
    }
    if(fout.is_open())
        fout.close();
}

//Get the strand from the XS aux tag
void JunctionsCreator::set_junction_strand(bam1_t *aln, Junction& j1) {
    uint8_t *p = bam_aux_get(aln, "XS");
    if(p != NULL) {
        char strand = bam_aux2A(p);
        strand ? j1.strand = string(1, strand) : j1.strand = string(1, '?');
    } else {
        j1.strand = string(1, '?');
        return;
    }
}

//Parse junctions from the read and store in junction map
int JunctionsCreator::parse_alignment_into_junctions(bam_hdr_t *header, bam1_t *aln) {
    const bam1_core_t *c = &aln->core;
    if (c->n_cigar <= 1) // max one cigar operation exists(likely all matches)
        return 0;

    int chr_id = aln->core.tid;
    int read_pos = aln->core.pos;
    string chr(header->target_name[chr_id]);
    uint32_t *cigar = bam_get_cigar(aln);
    int n_cigar = c->n_cigar;

    Junction j1;
    j1.chrom = chr;
    j1.start = read_pos; //maintain start pos of junction
    j1.thick_start = read_pos;
    set_junction_strand(aln, j1);
    bool started_junction = false;
    cerr << "\nread_pos " << read_pos;
    for (int i = 0; i < n_cigar; ++i) {
        char op =
               bam_cigar_opchr(cigar[i]);
        int len =
               bam_cigar_oplen(cigar[i]);
        cerr << "\ncigar " << op << " " << len;
        switch(op) {
            //Add first junction if read overlaps
            // two junctions
            case 'N':
                if(!started_junction) {
                    j1.end = j1.start + len;
                    j1.thick_end = j1.end;
                    started_junction = true;
                } else {
                    cerr << endl << "DEBUG " << read_pos << "\t" <<
                        j1.start << "\t" << j1.end << "\t" <<
                        j1.thick_start << "\t" << j1.thick_end << endl;
                    //Add the previous junction
                    try {
                        add_junction(j1);
                    } catch (const std::logic_error& e) {
                        cout << e.what() << '\n';
                    }
                    j1.added = true;
                    //Start the next junction
                    j1.added = false;
                    j1.thick_start = j1.end;
                    j1.start = j1.thick_end;
                    j1.end = j1.start + len;
                    j1.thick_end = j1.end;
                    started_junction = true;
                }
                break;
            case 'D':
            case '=':
            case 'X':
            case 'M':
                if(!started_junction)
                    j1.start += len;
                else
                    j1.thick_end += len;
                break;
            //SEQ not in reference genome - skip
            case 'I':
            case 'S':
            case 'H':
                break;
            default:
                cerr << "Unknown cigar " << op;
                break;
        }
    }
    if(!j1.added && started_junction) {
        cerr << endl << "DEBUG2 " << read_pos << "\t" <<
            j1.start << "\t" << j1.end << "\t" << j1.thick_start << "\t" << j1.thick_end << endl;
        try {
            add_junction(j1);
        } catch (const std::logic_error& e) {
            cout << e.what() << '\n';
        }
        j1.added = true;
    }
    return 0;
}

//The workhorse - identifies junctions from BAM
int JunctionsCreator::identify_junctions_from_BAM() {
    if(!bam_.empty()) {
        cerr << endl << "Opening BAM " << bam_ << endl;
        //open BAM for reading
        samFile *in = sam_open(bam_.c_str(), "r");
        if(in == NULL) {
            return 1;
        }
        //Load the index
        hts_idx_t *idx = sam_index_load(in, bam_.c_str());
        if(idx == NULL) {
            cerr << "Unable to load alignment index. Please index with Samtools.";
            return -1;
        }
        //Get the header
        bam_hdr_t *header = sam_hdr_read(in);
        //Initialize iterator
        hts_itr_t *iter = NULL;
        //Move the iterator to the region we are interested in
        if(region.empty())
            region = "."; //Default = entire file
        iter  = sam_itr_querys(idx, header, region.c_str());
        if(header == NULL || iter == NULL) {
            sam_close(in);
            return 1;
        }
        //Initiate the alignment record
        bam1_t *aln = bam_init1();
        while(sam_itr_next(in, iter, aln) >= 0) {
            parse_alignment_into_junctions(header, aln);
        }
        hts_itr_destroy(iter);
        bam_destroy1(aln);
        bam_hdr_destroy(header);
        sam_close(in);
    }
    return 0;
}

