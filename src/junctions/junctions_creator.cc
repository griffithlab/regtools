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
#include <sstream>
#include "junctions_creator.h"
#include "sam.h"
#include "faidx.h"
#include "kstring.h"

using namespace std;

//Parse the options passed to this tool
int JunctionsCreator::parse_options(int argc, char *argv[]) {
    static struct option long_options[] = {
        {"bam", required_argument, 0, 'b'},
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "b:",
                    long_options, &option_index);
    while(c != -1) {
        switch(c) {
            case 'b':
                bam_ = string(optarg);
                break;
            default:
                usage();
                exit(-1);
        }
        c = getopt_long(argc, argv, "b:",
                    long_options, &option_index);
    }
    if(optind < argc || bam_.empty()) {
        usage();
        exit(-1);
    }
    cerr << "\nBAM: " << bam_;
}

//Usage statement for this tool
int JunctionsCreator::usage() {
    cout << "\nUsage:\t\t" << "regtools junctions create [options] -b alignments.bam";
    cout << "\nOptions:\t" << "-b alignments.bam\tThe alignment file to extract junctions from.";
    cout << "\n";
    return 0;
}

//Add a junction to the junctions map
//The read_count field is the number of reads supporting the junction.
int JunctionsCreator::add_junction(Junction j1) {
    cerr << "Adding junction\n";
    cerr << "\nChr " << j1.chrom << "\tStart " << j1.start << "\tEnd " << j1.end;
    stringstream s1;
    string start, end;
    s1 << j1.start; start = s1.str();
    s1 << j1.end; end = s1.str();
    string key = j1.chrom + string(":") + start + "-" + end;
    if(!junctions.count(key)) {
        j1.read_count = 1;
    } else {
        j1 = junctions[key];
        j1.read_count++;
    }
    junctions[key] = j1;
    return 0;
}

//Print all the junctions
void JunctionsCreator::print_junctions() {
    for(map<string, Junction> :: iterator it = junctions.begin();
        it != junctions.end(); it++) {
        Junction j1 = it->second;
        cout << endl << j1.chrom << "\t" << j1.start << "\t" <<
                j1.end << "\t" << j1.read_count;
    }
}

//Parse junctions from the read and store in junction map
int JunctionsCreator::parse_cigar_into_junctions(string chr, int read_pos,
                                                 uint32_t *cigar, int n_cigar) {
    Junction j1;
    j1.chrom = chr;
    j1.start = read_pos; //maintain start pos of junction
    for (int i = 0; i < n_cigar; ++i) {
        char op =
               bam_cigar_opchr(cigar[i]);
        int len =
               bam_cigar_oplen(cigar[i]);
        switch(op) {
            case 'N':
                j1.end = j1.start + len;
                add_junction(j1);
                //Reset for next junction - long read case
                j1.start = j1.end;
                break;
            case 'D':
            case '=':
            case 'X':
            case 'M':
                j1.start += len;
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
}

//Pull out the cigar string from the read
int JunctionsCreator::parse_read(bam_hdr_t *header, bam1_t *aln) {
    const bam1_core_t *c = &aln->core;
    if (c->n_cigar) { // cigar
        int chr_id = aln->core.tid;
        int read_pos = aln->core.pos;
        string chr(header->target_name[chr_id]);
        uint32_t *cigar = bam_get_cigar(aln);
        int n_cigar = c->n_cigar;
        parse_cigar_into_junctions(chr, read_pos, cigar, n_cigar);
    }
    return 0;
}

//The workhorse - identifies junctions from BAM
int JunctionsCreator::identify_junctions_from_BAM() {
    if(!bam_.empty()) {
        cerr << "BAM is " << bam_;
        //open BAM for reading
        samFile *in = sam_open(bam_.c_str(), "r");
        //Get the header
        bam_hdr_t *header = sam_hdr_read(in);
        //Initiate the alignment record
        bam1_t *aln = bam_init1();
        while(sam_read1(in, header, aln) >= 0) {
            parse_read(header, aln);
        }
        bam_destroy1(aln);
        bam_hdr_destroy(header);
        sam_close(in);
    }
    return 0;
}

