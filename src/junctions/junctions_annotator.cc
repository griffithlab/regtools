/*  junctions_annotator.cc -- class methods for `junctions annotate`

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
#include <string>
#include "common.h"
#include "junctions_annotator.h"
#include "faidx.h"

using namespace std;

//Open junctions file
void JunctionsAnnotator::open_junctions() {
    junctions_.Open();
}

//Open junctions file
void JunctionsAnnotator::close_junctions() {
    junctions_.Close();
}

//Adjust the start and end of the junction
void JunctionsAnnotator::adjust_junction_ends(BED & line) {
    //Adjust the start and end with block sizes
    //The junction start is start + block_size1
    //The junction end is end + block_size1
    if(line.fields[10].empty()) {
        cerr << "\nBlock sizes not found. Invalid line.";
        cerr << "\n\tChr " << line.chrom << "\tPos " << line.start;
    }
    string blocksize_field = line.fields[10];
    vector<int> block_sizes;
    Tokenize(blocksize_field, block_sizes, ',');
    line.start += block_sizes[0];
    line.end -= block_sizes[1] - 1;
}

//Get a single line from the junctions file
void JunctionsAnnotator::get_single_junction(BED & line) {
    junctions_.GetNextBed(line);
    adjust_junction_ends(line);
}

//Get the anchor bases
bool JunctionsAnnotator::get_anchor_bases(const BED & line) {
    string position1 = line.chrom + ":" +
                      num_to_str(line.start + 1) + "-" + num_to_str(line.start + 2);
    string position2 = line.chrom + ":" +
                      num_to_str(line.start - 1) + "-" + num_to_str(line.start);
    string seq1, seq2;
    if(line.strand == "-") {
        seq1 = get_reference_sequence(position1);
        seq2 = get_reference_sequence(position2);
        seq1 = rev_comp(seq1);
        seq2 = rev_comp(seq2);
        cout << "\t" << seq2 + "-" + seq1;
    } else {
        cout << "\t" << seq1 + "-" + seq2;
    }
}

//Annotate a single line
bool JunctionsAnnotator::annotate_single_line(BED & line) {
    cout << endl << line.chrom << "\t" << line.start <<
            "\t" << line.end << "\t" << line.strand <<
            "\t" << line.name << "\t" << line.score;
    get_anchor_bases(line);
    cout << endl;
    return true;
}

//Extract gtf info
bool JunctionsAnnotator::read_gtf() {
    //the parse_options method sets the gtf filename
    gtf_.create_transcript_map();
    gtf_.construct_junctions();
    gtf_.sort_exons_within_transcripts();
    gtf_.annotate_transcript_with_bins();
    //gtf_.print_transcripts();
    return true;
}

//Annotate with gtf
//Takes a single junction BED and annotates with GTF
void JunctionsAnnotator::annotate_junction_with_gtf(BED j1) {
    BIN junction_bin = getBin(j1.start, j1.end);

    //From BedTools
    BIN start_bin, end_bin;
    start_bin = (j1.start >> _binFirstShift);
    end_bin = ((j1.end-1) >> _binFirstShift);
    CHRPOS j1_length = (j1.end - j1.start);

    cout << endl << "junction_bin " << junction_bin;
    /*
      SYNOPSIS:
         1. We loop through each UCSC BIN level for feature A's chrom.
         2. For each BIN, we loop through each B feature and add it to
            hits if it meets all of the user's requests.
    */
    cout << "BINS\n";
    for (BINLEVEL i = 0; i < _binLevels; ++i) {
        cout << endl << "BINLEVEL " << i;
        BIN offset = _binOffsetsExtended[i];
        //rename j to b
        cout << endl << "start_bin " << start_bin;
        cout << endl << "end_bin " << end_bin;
        for (BIN j = (start_bin + offset); j <= (end_bin + offset); ++j) {
            vector<string> transcripts = gtf_.transcripts_from_bin(j1.chrom, j);
            if(transcripts.size())
                for(int i = 0; i < transcripts.size(); i++)
                    cout << transcripts[i] << endl;
        }
        start_bin >>= _binNextShift;
        end_bin >>= _binNextShift;
    }
    vector<BIN> binv = gtf_.bin_from_transcript("ENST00000447898");
    cout << "\nBINS for this transcript.";
    for(int i = 0; i < binv.size(); i++) {
        cout << "\t" << binv[i];
    }

    cout << endl;
}

//Get the reference sequence at a particular coordinate
string JunctionsAnnotator::get_reference_sequence(string position) {
    int len;
    faidx_t *fai = fai_load(ref_.c_str());
    char *s = fai_fetch(fai, position.c_str(), &len);
    std::string seq(s);
    free(s);
    fai_destroy(fai);
    return seq;
}

//Get the name of the GTF file
string JunctionsAnnotator::gtf_file() {
    return gtf_.gtffile();
}

//Parse the options passed to this tool
int JunctionsAnnotator::parse_options(int argc, char *argv[]) {
    static struct option long_options[] = {
        {"ref", required_argument, 0, 'r'},
        {"junctions", required_argument, 0, 'j'},
        {"gtf", required_argument, 0, 'g'},
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "r:g:j:",
                    long_options, &option_index);
    while(c != -1) {
        switch(c) {
            case 'r':
                ref_ = string(optarg);
                break;
            case 'g':
                gtf_.set_gtffile(string(optarg));
                break;
            case 'j':
                junctions_.bedFile = string(optarg);
                break;
            default:
                usage();
                exit(-1);
        }
        c = getopt_long(argc, argv, "r:g:j:",
                    long_options, &option_index);
    }
    if(optind < argc || ref_.empty() || junctions_.bedFile.empty() || gtf_.gtffile().empty()) {
        usage();
        exit(-1);
    }
    cerr << "\nReference: " << ref_;
    cerr << "\nGTF: " << gtf_.gtffile();
    cerr << "\nJunctions: " << junctions_.bedFile;
}

//Usage statement for this tool
int JunctionsAnnotator::usage() {
    cout << "\nUsage:\t\t" << "regtools junctions annotate [options] -r ref.fa -j junctions.bed -g annotations.gtf";
    cout << "\nCommand:\t" << "annotate\tAnnotate the junctions.";
    cout << "\n";
    return 0;
}


