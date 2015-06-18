/*  gtf_parser.cc methods for gtf parsing

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

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm>
#include "bedFile.h"
#include "gtf_parser.h"
#include "lineFileUtilities.h"

using namespace std;

//Open the GTF file.
bool GtfParser::open() {
    gtf_fh_.open(gtffile_.c_str());
    if(!gtf_fh_.is_open()) {
        cerr << "\nUnable to open GTF file.";
        exit(1);
    }
    return true;
}

//Close the GTF filehandle
bool GtfParser::close() {
    if(gtf_fh_.is_open())
        gtf_fh_.close();
}

//parse an exon single line into a Gtf struct
Gtf GtfParser::parse_exon_line(string line) {
    Gtf gtf1;
    vector<string> fields;
    Tokenize(line, fields);
    assert(fields.size() == 9);
    if(fields[2] != "exon") {
        gtf1.is_exon = false;
        return gtf1;
    }
    gtf1.is_exon = true;
    gtf1.seqname = fields[0];
    gtf1.source = fields[1];
    gtf1.feature = fields[2];
    gtf1.start = atol(fields[3].c_str());
    gtf1.end = atol(fields[4].c_str());
    gtf1.score = fields[5];
    gtf1.strand = fields[6];
    gtf1.frame = fields[7][0];
    gtf1.attributes = fields[8];
    return gtf1;
}

//Parse the transcript name from attributes column
string parse_transcript_id(vector<string> attributes1) {
    for (int i = 0; i < attributes1.size(); i++) {
        vector<string> tokens;
        //some attributes have a leading whitespace
        if(attributes1[i][0] == ' ') {
            attributes1[i].erase(0, 1);
        }
        Tokenize(attributes1[i], tokens, ' ');
        if(tokens[0] == "transcript_id") {
            return tokens[1];
        }
    }
    return string("NA");
}

//Add an exon to the transcript map
bool GtfParser::add_exon_to_transcript_map(Gtf gtf1) {
    vector<string> attributes;
    Tokenize(gtf1.attributes, attributes, ';');
    string transcript_id = parse_transcript_id(attributes);
    //create a BED6 object
    BED exon = BED(gtf1.seqname, gtf1.start,
                   gtf1.end, gtf1.feature,
                   gtf1.score, gtf1.strand);
    if(transcript_id != "NA") {
        transcript_map_[transcript_id].exons.push_back(exon);
    }
    return true;
}

//Return vector of transcripts in a bin
vector<string> GtfParser::transcripts_from_bin(string chr, BIN bin1) {
    return chrbin_to_transcripts_[chr][bin1];
}

//Return the BIN that the transcript falls in
//This is formed by using the ends of the transcript
vector<BIN> GtfParser::bin_from_transcript(string transcript_id) {
    return transcript_to_bin_[transcript_id];
}

//Annotate each transcript with its bin
//Maintain a map that allows one to jump from chr,bin to a vector of
//transcript ids
//Also maintain another hash that allows one to jump from transcriptID
//to the bins that the exon-exon junctions for that transcript fall in.
bool GtfParser::annotate_transcript_with_bins() {
    if(!transcripts_sorted_) {
        sort_exons_within_transcripts();
    }
    for (std::map<string, Transcript>::iterator it = transcript_map_.begin(); \
            it != transcript_map_.end(); it++) {
        string transcript_id = it->first;
        vector<BED> & exons = (it->second).exons;
        cout << endl << transcript_id << exons.size();
        for (int i = 0; i< exons.size() - 1; i++) {
            string chr = exons[i].chrom;
            CHRPOS start = exons[i].end;
            CHRPOS end = exons[i + 1].start;
            BIN bin1 = getBin(start, end);
            chrbin_to_transcripts_[chr][bin1].push_back(transcript_id);
            transcript_to_bin_[transcript_id].push_back(bin1);
            cout << endl << transcript_id << "\t" << start << "\t" << end << bin1;
        }
    }
}

//Construct the junctions using exon information
bool GtfParser::construct_junctions() {
    if(!transcripts_sorted_) {
        sort_exons_within_transcripts();
    }
    for (std::map<string, Transcript>::iterator it = transcript_map_.begin(); \
            it != transcript_map_.end(); it++) {
        string transcript_id = it->first;
        vector<BED> exon_vector = it->second.exons;
        for (int i = 0; i < exon_vector.size() - 1; i++) {
           //Create the junction
           BED j1 =
               BED(exon_vector[i].chrom, exon_vector[i].end,
                   exon_vector[i+1].start);
           transcript_map_[transcript_id].junctions.push_back(j1);
        }
    }
}

//Sort the exons within transcripts by start position
bool GtfParser::sort_exons_within_transcripts() {
    for (std::map<string, Transcript>::iterator it = transcript_map_.begin(); \
            it != transcript_map_.end(); it++) {
        sort(it->second.exons.begin(), it->second.exons.end(), sortByStart);
    }
    transcripts_sorted_ = true;
}

//Print out transcripts - exons and junctions
bool GtfParser::print_transcripts() {
     for (std::map<string, Transcript>::iterator it = transcript_map_.begin(); \
          it != transcript_map_.end(); it++) {
            std::cout << it->first << " => \n";
            cout << "\tExons\n";
            vector<BED> & exons = (it->second).exons;
            for(vector<BED>::iterator it2 = exons.begin(); it2 != exons.end(); it2++) {
                cout << "\t" << it2->chrom << "\t" << it2-> start << "\t" << it2->end << "\n";
            }
            cout << "\tJunctions\n";
            vector<BED> & junctions = (it->second).junctions;
            for(vector<BED>::iterator it2 = junctions.begin(); it2 != junctions.end(); it2++) {
                cout << "\t" << it2->chrom << "\t" << it2-> start << "\t" << it2->end << "\n";
            }
     }
}

//Create a transcript map from the GTF
//This is a <map> data structure
//The key is transcript_id
//The value is a vector<BED> corresponding to exons
bool GtfParser::create_transcript_map() {
    if(!gtf_fh_.is_open()) {
        GtfParser::open();
    }
    string line;
    while(getline(gtf_fh_, line)) {
        Gtf gtf_l = parse_exon_line(line);
        if(gtf_l.is_exon) {
            add_exon_to_transcript_map(gtf_l);
            n_exons_++;
        }
    }
    GtfParser::close();
}

//Set the gtf file
bool GtfParser::set_gtffile(string filename) {
    gtffile_ = filename;
}

