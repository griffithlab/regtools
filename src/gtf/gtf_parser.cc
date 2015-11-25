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
#include <stdexcept>
#include <vector>
#include <algorithm>
#include "common.h"
#include "bedFile.h"
#include "gtf_parser.h"
#include "lineFileUtilities.h"

using namespace std;

//Sort exons by start - positive strand
bool sort_by_start_ps(const BED & a, const BED & b) {
    return a.start < b.start;
}

//Sort exons by start - negative strand
bool sort_by_start_ns(const BED & a, const BED & b) {
    return a.start > b.start;
}

//Open the GTF file.
void GtfParser::open() {
    gtf_fh_.open(gtffile_.c_str());
    if(!gtf_fh_.is_open()) {
        cerr << "\nUnable to open GTF file.";
        exit(1);
    }
}

//Close the GTF filehandle
void GtfParser::close() {
    if(gtf_fh_.is_open())
        gtf_fh_.close();
}

//parse an exon single line into a Gtf struct
Gtf GtfParser::parse_exon_line(string line) {
    Gtf gtf1;
    vector<string> fields;
    Tokenize(line, fields);
    if(fields.size() != 9) {
        cerr << line << endl << fields.size();
        throw runtime_error("Expected 9 fields in GTF line.");
    }
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

//Parse the required field from attributes column
string GtfParser::parse_attribute(vector<string> attributes1,
                           string field_name) {
    for (std::size_t i = 0; i < attributes1.size(); i++) {
        vector<string> tokens;
        //some attributes have a leading whitespace
        if(attributes1[i][0] == ' ') {
            attributes1[i].erase(0, 1);
        }
        Tokenize(attributes1[i], tokens, ' ');
        if(tokens[0] == field_name) {
            unquote(tokens[1]);
            return tokens[1];
        }
    }
    return string("NA");
}

//Add an exon to the transcript map
void GtfParser::add_exon_to_transcript_map(Gtf gtf1) {
    vector<string> attributes;
    Tokenize(gtf1.attributes, attributes, ';');
    string transcript_id = parse_attribute(attributes, "transcript_id");
    string gene_name = parse_attribute(attributes, "gene_name");
    //create a BED6 object
    BED exon = BED(gtf1.seqname, gtf1.start,
                   gtf1.end, gtf1.feature,
                   gtf1.score, gtf1.strand);
    if(transcript_id != string("NA")) {
        transcript_map_[transcript_id].exons.push_back(exon);
        set_transcript_gene(transcript_id, gene_name);
    }
}

//Return the exons corresponding to a transcript
//The return value is a vector of BEDs
const vector<BED> & GtfParser::get_exons_from_transcript(string transcript_id) {
    return transcript_map_[transcript_id].exons;
}

//Return vector of transcripts in a bin
vector<string> GtfParser::transcripts_from_bin(string chr, BIN bin1) {
    return chrbin_to_transcripts_[chr][bin1];
}

//Return the BIN that the transcript falls in
//This is formed by using the ends of the transcript
BIN GtfParser::bin_from_transcript(string transcript_id) {
    return transcript_to_bin_[transcript_id];
}

//Annotate each transcript with its bin
//Maintain a map that allows one to jump from chr,bin to a vector of
//transcript ids
//Also maintain another hash that allows one to jump from transcriptID
//to the bin that the entire transcript falls in.
void GtfParser::annotate_transcript_with_bins() {
    //make sure exons are sorted
    if(!transcripts_sorted_) {
        sort_exons_within_transcripts();
    }
    for (std::map<string, Transcript>::iterator it = transcript_map_.begin(); \
            it != transcript_map_.end(); it++) {
        string transcript_id = it->first;
        vector<BED> & exons = (it->second).exons;
        string chr = exons[0].chrom;
        //start of first exon
        CHRPOS start = exons[0].start;
        //end of last exon
        CHRPOS end = exons[exons.size() - 1].end;
        BIN bin1 = getBin(start, end);
        chrbin_to_transcripts_[chr][bin1].push_back(transcript_id);
        transcript_to_bin_[transcript_id] = bin1;
    }
}

//Construct the junctions using exon information
void GtfParser::construct_junctions() {
    if(!transcripts_sorted_) {
        sort_exons_within_transcripts();
    }
    for (std::map<string, Transcript>::iterator it = transcript_map_.begin(); \
            it != transcript_map_.end(); it++) {
        string transcript_id = it->first;
        vector<BED> exon_vector = it->second.exons;
        for (std::size_t i = 0; i < exon_vector.size() - 1; i++) {
           //Create the junction
           BED j1 =
               BED(exon_vector[i].chrom, exon_vector[i].end,
                   exon_vector[i+1].start);
           transcript_map_[transcript_id].junctions.push_back(j1);
        }
    }
}

//Sort the exons within transcripts by start position
void GtfParser::sort_exons_within_transcripts() {
    for (std::map<string, Transcript>::iterator it = transcript_map_.begin(); \
            it != transcript_map_.end(); it++) {
        if(it->second.exons[0].strand == "+")
            sort(it->second.exons.begin(), it->second.exons.end(), sort_by_start_ps);
        else if(it->second.exons[0].strand == "-")
            sort(it->second.exons.begin(), it->second.exons.end(), sort_by_start_ns);
        else {
            cerr << "Undefined strand for exon ";
            cerr << it->second.exons[0].start << it->second.exons[0].end;
            exit(1);
        }
    }
    transcripts_sorted_ = true;
}

//Print out transcripts - exons and junctions
void GtfParser::print_transcripts() {
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
void GtfParser::create_transcript_map() {
    if(!gtf_fh_.is_open()) {
        GtfParser::open();
    }
    string line;
    while(getline(gtf_fh_, line)) {
        Gtf gtf_l = parse_exon_line(line);
        if(gtf_l.is_exon) {
            add_exon_to_transcript_map(gtf_l);
        }
    }
    GtfParser::close();
}

//Set the gtf file
void GtfParser::set_gtffile(string filename) {
    gtffile_ = filename;
}

//Get the gene ID using the trancript ID
string GtfParser::get_gene_from_transcript(string transcript_id) {
    if(transcript_to_gene_.count(transcript_id)) {
        return transcript_to_gene_[transcript_id];
    } else {
        return "NA";
    }
}

//Load all the necessary objects into memory
void GtfParser::load() {
    create_transcript_map();
    construct_junctions();
    sort_exons_within_transcripts();
    annotate_transcript_with_bins();
    //print_transcripts();
}

//Set the gene ID for a trancript ID
inline void GtfParser::set_transcript_gene(string transcript_id, string gene_id) {
    //check if key already exists
    if(transcript_to_gene_.count(transcript_id) == 0)
        transcript_to_gene_[transcript_id] = gene_id;
}
