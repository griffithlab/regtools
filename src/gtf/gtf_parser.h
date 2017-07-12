/*  gtf_parser.h -- interface for gtf parsing

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
#ifndef GTF_PARSER_H_
#define GTF_PARSER_H_

#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include "bedFile.h"
#include "lineFileUtilities.h"

using namespace std;

//BedVector to hold BED entries.
typedef vector<BED> BedVector;

//Each transcript is a map
//key is the transcript_id
//value is a BedVector corresponding to exons
typedef map<string, bedVector> TranscriptToExon;

//A vector of transcripts
typedef vector<string> TranscriptVector;

//Jump from a bin to a vector of transcripts
//The index for this vector is the UCSC bin
typedef map<int, TranscriptVector> BinToTranscripts;

//Jump from a chromosome and bin to transcript
typedef map<string, BinToTranscripts> ChrBinToTranscripts;

//Jump from a transcript ID to all the bins its exons fall in
typedef map<string, BIN> TranscriptToBin;

//Struct to hold each Transcript
struct Transcript {
    vector<BED> exons;
    vector<BED> junctions;
};

//Struct to hold each GTF line
class Gtf {
    public:
        string seqname; //Name of chromosome
        string source; //Program name
        string feature; //Feature type e.g gene
        CHRPOS start; //Feature start position, this type is from bedtools
        CHRPOS end; //Feature end position
        string score; //Feature score
        string strand; //Feature strand
        char frame; //Frame position of feature-start 0/1/2
        string attributes; //semi-colon delimited tag-value pairs
        bool is_exon; //Is this feature an exon?
        //Equals operator
        bool operator== (const Gtf &other) const {
            return seqname == other.seqname &&
                source == other.source &&
                feature == other.feature &&
                start == other.start &&
                end == other.end &&
                score == other.score &&
                strand == other.strand &&
                frame == other.frame &&
                attributes == other.attributes &&
                is_exon == other.is_exon;
        }
        //Insert operator
        friend ostream& operator<< (ostream &out, const Gtf& gtf1);
};

inline ostream& operator<< (ostream &out, const Gtf& gtf1) {
    out << gtf1.seqname << "\t";
    out << gtf1.source << "\t";
    out << gtf1.feature << "\t";
    out << gtf1.start << "\t";
    out << gtf1.end << "\t";
    out << gtf1.score << "\t";
    out << gtf1.strand << "\t";
    out << gtf1.frame << "\t";
    out << gtf1.attributes << "\t";
    out << gtf1.is_exon << endl;
    return out;
}

class GtfParser {
    private:
        //Name of the gtf file
        string gtffile_;
        //GTF filehandle
        ifstream gtf_fh_;
        //Are exons within transcripts sorted
        bool transcripts_sorted_;
        //Jump from transcript-id to {gene name, gene id}
        map<string, string*> transcript_to_gene_;
        //Store transcripts as a vector of exon BEDs
        //keyed by transcript_id
        map<string, Transcript> transcript_map_;
        //Bin for transcript
        TranscriptToBin transcript_to_bin_;
        //keyed by transcript_id
        ChrBinToTranscripts chrbin_to_transcripts_;
    public:
        //Constructor
        GtfParser()
            : transcripts_sorted_(false)
        {}
        //Constructor
        GtfParser(string gtf1)
            : gtffile_(gtf1)
            , transcripts_sorted_(false)
        {}
        //Copy constructor
        GtfParser(const GtfParser &gp1) {
            gtffile_ = gp1.gtffile_;
            transcripts_sorted_ = gp1.transcripts_sorted_;
            transcript_to_gene_ = gp1.transcript_to_gene_;
            transcript_map_ = gp1.transcript_map_;
            transcript_to_bin_ = gp1.transcript_to_bin_;
            chrbin_to_transcripts_ = gp1.chrbin_to_transcripts_;
        }
        //Parse an exon line into a gtf struct
        Gtf parse_exon_line(string line);
        //Parse the required field from attributes column
        string parse_attribute(vector<string> attributes1,
                           string field_name);
        //Close the gtf filehandle
        void close();
        //Assemble transcripts into a map
        //this is an associative container,
        //key is the transcript name, the value
        //is a vector of BEDs corresponding to exons.
        void create_transcript_map();
        //Add an exon to a transcript map
        void add_exon_to_transcript_map(Gtf gtf1);
        //Open the gtf file
        void open();
        //Set the gtf filename
        void set_gtffile(string filename);
        //Get the gtf filename
        string gtffile() { return gtffile_; }
        //Sort the exons within transcripts by start position
        void sort_exons_within_transcripts();
        //Construct junction information using exons
        void construct_junctions();
        //Annotate each transcript with its bin
        void annotate_transcript_with_bins();
        //Print out transcripts
        void print_transcripts();
        //Return vector of transcripts in a bin
        vector<string> transcripts_from_bin(string chr, BIN b1);
        //Return the bins that the exon-exon junctions
        //of a transcript fall in
        BIN bin_from_transcript(string transcript_id);
        //Return the exons corresponding to a transcript
        //The return value is a vector of BEDs
        const vector<BED> & get_exons_from_transcript(string transcript_id);
        //Get the gene ID using the trancript ID
        string* get_gene_from_transcript(string transcript_id);
        //Set the gene ID for a trancript ID
        void set_transcript_gene(string transcript_id, string gene_name, string gene_id);
        //Load all the necessary objects into memory
        void load();
        //Assignment operator
        GtfParser& operator= (const GtfParser& gtf1);
};

#endif
