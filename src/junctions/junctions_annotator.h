/*  junctions_annotator.h -- class definitions for `junctions annotate`

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

#ifndef JUNCTIONS_ANNOTATOR_H_
#define JUNCTIONS_ANNOTATOR_H_

#include <iostream>
#include <iterator>
#include "bedFile.h"
#include "gtf_parser.h"

using namespace std;

struct AnnotatedJunction : BED {
    //set of transcripts that
    //the junction overlaps
    set<string> transcripts_overlap;
    //set of genes that
    //the junction overlaps
    set<string> genes_overlap;
    //set of exons that the junction
    //overlaps
    set<string> exons_skipped;
    //set of acceptor positions junction overlaps
    set<CHRPOS> acceptors_skipped;
    //set of donor positions junction overlaps
    set<CHRPOS> donors_skipped;
    //five prime reference seq
    string anchor_seq;
    //splice site annotation (D/DA/NA etc)
    string splice_site;
    void print_header() {
        cout << "chrom" << "\t" << "start" <<
                "\t" << "end" << "\t" << "name" <<
                "\t" << "score" << "\t" << "strand" <<
                "\t" << "anchor_seq" << "\t" << "acceptors_skipped" <<
                "\t" << "exons_skipped" << "\t" << "donors_skipped" <<
                "\t" << "splice_site" <<
                "\t" << "transcripts" << "\t" << "genes";
    }
    void print() {
        //remove trailing comma
        if(!splice_site.empty())
            splice_site.erase(splice_site.end() - 1);
        else
            splice_site = "unknown";
        cout << endl << chrom << "\t" << start <<
                "\t" << end << "\t" << name <<
                "\t" << score << "\t" << strand <<
                "\t" << anchor_seq << "\t" << acceptors_skipped.size()<<
                "\t" << exons_skipped.size() << "\t" << donors_skipped.size() <<
                "\t" << splice_site;
        if(transcripts_overlap.size()) {
            cout << "\t";
            for(set<string>::iterator it = transcripts_overlap.begin(); it != transcripts_overlap.end(); ++it) {
                if(it != transcripts_overlap.begin())
                    cout << ",";
                cout << *it;
            }
        }
        if(genes_overlap.size()) {
            cout << "\t";
            for(set<string>::iterator it = genes_overlap.begin(); it != genes_overlap.end(); ++it) {
                if(it != genes_overlap.begin())
                    cout << ",";
                cout << *it;
            }
        }
    }
    void reset() {
        anchor_seq = "";
        splice_site = "";
        exons_skipped.clear();
        acceptors_skipped.clear();
        donors_skipped.clear();
        transcripts_overlap.clear();
        genes_overlap.clear();
    }
};

class JunctionsAnnotator {
    private:
        //Junctions file to be annotated
        BedFile junctions_;
        //Reference FASTA file
        string ref_;
        //skip single exon genes
        bool skip_single_exon_genes_;
        //GTF file object
        GtfParser gtf_;
        //Check for overlap between a transcript and junctions
        //See if the junction we saw is a known junction
        void check_for_overlap(string transcript_id,
                               AnnotatedJunction & junction);
        //Find overlap for transcripts on the positive strand
        bool positive_overlap(const vector<BED> & exons,
                              AnnotatedJunction & j1);
        //Find overlap for transcripts on the positive strand
        bool negative_overlap(const vector<BED> & exons,
                              AnnotatedJunction & j1);
    public:
        //Default constructor
        JunctionsAnnotator() {};
        //Destructor
        ~JunctionsAnnotator() {};
        //Get the GTF file
        string gtf_file();
        //Parse command-line options for this tool
        int parse_options(int argc, char *argv[]);
        //Print default usage
        int usage();
        //Get the reference bases at a position
        string get_reference_sequence(string position);
        //Get a single line from the junctions file
        bool get_single_junction(BED & line);
        //Get the anchor bases
        bool get_anchor_seq(AnnotatedJunction & line);
        //Open junctions file
        void open_junctions();
        //Close junctions file
        void close_junctions();
        //Extract gtf info
        bool read_gtf();
        //Annotate with gtf
        void annotate_junction_with_gtf(AnnotatedJunction & j1);
        //Adjust the start and end of the junction
        void adjust_junction_ends(BED & line);
};

#endif

