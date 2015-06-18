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
#include "bedFile.h"
#include "gtf_parser.h"

using namespace std;

struct AnnotatedJunctions : BED {
    int n_exons_skipped;
    int n_acceptors_skipped;
    int n_donors_skipped;
    string five_prime_seq;
    string three_prime_seq;
};

class JunctionsAnnotator {
    private:
        //Junctions file to be annotated
        BedFile junctions_;
        //Reference FASTA file
        string ref_;
        //GTF file object
        GtfParser gtf_;
        //Get the anchor bases
        bool get_anchor_bases(const BED & line);
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
        //Annotate a single line
        bool annotate_single_line(BED & line);
        //Get a single line from the junctions file
        void get_single_junction(BED & line);
        //Open junctions file
        void open_junctions();
        //Close junctions file
        void close_junctions();
        //Extract gtf info
        bool read_gtf();
        //Annotate with gtf
        void annotate_junction_with_gtf(BED j1);
        //Adjust the start and end of the junction
        void adjust_junction_ends(BED & line);
};

#endif

