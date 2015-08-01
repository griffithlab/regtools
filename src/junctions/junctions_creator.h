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

#ifndef JUNCTIONS_CREATOR_H
#define JUNCTIONS_CREATOR_H

#include <iostream>
#include "bedFile.h"
#include "sam.h"
//#include "faidx.h"
//#include "kstring.h"

using namespace std;

//Format of an junction
struct Junction : BED {
    int read_count;
};

//The class that deals with creating the junctions
class JunctionsCreator {
    private:
        //Alignment file
        string bam_;
        //Map to store the junctions
        //The key is "chr:start-end"
        //The value is an object of type Junction(see above)
        map<string, Junction> junctions;
        //Pull out the cigar string from the read
        int parse_read(bam_hdr_t *header, bam1_t *aln);
        //Parse junctions from the read and store in junction map
        int parse_cigar_into_junctions(string chr, int read_pos,
                                       uint32_t *cigar, int n_cigar);
        //Add a junction to the junctions map
        int add_junction(Junction j1);
    public:
        //Parse command-line options for this tool
        int parse_options(int argc, char *argv[]);
        //Print default usage
        int usage();
        //Identify exon-exon junctions
        int identify_junctions_from_BAM();
        //Print all the junctions
        void print_junctions();
};

#endif

