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

#ifndef JUNCTIONS_EXTRACTOR_H
#define JUNCTIONS_EXTRACTOR_H

#include <iostream>
#include "bedFile.h"
#include "sam.h"

using namespace std;

//Format of an junction
struct Junction : BED {
    //Number of reads supporting the junction
    unsigned int read_count;
    //This is the start - max overhang
    CHRPOS thick_start;
    //This is the end + max overhang
    CHRPOS thick_end;
    //Has the junction been added to the map
    bool added;
    //Does the junction satisfy the min anchor requirement
    //Minimum anchor on either side can be supported by different
    //reads, only junctions anchored on both sides are reported.
    bool has_left_min_anchor;
    bool has_right_min_anchor;
    //Color for the BED line
    string color;
    //Number of blocks
    int nblocks;
    Junction() {
        start = 0;
        end = 0;
        thick_start = 0;
        thick_end = 0;
        read_count = 0;
        added = false;
        has_left_min_anchor = false;
        has_right_min_anchor = false;
        name = "NA";
        color = "255,0,0";
        nblocks = 2;
    }
    Junction(string chrom1, CHRPOS start1, CHRPOS end1,
             CHRPOS thick_start1, CHRPOS thick_end1,
             string strand1) {
        chrom = chrom1;
        start = start1;
        end = end1;
        thick_start = thick_start1;
        thick_end = thick_end1;
        name = "NA";
        read_count = 0;
        strand = strand1;
        added = false;
        has_left_min_anchor = false;
        has_right_min_anchor = false;
        color = "255,0,0";
        nblocks = 2;
    }
};

//Compare two junctions
//Return true if j1.start < j2.start
//If j1.start == j2.start, return true if j1.end < j2.end
static inline bool compare_junctions(const Junction &j1,
                       const Junction &j2) {
    //Different chromosome
    if(j1.chrom < j2.chrom)
        return true;
    if(j1.chrom > j2.chrom)
        return false;
    //Same chromosome
    if(j1.thick_start == j2.thick_start) {
        if(j1.thick_end < j2.thick_end)
            return true;
        else
            return false;
    }
    return j1.thick_start < j2.thick_start;
}

//The class that deals with creating the junctions
class JunctionsExtractor {
    private:
        //Alignment file
        string bam_;
        //Minimum anchor length for junctions
        //Junctions need atleast this many bp overlap
        // on both ends.
        uint32_t min_anchor_length_;
        //Minimum length of an intron, i.e min junction width
        uint32_t min_intron_length_;
        //Maximum length of an intron, i.e max junction width
        uint32_t max_intron_length_;
        //Map to store the junctions
        //The key is "chr:start-end"
        //The value is an object of type Junction(see above)
        map<string, Junction> junctions_;
        //Maintain a sorted list of junctions
        vector<Junction> junctions_vector_;
        //Are the junctions sorted
        bool junctions_sorted_;
        //File to write output to - optional, write to STDOUT by default
        string output_file_;
        //Region to identify junctions, in "chr:start-end" format
        string region_;
    public:
        //Default constructor
        JunctionsExtractor() {
            min_anchor_length_ = 8;
            min_intron_length_ = 70;
            max_intron_length_ = 500000;
            junctions_sorted_ = false;
            output_file_ = "NA";
            region_ = ".";
        };
        //Name the junction based on the number of junctions
        // in the map.
        string get_new_junction_name();
        //Parse command-line options for this tool
        int parse_options(int argc, char *argv[]);
        //Print default usage
        int usage(ostream& out = cerr);
        //Identify exon-exon junctions
        int identify_junctions_from_BAM();
        //Print one junction
        void print_one_junction(const Junction j1, ostream& out = cout);
        //Print all the junctions
        void print_all_junctions(ostream& out = cout);
        //Get a vector of all the junctions
        vector<Junction> get_all_junctions();
        //Get the BAM filename
        string get_bam();
        //Parse the alignment into the junctions map
        int parse_alignment_into_junctions(bam_hdr_t *header, bam1_t *aln);
        //Check if junction satisfies qc
        bool junction_qc(Junction &j1);
        //Sort all the junctions by their position
        void sort_junctions();
        //Create the junctions vector from the map
        void create_junctions_vector();
        //Pull out the cigar string from the read
        int parse_read(bam_hdr_t *header, bam1_t *aln);
        //Parse junctions from the read and store in junction map
        int parse_cigar_into_junctions(string chr, int read_pos,
                                       uint32_t *cigar, int n_cigar);
        //Add a junction to the junctions map
        int add_junction(Junction j1);
        //Get the strand from the XS aux tag
        void set_junction_strand(bam1_t *aln, Junction& j1);
};

#endif
