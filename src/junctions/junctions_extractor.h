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

#include <algorithm>
#include <iomanip>
#include <iostream>
#include "bedFile.h"
#include "htslib/sam.h"
#include <unordered_map>
#include <unordered_set>

using namespace std;

//Format of an junction
struct Junction : BED {
    //Number of reads supporting the junction
    unsigned int read_count;
    //Reads supporting the junction
    unordered_set<string> reads;
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
    // single cell, key is barcode, val is count of that barcode
    unordered_map<string, int> barcodes;
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
    //Print junction
    void print(ostream& out) const {
        out << chrom <<
            "\t" << thick_start << "\t" << thick_end <<
            "\t" << name << "\t" << read_count << "\t" << strand <<
            "\t" << thick_start << "\t" << thick_end <<
            "\t" << color << "\t" << nblocks <<
            "\t" << start - thick_start << "," << thick_end - end <<
            "\t" << "0," << end - thick_start << endl;
    }
    void print_barcodes(ostream& out) const {
        // map - I'm not sure why .begin() is giving me a const_iterator (suggesting barcodes is const-qualified)
        //  but see map<string, Junction> :: iterator it = junctions_.begin() in JunctionsExtractor::create_junctions_vector()
        //  ... idk maybe it has something to do with the fact that junctions_ is a private member of a class whereas barcodes is just a struct field?
        out << barcodes.size() << "\t";
        for (unordered_map<string, int>::const_iterator it = barcodes.begin(); it != barcodes.end(); it++){
             if (it != barcodes.begin()){
                 out << ",";
             }
             out << it->first << ":" << it->second;
        }
        out << endl;
    }
};

//Compare two junctions
//Return true if j1.start < j2.start
//If j1.start == j2.start, return true if j1.end < j2.end
static inline bool compare_junctions(const Junction &j1,
                       const Junction &j2) {
    //Different chromosome
    if(j1.chrom < j2.chrom){
        return true;
    }
    if(j1.chrom > j2.chrom){
        return false;
    }
    //Same chromosome
    if(j1.thick_start < j2.thick_start) {
        return true;
    }
    if(j1.thick_start > j2.thick_start) {
        return false;
    }
    if(j1.thick_end < j2.thick_end) {
        return true;
    }
    if(j1.thick_end > j2.thick_end) {
        return false;
    }
    return j1.name < j2.name;
}

//Sort a vector of junctions
template <class CollectionType>
inline void sort_junctions(CollectionType &junctions) {
    sort(junctions.begin(), junctions.end(), compare_junctions);
}

//The class that deals with creating the junctions
class JunctionsExtractor {
    private:
        //Alignment file
        string bam_;
        //Reference FASTA file
        string ref_;
        //Minimum anchor length for junctions
        //Junctions need at least this many bp overlap on both ends.
        uint32_t min_anchor_length_;
        //Reads need at least this many bp overlap to support a junction
        uint32_t min_read_anchor_length_;
        //Minimum length of an intron, i.e min junction width
        uint32_t min_intron_length_;
        //Maximum length of an intron, i.e max junction width
        uint32_t max_intron_length_;
        //Map to store the junctions
        //The key is "chr:start-end:strand"
        //The value is an object of type Junction(see above)
        map<string, Junction> junctions_;
        //Maintain a sorted list of junctions
        vector<Junction> junctions_vector_;
        //Are the junctions sorted
        bool junctions_sorted_;
        //File to write output to - optional, write to STDOUT by default
        string output_file_;
        //File to write barcodes to
        string output_barcodes_file_;
        //Region to identify junctions, in "chr:start-end" format
        string region_;
        //strandness of data; 0 = unstranded, 1 = RF, 2 = FR, 3 = intron-motif
        int strandness_;
        //tag used in BAM to denote strand, default "XS"
        string strand_tag_;
        //tag used in BAM to denote single cell barcode
        string barcode_tag_;
        //filter reads containing any of these flags
        uint16_t filter_flags_;
        // filter reads not containing all of these flags
        uint16_t require_flags_;
        // filter reads below the minimum mapping quality
        uint8_t min_map_qual_;
    public:
        //Default constructor
        JunctionsExtractor() {
            min_anchor_length_ = 8;
            min_read_anchor_length_ = 0;
            min_intron_length_ = 70;
            max_intron_length_ = 500000;
            filter_flags_ = 0;
            require_flags_ = 0;
            min_map_qual_ = 0;
            junctions_sorted_ = false;
            strandness_ = -1;
            strand_tag_ = "XS";
            barcode_tag_ = "CB";
            bam_ = "NA";
            output_file_ = "NA";
            output_barcodes_file_ = "NA";
            region_ = ".";
            ref_ = "NA";
        }
        JunctionsExtractor(
                string bam1, 
                string region1, 
                int strandness1, 
                string strand_tag1, 
                uint32_t min_anchor_length1, 
                uint32_t min_read_anchor_length1, 
                uint32_t min_intron_length1, 
                uint32_t max_intron_length1, 
                uint16_t filter_flags, 
                uint16_t require_flags, 
                uint8_t min_map_qual, 
                string ref1) : 
            bam_(bam1), 
            region_(region1), 
            strandness_(strandness1), 
            strand_tag_(strand_tag1), 
            min_anchor_length_(min_anchor_length1), 
            min_read_anchor_length_(min_read_anchor_length1), 
            min_intron_length_(min_intron_length1), 
            max_intron_length_(max_intron_length1), 
            filter_flags_(filter_flags), 
            require_flags_(require_flags), 
            min_map_qual_(min_map_qual), 
            ref_(ref1) {
            junctions_sorted_ = false;
            output_file_ = "NA";
            output_barcodes_file_ = "NA";
            barcode_tag_ = "CB";
        }
        //Name the junction based on the number of junctions
        // in the map.
        string get_new_junction_name();
        //Parse command-line options for this tool
        int parse_options(int argc, char *argv[]);
        //Print default usage
        int usage(ostream& out = cerr);
        //Identify exon-exon junctions
        int identify_junctions_from_BAM();
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
        //Create the junctions vector from the map
        void create_junctions_vector();
        //Pull out the cigar string from the read
        int parse_read(bam_hdr_t *header, bam1_t *aln);
        //Returns whether alignment should be filtered from junction analysis
        bool filter_alignment(bam_hdr_t *header, bam1_t *aln);
        //Parse junctions from the read and store in junction map
        int parse_cigar_into_junctions(string chr, int read_pos,
                                       uint32_t *cigar, int n_cigar);
        //Add a junction to the junctions map
        int add_junction(Junction j1);
        //Get the strand from the XS aux tag
        void set_junction_strand_XS(bam1_t *aln, Junction& j1);
        //Get the strand from bitwise flag
        void set_junction_strand_flag(bam1_t *aln, Junction& j1);
        //Infer strand from canonical-motifs
        void set_junction_strand_intron_motif(string intron_motif, Junction& j1);
        //Get the strand
        void set_junction_strand(bam1_t *aln, Junction& j1, string intron_motif);
        //Get the barcode
        void set_junction_barcode(bam1_t *aln, Junction& j1);
        //Get the reference bases at a position
        string get_reference_sequence(string position);
        //Get the anchor bases
        string get_splice_site(Junction & line);
};

#endif
