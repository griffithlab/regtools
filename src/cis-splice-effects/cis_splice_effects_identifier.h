/*  cis_splice_effects_identifier.h -- 'cis-splice-effects identify'

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
#ifndef CIS_SPLICE_EFFECTS_IDENTIFIER_
#define CIS_SPLICE_EFFECTS_IDENTIFIER_

#include "junctions_annotator.h"
#include "junctions_extractor.h"
#include "variants_annotator.h"

//Workhorse for "cis-splice-effects identify"
class CisSpliceEffectsIdentifier {
    private:
        //Object that annotates the variants
        VariantsAnnotator va;
        //Object that extracts the junctions
        JunctionsExtractor je;
        //Object that annotates the junctions
        JunctionsAnnotator ja;
        //VCF file with variants
        string vcf_;
        //RNAseq alignments
        string bam_;
        //Reference sequence FASTA
        string ref_;
        //GTF file with annotations
        string gtf_;
        //File to write output to
        string output_file_;
        //Aberrant splice junctions in BED12 format
        string output_junctions_bed_;
        //File to write barcodes to
        string output_barcodes_file_;
        //File to write output to
        string annotated_variant_file_;
        //Flag to indicate whether to write output vcf
        bool write_annotated_variants_;
        //Window size to look in
        //Looks at variant.pos +/- window_size
        uint32_t window_size_;
        //output stream to output annotated junctions file
        ofstream ofs_;
        //output stream to output BED12 junctions file
        ofstream ofs_junctions_bed_;
        //output stream to output barcodes file
        ofstream ofs_barcode_counts_;
        //Unique set of junctions near splicing variants
        set<Junction> unique_junctions_;
        //Store the variants that mark a junction
        map<Junction, set<AnnotatedVariant> > junction_to_variant_;
        //Minimum distance of a variant from
        //edge of an exon(intronic) to be considered
        //a splicing variant
        uint32_t intronic_min_distance_;
        //Minimum distance of a variant from
        //edge of an exon(exonic) to be considered
        //a splicing variant
        uint32_t exonic_min_distance_;
        //Flag set by the -I option
        bool all_intronic_space_;
        //Flag set by the -E option
        bool all_exonic_space_;
        //Option to skip single exon genes
        bool skip_single_exon_genes_;
        //strandness of data
        int strandness_;
        //tag used in BAM to denote strand, default "XS"
        string strand_tag_;
        //Minimum anchor length for junctions
        //Junctions need at least this many bp overlap on both ends.
        uint32_t min_anchor_length_;
        //Reads need at least this many bp overlap on both ends to support junction.
        uint32_t min_read_anchor_length_;
        //Minimum length of an intron, i.e min junction width
        uint32_t min_intron_length_;
        //Maximum length of an intron, i.e max junction width
        uint32_t max_intron_length_;
        //filter reads containing any of these flags
        uint16_t filter_flags_;
        // filter reads not containing all of these flags
        uint16_t require_flags_;
        // filter reads below the minimum mapping quality
        uint8_t min_map_qual_;
        //whether to override strand of extracted junctions using intron-motif method
        bool override_strand_with_canonical_intron_motif_;
    public:
        //Constructor
        CisSpliceEffectsIdentifier() : vcf_("NA"), output_file_("NA"),
                                       output_junctions_bed_("NA"),
                                       output_barcodes_file_("NA"),
                                       annotated_variant_file_("NA"),
                                       write_annotated_variants_(false),
                                       window_size_(0),
                                       intronic_min_distance_(2),
                                       exonic_min_distance_(3),
                                       all_intronic_space_(false),
                                       all_exonic_space_(false),
                                       skip_single_exon_genes_(true),
                                       strandness_(-1),
                                       strand_tag_("XS"),
                                       min_anchor_length_(8),
                                       min_read_anchor_length_(0),
                                       min_intron_length_(70),
                                       max_intron_length_(500000),
                                       filter_flags_(0),
                                       require_flags_(0),
                                       min_map_qual_(0),
                                       override_strand_with_canonical_intron_motif_(false) {}
        //Destructor
        ~CisSpliceEffectsIdentifier() {
            if(ofs_.is_open()) {
                ofs_.close();
            }
            if(ofs_junctions_bed_.is_open()) {
                ofs_junctions_bed_.close();
            }
            if(ofs_barcode_counts_.is_open()) {
                ofs_barcode_counts_.close();
            }
        }
        //Parse command line arguments
        void parse_options(int argc, char* argv[]);
        //Check if files exist
        void file_qc();
        //Identify cis splicing effects
        void identify();
        //Usage for this tool
        void usage(ostream &out);
        //Set ofstream object to appropriate value
        //This could write to output_file_ or to std::cout
        void set_ostream();
        //Close ofstream object
        void close_ostream();
        //Get the window size
        uint32_t window_size() { return window_size_; }
        //Get the annotated variants file(VCF)
        string annotated_variant_file() { return annotated_variant_file_; }
        //Get the barcodes file
        string output_barcodes_file() { return output_barcodes_file_; }
        //Get the file with splice junctions(BED)
        string output_file() { return output_file_; }
        //Get the Input VCF
        string vcf() { return vcf_; }
        //Call the junctions annotator
        void annotate_junctions(const GtfParser& gtf_p1);
        //Get junction name given an index
        string get_junction_name(int i);
};

#endif
