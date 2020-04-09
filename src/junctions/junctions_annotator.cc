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
#include <stdexcept>
#include <string>
#include "common.h"
#include "junctions_annotator.h"
#include "htslib/faidx.h"

using namespace std;

//Return stream to write output to
void JunctionsAnnotator::close_ofstream() {
    if(ofs_.is_open())
        ofs_.close();
}

//Return stream to write output to
//If output file is not empty, attempt to open
//If output file is empty, set to cout
void JunctionsAnnotator::set_ofstream_object(ofstream &out) {
    if(output_file_ == "NA") {
        common::copy_stream(cout, out);
        return;
    }
    ofs_.open(output_file_.c_str());
    if(!ofs_.is_open())
        throw runtime_error("Unable to open " +
                            output_file_);
    common::copy_stream(ofs_, out);
}

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
    //The junction start is thick_start + block_size1
    //The junction end is thick_end - block_size2 + 1
    if(line.fields.size() != 12  || line.fields[10].empty()) {
        stringstream position;
        position << line.chrom << ":" << line.start;
        throw runtime_error("BED line not in BED12 format. start: " +
                            position.str());
    }
    string blocksize_field = line.fields[10];
    vector<int> block_sizes;
    Tokenize(blocksize_field, block_sizes, ',');
    line.start += block_sizes[0];
    line.end -= block_sizes[1] - 1;
}

//Get a single line from the junctions file
bool JunctionsAnnotator::get_single_junction(BED & line) {
    junctions_._status = BED_INVALID;
    if(junctions_.GetNextBed(line) && junctions_._status == BED_VALID) {
        return true;
    } else {
        return false;
    }
}

//Get the splice_site bases
void JunctionsAnnotator::get_splice_site(AnnotatedJunction & line) {
    string position1 = line.chrom + ":" +
                      common::num_to_str(line.start + 1) + "-" + common::num_to_str(line.start + 2);
    string position2 = line.chrom + ":" +
                      common::num_to_str(line.end - 2) + "-" + common::num_to_str(line.end - 1);
    string seq1, seq2;
    try {
        seq1 = get_reference_sequence(position1);
        seq2 = get_reference_sequence(position2);
    } catch (const runtime_error& e) {
        throw e;
    }
    if(line.strand == "-") {
        seq1 = common::rev_comp(seq1);
        seq2 = common::rev_comp(seq2);
        line.splice_site = seq2 + "-" + seq1;
    } else {
        line.splice_site = seq1 + "-" + seq2;
    }
    return;
}

//Extract gtf info
bool JunctionsAnnotator::load_gtf() {
    try {
        gtf_.load();
    } catch (runtime_error e) {
        throw e;
    }
    return true;
}

//Find overlap between transcript and junction on the negative strand,
//function returns true if either the acceptor or the donor is known
bool JunctionsAnnotator::overlap_ps(const vector<BED>& exons,
                                          AnnotatedJunction & junction) {
    //skip single exon genes
    if(skip_single_exon_genes_ && exons.size() == 1) return false;
    bool junction_start = false;
    bool known_junction = false;
    //check if transcript overlaps with junction
    if(exons[0].start > junction.end ||
            exons[exons.size() - 1].end < junction.start)
        return known_junction;
    for(std::size_t i = 0; i < exons.size(); i++) {
        if(exons[i].start > junction.end) {
            //No need to look any further
            //the rest of the exons are outside the junction
            break;
        }
        //known junction
        if(exons[i].end == junction.start &&
                exons[i + 1].start == junction.end) {
            junction.known_acceptor = true;
            junction.known_donor = true;
            junction.known_junction = true;
            known_junction = true;
        }
        else {
            // YY NOTE: if we have not yet reached the junction on this transcript,
            if(!junction_start) {
                // then check if we have with this exon (i.e. the 
                // exon end is past or at the junction start)
                if(exons[i].end >= junction.start) {
                    junction_start = true;
                }
            }
            if(junction_start) {
                // YY NOTE: if the exon lies completely within the junction,
                //              count it as skipped
                if(exons[i].start > junction.start &&
                        exons[i].end < junction.end &&
                        i > 0 && i < (exons.size()-1)) {
                    string exon_coords = common::num_to_str(exons[i].start) + "-" + common::num_to_str(exons[i].end);
                    junction.exons_skipped.insert(exon_coords);
                }
                // YY NOTE: if the exon ends after the junction starts
                //              (and junction_start == true), then 
                //              count the donor as skipped
                        // YY NOTE: maybe it should be junction.end-1? since
                        //              if the "donor" and "acceptor" are already
                        //              contiguous in the DNA it would get counted
                        //              as "skipped" 
                if(exons[i].end > junction.start &&
                        exons[i].end < junction.end && 
                        i < (exons.size()-1)) {
                    junction.donors_skipped.insert(exons[i].end);
                }
                // YY NOTE: if the exon starts before the junction ends
                //              (and junction_start == true), then
                //              count the acceptor as skipped
                if(exons[i].start < junction.end &&
                        exons[i].start > junction.start &&
                        i > 0) {
                    junction.acceptors_skipped.insert(exons[i].start);
                }
                if(exons[i].end == junction.start) {
                    junction.known_donor = true;
                }
                if(exons[i].start == junction.end) {
                    junction.known_acceptor = true;
                }
            }
        }
    }
    annotate_anchor(junction);
    return (junction.anchor != "N");
}

//YY NOTES
/*

So based on my understanding of the bin search, we should be looking at all the 
possible transcripts that might overlap with a junction. That is to say, it 
would seem that any discrepancy in acceptors/exons/donors skipped would come 
from a problem with the actual counting.

The overlap function will be run on each transcript, and you can see from the
if block how it judges an acceptor/exon/donor as skipped. That logic appears
correct to me (see notes above). Doesn't immediately seem like a problem
with the actual recording of acceptors/exons/donors skipped (but probably
warrants a closer look).

Then the final place the count could get messed up is in the actual
printing to the tsv. The way this works is the skipped elements are 
held in sets. Acceptors/donors skipped are held in sets based on their 
coordinates while exons skipped are held in sets based on their names.
The number of elements skipped is simply the size of the set, so we only 
count unique elements. Also don't see a problem immediately with how this works.

*/

//Find overlap between transcript and junction on the negative strand,
//function returns true if either the acceptor or the donor is known
bool JunctionsAnnotator::overlap_ns(const vector<BED> & exons,
                                          AnnotatedJunction & junction) {
    //skip single exon genes
    if(skip_single_exon_genes_ && exons.size() == 1) return false;
    bool junction_start = false;
    bool known_junction = false;
    //check if transcript overlaps with junction
    if(exons[0].end < junction.start ||
            exons[exons.size() - 1].start > junction.end) {
        return known_junction;
    }
    for(std::size_t i = 0; i < exons.size(); i++) {
        if(exons[i].end < junction.start) {
            //No need to look any further
            //the rest of the exons are outside the junction
            break;
        }
        //Check if this is a known junction
        if(exons[i].start == junction.end &&
                exons[i + 1].end == junction.start) {
            junction.known_acceptor = true;
            junction.known_donor = true;
            junction.known_junction = true;
            known_junction = true;
        }
        else {
            if(!junction_start) {
                if(exons[i].start <= junction.end) {
                    junction_start = true;
                }
            }
            if(junction_start) {
                if(exons[i].start > junction.start &&
                        exons[i].end < junction.end &&
                        i > 0 && i < (exons.size()-1)) {
                    string exon_coords = common::num_to_str(exons[i].start) + "-" + common::num_to_str(exons[i].end);
                    junction.exons_skipped.insert(exon_coords);
                }
                // YY NOTE: if the exon ends after the junction starts
                //              (and junction_start == true), then 
                //              count the donor as skipped
                if(exons[i].end > junction.start &&
                        exons[i].end < junction.end && 
                        i < (exons.size()-1)) {
                    junction.acceptors_skipped.insert(exons[i].end);
                }
                // YY NOTE: if the exon starts before the junction ends
                //              (and junction_start == true), then
                //              count the acceptor as skipped
                if(exons[i].start < junction.end &&
                        exons[i].start > junction.start) {
                    junction.donors_skipped.insert(exons[i].start);
                }
                if(exons[i].end == junction.start) {
                    junction.known_acceptor = true;
                }
                if(exons[i].start == junction.end) {
                    junction.known_donor = true;
                }
            }
        }
    }
    annotate_anchor(junction);
    return (junction.anchor != "N");
}

//Annotate the anchor i.e is this a known/novel donor-acceptor pair
void JunctionsAnnotator::annotate_anchor(AnnotatedJunction & junction) {
    junction.anchor = string("N");
    if(junction.known_junction) {
        junction.anchor = "DA";
    } else {
        if(junction.known_donor)
            if(junction.known_acceptor)
                junction.anchor = string("NDA");
            else
                junction.anchor = string("D");
        else if(junction.known_acceptor)
            junction.anchor = string("A");
    }
}

//Check for overlap between a transcript and junctions
//Check if the junction we saw is a known junction
//Calculate exons_skipped, donors_skipped, acceptors_skipped
void JunctionsAnnotator::check_for_overlap(string transcript_id, AnnotatedJunction & junction) {
    const vector<BED> & exons =
        gtf_.get_exons_from_transcript(transcript_id);
    if(!exons.size()) {
        throw runtime_error("Unexpected error. No exons for transcript "
                            + transcript_id + "\n\n");
    }
    string transcript_strand = exons[0].strand;
    //Make sure the strands of the junction and transcript match
    if(junction.strand != transcript_strand)
        return;
    //Remember exons are sorted from exon1 to last exon
    if(junction.strand == "+") {
        if(overlap_ps(exons, junction)) {
            junction.transcripts_overlap.insert(transcript_id);
            junction.genes_overlap.insert(
                    gtf_.get_gene_from_transcript(transcript_id));
        }
    } else if(junction.strand == "-") {
        if(overlap_ns(exons, junction)) {
            junction.transcripts_overlap.insert(transcript_id);
            junction.genes_overlap.insert(
                    gtf_.get_gene_from_transcript(transcript_id));
        }
    } else {
        throw runtime_error("Unknown strand " + junction.strand + "\n\n");
    }
}

//Annotate with gtf
//Takes a single junction BED and annotates with GTF
void JunctionsAnnotator::annotate_junction_with_gtf(AnnotatedJunction & j1) {
    //From BedTools
    BIN start_bin, end_bin;
    start_bin = (j1.start >> _binFirstShift);
    end_bin = ((j1.end-1) >> _binFirstShift);
    //We loop through each UCSC BIN level for feature A's chrom.
    //For each BIN, we loop through each B feature and add it to
    // hits if it meets all of the user's requests.
    for (BINLEVEL i = 0; i < _binLevels; ++i) {
        BIN offset = _binOffsetsExtended[i];
        for (BIN b = (start_bin + offset); b <= (end_bin + offset); ++b) {
            vector<string> transcripts = gtf_.transcripts_from_bin(j1.chrom, b);
            if(transcripts.size())
                for(std::size_t i = 0; i < transcripts.size(); i++)
                    check_for_overlap(transcripts[i], j1);
        }
        start_bin >>= _binNextShift;
        end_bin >>= _binNextShift;
    }
}

//Get the reference sequence at a particular coordinate
string JunctionsAnnotator::get_reference_sequence(string position) {
    int len;
    faidx_t *fai = fai_load(ref_.c_str());
    char *s = fai_fetch(fai, position.c_str(), &len);
    cerr << "position = " << position << endl;
    if(s == NULL)
        throw runtime_error("Unable to extract FASTA sequence "
                             "for position " + position + "\n\n");
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
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    while((c = getopt(argc, argv, "So:h")) != -1) {
        switch(c) {
            case 'S':
                skip_single_exon_genes_ = false;
                break;
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'h':
                usage(help_ss);
                throw common::cmdline_help_exception(help_ss.str());
            default:
                usage();
                throw runtime_error("Error parsing inputs!(1)\n\n");
        }
    }
    if(argc - optind >= 3) {
        junctions_.bedFile = string(argv[optind++]);
        ref_ = string(argv[optind++]);
        gtf_.set_gtffile(string(argv[optind++]));
    }
    if(optind < argc ||
       ref_ == "NA" ||
       junctions_.bedFile.empty() ||
       gtf_.gtffile().empty()) {
        usage();
        throw runtime_error("Error parsing inputs!(2)\n\n");
    }
    cerr << "Reference: " << ref_ << endl;
    cerr << "GTF: " << gtf_.gtffile() << endl;
    cerr << "Junctions: " << junctions_.bedFile << endl;
    if(skip_single_exon_genes_)
        cerr << "Skipping single exon genes." << endl;
    if(output_file_ != "NA")
        cerr << "Output file: " << output_file_ << endl;
    cerr << endl;
    return 0;
}

//Usage statement for this tool
int JunctionsAnnotator::usage(ostream& out) {
    out << "Usage:\t\t" << "regtools junctions annotate [options] junctions.bed ref.fa annotations.gtf" << endl;
    out << "Options:\t" << "-S include single exon genes" << endl;
    out << "\t\t" << "-o FILE\tThe file to write output to. [STDOUT]" << endl;
    out << endl;
    return 0;
}
