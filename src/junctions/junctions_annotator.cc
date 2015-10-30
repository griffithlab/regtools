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
#include "faidx.h"

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
        copy_stream(cout, out);
        return;
    }
    ofs_.open(output_file_.c_str());
    if(!ofs_.is_open())
        throw runtime_error("Unable to open " +
                            output_file_);
    copy_stream(ofs_, out);
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
    if(!line.fields.size() || line.fields[10].empty()) {
        stringstream position;
        position << line.chrom << ":" << line.start;
        throw runtime_error("Block sizes not found. Invalid line. " +
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
        adjust_junction_ends(line);
        return true;
    } else {
        return false;
    }
}

//Get the splice_site bases
bool JunctionsAnnotator::get_splice_site(AnnotatedJunction & line) {
    string position1 = line.chrom + ":" +
                      num_to_str(line.start + 1) + "-" + num_to_str(line.start + 2);
    string position2 = line.chrom + ":" +
                      num_to_str(line.end - 2) + "-" + num_to_str(line.end - 1);
    string seq1, seq2;
    try {
        seq1 = get_reference_sequence(position1);
        seq2 = get_reference_sequence(position2);
    } catch (const runtime_error& e) {
        throw e;
    }
    if(line.strand == "-") {
        seq1 = rev_comp(seq1);
        seq2 = rev_comp(seq2);
        line.splice_site = seq2 + "-" + seq1;
    } else {
        line.splice_site = seq1 + "-" + seq2;
    }
    return true;
}

//Extract gtf info
bool JunctionsAnnotator::read_gtf() {
    try {
        gtf_.create_transcript_map();
        gtf_.construct_junctions();
        gtf_.sort_exons_within_transcripts();
        gtf_.annotate_transcript_with_bins();
        //gtf_.print_transcripts();
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
    if(exons[0].start > junction.end &&
            exons[exons.size() - 1].end < junction.start)
        return false;
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
            if(!junction_start) {
                if(exons[i].end >= junction.start) {
                    junction_start = true;
                }
            }
            if(junction_start) {
                if(exons[i].start > junction.start &&
                        exons[i].end < junction.end) {
                    junction.exons_skipped.insert(exons[i].name);
                }
                if(exons[i].start > junction.start) {
                    junction.donors_skipped.insert(exons[i].start);
                }
                if(exons[i].end < junction.end) {
                    junction.acceptors_skipped.insert(exons[i].end);
                }
                if(exons[i].end == junction.start) {
                    junction.known_donor = true;
                }
                //TODO - check for last exon
                if(exons[i].start == junction.end) {
                    junction.known_acceptor = true;
                }
            }
        }
    }
    annotate_anchor(junction);
    return (junction.anchor != "N");
}

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
                        exons[i].end < junction.end) {
                    junction.exons_skipped.insert(exons[i].name);
                }
                if(exons[i].start > junction.start) {
                    junction.donors_skipped.insert(exons[i].start);
                }
                if(exons[i].end < junction.end) {
                    junction.acceptors_skipped.insert(exons[i].end);
                }
                if(exons[i].start == junction.end) {
                    junction.known_donor = true;
                }
                if(exons[i].end == junction.start) {
                    junction.known_acceptor = true;
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
                            + transcript_id);
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
        throw runtime_error("\nUnknown strand " + junction.strand);
    }
    cerr << "\n\tDonors skipped " << junction.donors_skipped.size();
    cerr << "\n\tExons skipped " << junction.exons_skipped.size();
    cerr << "\n\tAcceptors skipped " << junction.acceptors_skipped.size();
    cerr << "\n\tSplice site " << junction.anchor;
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
    if(s == NULL)
        throw runtime_error("Unable to extract FASTA sequence "
                             "for position " + position);
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
    while((c = getopt(argc, argv, "Eo:")) != -1) {
        switch(c) {
            case 'E':
                skip_single_exon_genes_ = false;
                break;
            case 'o':
                output_file_ = string(optarg);
                break;
            default:
                usage();
                throw runtime_error("\nError parsing inputs!(1)");
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
        throw runtime_error("\nError parsing inputs!(2)");
    }
    cerr << "\nReference: " << ref_;
    cerr << "\nGTF: " << gtf_.gtffile();
    cerr << "\nJunctions: " << junctions_.bedFile;
    if(skip_single_exon_genes_)
        cerr << "\nSkip single exon genes.";
    if(output_file_ != "NA")
        cerr << "\nOutput file: " << output_file_;
    return 0;
}

//Usage statement for this tool
int JunctionsAnnotator::usage() {
    cout << "\nUsage:\t\t" << "regtools junctions annotate [options] junctions.bed ref.fa annotations.gtf";
    cout << "\nOptions:\t" << "-E include single exon genes";
    cout << "\n\t\t" << "-o Output file";
    cout << "\n";
    return 0;
}
