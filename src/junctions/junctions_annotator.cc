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
#include <string>
#include "common.h"
#include "junctions_annotator.h"
#include "faidx.h"

using namespace std;

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
    //The junction start is start + block_size1
    //The junction end is end + block_size1
    if(line.fields[10].empty()) {
        cerr << "\nBlock sizes not found. Invalid line.";
        cerr << "\n\tChr " << line.chrom << "\tPos " << line.start;
    }
    string blocksize_field = line.fields[10];
    vector<int> block_sizes;
    Tokenize(blocksize_field, block_sizes, ',');
    line.start += block_sizes[0];
    line.end -= block_sizes[1] - 1;
}

//Get a single line from the junctions file
bool JunctionsAnnotator::get_single_junction(BED & line) {
    junctions_.GetNextBed(line);
    adjust_junction_ends(line);
    return true;
}

//Get the splice_site bases
bool JunctionsAnnotator::get_splice_site(AnnotatedJunction & line) {
    string position1 = line.chrom + ":" +
                      num_to_str(line.start + 1) + "-" + num_to_str(line.start + 2);
    string position2 = line.chrom + ":" +
                      num_to_str(line.end - 2) + "-" + num_to_str(line.end - 1);
    string seq1 = get_reference_sequence(position1);
    string seq2 = get_reference_sequence(position2);
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
    gtf_.create_transcript_map();
    gtf_.construct_junctions();
    gtf_.sort_exons_within_transcripts();
    gtf_.annotate_transcript_with_bins();
    //gtf_.print_transcripts();
    return true;
}

//Find overlap between transcript and junction on the negative strand,
//function returns true if this is a known junction in the transcript
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
        cerr << endl << "exon number " << i << "\t" << exons[i].start
             << "\t" << exons[i].end;
        if(exons[i].start > junction.end) {
            cerr << endl << "-1";
            //No need to look any further
            //the rest of the exons are outside the junction
            break;
        }
        //known junction
        if(exons[i].end == junction.start &&
                exons[i + 1].start == junction.end) {
            cerr << endl << "DA";
            junction.anchor = "DA";
            junction.known_acceptor = true;
            junction.known_donor = true;
            junction.known_junction = true;
            known_junction = true;
        }
        else {
            cerr << endl << "0";
            if(!junction_start) {
                if(exons[i].end >= junction.start) {
                    junction_start = true;
                    cerr << endl << "1";
                }
            }
            if(junction_start) {
                if(exons[i].start > junction.start &&
                        exons[i].end < junction.end) {
                    cerr << endl << "2";
                    junction.exons_skipped.insert(exons[i].name);
                }
                if(exons[i].start > junction.start) {
                    cerr << endl << "3";
                    junction.donors_skipped.insert(exons[i].start);
                }
                if(exons[i].end < junction.end) {
                    cerr << endl << "4";
                    junction.acceptors_skipped.insert(exons[i].end);
                }
                if(exons[i].end == junction.start) {
                    cerr << endl << "5";
                    junction.known_donor = true;
                }
                //TODO - check for last exon
                if(exons[i].start == junction.end) {
                    cerr << endl << "6";
                    junction.known_acceptor = true;
                }
            }
        }
    }
    annotate_anchor(junction);
    return known_junction;
}

//Find overlap between transcript and junction on the negative strand,
//function returns true if this is a known junction in the transcript
bool JunctionsAnnotator::overlap_ns(const vector<BED> & exons,
                                          AnnotatedJunction & junction) {
    cerr << endl << "in negative overlap";
    //skip single exon genes
    if(skip_single_exon_genes_ && exons.size() == 1) return false;
    bool junction_start = false;
    bool known_junction = false;
    //check if transcript overlaps with junction
    if(exons[0].end < junction.start ||
            exons[exons.size() - 1].start > junction.end) {
        cerr << endl << "transcript outside junction";
        return known_junction;
    }
    for(std::size_t i = 0; i < exons.size(); i++) {
        if(exons[i].end < junction.start) {
            cerr << endl << "-1";
            //No need to look any further
            //the rest of the exons are outside the junction
            break;
        }
        //Check if this is a known junction
        if(exons[i].start == junction.end &&
                exons[i + 1].end == junction.start) {
            cerr << endl << "DA";
            junction.anchor = "DA";
            junction.known_acceptor = true;
            junction.known_donor = true;
            junction.known_junction = true;
            known_junction = true;
        }
        else {
            cerr << endl << "0\t" << i;
            if(!junction_start) {
                if(exons[i].start <= junction.end) {
                    cerr << endl << "1";
                    junction_start = true;
                }
            }
            if(junction_start) {
                if(exons[i].start > junction.start &&
                        exons[i].end < junction.end) {
                    cerr << endl << "2";
                    junction.exons_skipped.insert(exons[i].name);
                }
                if(exons[i].start > junction.start) {
                    cerr << endl << "3";
                    junction.donors_skipped.insert(exons[i].start);
                }
                if(exons[i].end < junction.end) {
                    cerr << endl << "4";
                    junction.acceptors_skipped.insert(exons[i].end);
                }
                if(exons[i].start == junction.end) {
                    cerr << endl << "5";
                    junction.known_donor = true;
                }
                if(exons[i].end == junction.start) {
                    cerr << endl << "6";
                    junction.known_acceptor = true;
                }
            }
        }
    }
    annotate_anchor(junction);
    return known_junction;
}

//Annotate the anchor
void JunctionsAnnotator::annotate_anchor(AnnotatedJunction & junction) {
    //check if known junction
    if(!junction.known_junction) {
        if(junction.known_donor)
            if(junction.known_acceptor)
                junction.anchor = string("NDA");
            else
                junction.anchor = string("D");
        else if(junction.known_acceptor)
            junction.anchor = string("A");
        else
            junction.anchor = string("N");
    }
}

//Check for overlap between a transcript and junctions
//Check if the junction we saw is a known junction
//Calculate exons_skipped, donors_skipped, acceptors_skipped
void JunctionsAnnotator::check_for_overlap(string transcript_id, AnnotatedJunction & junction) {
    cerr << "\nTranscript id: " << transcript_id;
    const vector<BED> & exons =
        gtf_.get_exons_from_transcript(transcript_id);
    if(!exons.size()) {
        cerr << "Unexpected error. No exons for transcript "
             << transcript_id;
        exit(1);
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
        cerr << "\nUnknown strand " << junction.strand;
        exit(1);
    }
    cerr << "\n\tDonors skipped " << junction.donors_skipped.size();
    cerr << "\n\tExons skipped " << junction.exons_skipped.size();
    cerr << "\n\tAcceptors skipped " << junction.acceptors_skipped.size();
    cerr << "\n\tSplice site " << junction.anchor;
}

//Annotate with gtf
//Takes a single junction BED and annotates with GTF
void JunctionsAnnotator::annotate_junction_with_gtf(AnnotatedJunction & j1) {
    BIN junction_bin = getBin(j1.start, j1.end);
    //From BedTools
    BIN start_bin, end_bin;
    start_bin = (j1.start >> _binFirstShift);
    end_bin = ((j1.end-1) >> _binFirstShift);
    cerr << endl << "junction_bin " << junction_bin;
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
    static struct option long_options[] = {
        {"ref", required_argument, 0, 'r'},
        {"junctions", required_argument, 0, 'j'},
        {"gtf", required_argument, 0, 'g'},
    };
    int option_index = 0;
    int c = getopt_long(argc, argv, "r:g:j:",
                    long_options, &option_index);
    while(c != -1) {
        switch(c) {
            case 'r':
                ref_ = string(optarg);
                break;
            case 'g':
                gtf_.set_gtffile(string(optarg));
                break;
            case 'j':
                junctions_.bedFile = string(optarg);
                break;
            default:
                usage();
                exit(-1);
        }
        c = getopt_long(argc, argv, "r:g:j:",
                    long_options, &option_index);
    }
    //This could be made an option if need be
    skip_single_exon_genes_ = true;
    if(optind < argc || ref_.empty() || junctions_.bedFile.empty() || gtf_.gtffile().empty()) {
        usage();
        exit(-1);
    }
    cerr << "\nReference: " << ref_;
    cerr << "\nGTF: " << gtf_.gtffile();
    cerr << "\nJunctions: " << junctions_.bedFile;
    return 0;
}

//Usage statement for this tool
int JunctionsAnnotator::usage() {
    cout << "\nUsage:\t\t" << "regtools junctions annotate [options] -r ref.fa -j junctions.bed -g annotations.gtf";
    cout << "\nOptions:\t" << "-r ref.fa\tThe reference FASTA file.";
    cout << "\nOptions:\t" << "-j junctions.bed\tThe junctions to be annotated.";
    cout << "\nOptions:\t" << "-g annotations.gtf\tThe transcripts that we want to annotate the junctions with.";
    cout << "\n";
    return 0;
}


