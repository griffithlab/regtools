/*  variants_annotator.cc -- `variants annotate`

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

#include "bedFile.h"
#include "common.h"
#include "hts.h"
#include "variants_annotator.h"
#include <cstdlib>
#include <stdexcept>

//Usage statement for this tool
int VariantsAnnotator::usage(ostream& out) {
    out << "\nUsage:\t\t" << "regtools variants annotate [options] variants.vcf annotations.gtf annotated_output.vcf";
    out << "\n\t\t" << "-e INT\tMinimum distance from the start/end of an exon "
                       "\n\t\t\tto annotate a variant as relevant to splicing, the variant "
                       "\n\t\t\tis in exonic space, i.e a coding variant. [3]";
    out << "\n\t\t" << "-i INT\tMinimum distance from the start/end of an exon "
                       "\n\t\t\tto annotate a variant as relevant to splicing, the variant "
                       "\n\t\t\tis in intronic space. [2]";
    out << "\n\t\t" << "-o\tFile to write output to. [STDOUT]";
    out << "\n\t\t" << "-S\tDon't skip single exon transcripts.";
    out << "\n";
    return 0;
}

//Parse command line options
int VariantsAnnotator::parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int16_t c;
    stringstream help_ss;
    while((c = getopt(argc, argv, "e:hi:o:S")) != -1) {
        switch(c) {
            case 'i':
                intronic_min_distance_ = atoi(optarg);
                break;
            case 'e':
                exonic_min_distance_ = atoi(optarg);
                break;
            case 'o':
                vcf_out_ = string(optarg);
                break;
            case 'S':
                skip_single_exon_genes_ = false;
                break;
            case 'h':
                usage(help_ss);
                throw common::cmdline_help_exception(help_ss.str());
            default:
                usage(std::cout);
                throw runtime_error("\nError parsing inputs!(1)\n");
        }
    }
    if(argc - optind >= 2) {
        vcf_ = string(argv[optind++]);
        gtffile_ = string(argv[optind++]);
        gtf_.set_gtffile(gtffile_);
    }
    if(optind < argc ||
       vcf_ == "NA" ||
       gtffile_ == "NA") {
        usage(std::cout);
        throw runtime_error("\nError parsing inputs!(2)\n");
    }
    cerr << "\nVariant file: " << vcf_;
    cerr << "\nGTF file: " << gtffile_;
    cerr << "\nOutput vcf file: " << vcf_out_;
    cerr << "\nIntronic min distance: " << intronic_min_distance_;
    cerr << "\nExonic min distance: " << exonic_min_distance_;
    if(!skip_single_exon_genes_)
        cerr << "\nNot skipping single exon genes.";
    if(vcf_out_ != "NA")
        cerr << "\nOutput file: " << vcf_out_;
    cerr << endl;
    return 0;
}

//Read gtf info into gtf_
void VariantsAnnotator::load_gtf() {
    gtf_.load();
}

//Open input VCF file
void VariantsAnnotator::open_vcf_in() {
    vcf_fh_in_ = bcf_open(vcf_.c_str(), "r");
    if(vcf_fh_in_ == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    vcf_header_in_ = bcf_hdr_read(vcf_fh_in_);
    if(vcf_header_in_ == NULL) {
        throw std::runtime_error("Unable to read header.");
    }
}

//Open output VCF file
void VariantsAnnotator::open_vcf_out() {
    vcf_fh_out_ =  hts_open(vcf_out_ == "NA" ? "-" : vcf_out_.c_str(),
                            "w");
    if(vcf_fh_out_ == NULL) {
        throw runtime_error("Unable to open output VCF file");
    }
    vcf_header_out_ = bcf_hdr_dup(vcf_header_in_);
    bcf_hdr_append(vcf_header_out_,
                   "##INFO=<ID=genes,Number=1,Type=String,"
                   "Description=\"The Variant falls in the splice "
                   "region of these genes\"");
    bcf_hdr_append(vcf_header_out_,
                   "##INFO=<ID=transcripts,Number=1,Type=String,"
                   "Description=\"The Variant falls in the splice "
                   "region of these transcripts\"");
    bcf_hdr_append(vcf_header_out_,
                   "##INFO=<ID=distances,Number=1,Type=String,"
                   "Description=\"Vector of Min(Distance from start/end of exon in the transcript.)\"");
    bcf_hdr_append(vcf_header_out_,
                   "##INFO=<ID=annotations,Number=1,Type=String,"
                   "Description=\"Does the variant fall in exonic/intronic splicing "
                   "related space in the transcript.\"");
    bcf_hdr_sync(vcf_header_out_);
    bcf_hdr_write(vcf_fh_out_, vcf_header_out_);
}

//Free relevant pointers
void VariantsAnnotator::cleanup_vcf() {
    if(vcf_header_in_)
        bcf_hdr_destroy(vcf_header_in_);
    if(vcf_fh_in_)
        bcf_close(vcf_fh_in_);
    if(vcf_header_out_)
        bcf_hdr_destroy(vcf_header_out_);
    if(vcf_fh_out_)
        bcf_close(vcf_fh_out_);
    if(vcf_record_)
        bcf_destroy(vcf_record_);
}

//Given a transcript ID and variant position,
//check if the variant is in a splice relevant region
//relevance depends on the user params
//intronic_min_distance_ and exonic_min_distance_
//The zero-based arithmetic is always fun.
//The variant object is one-based.
//GTF is one based
void VariantsAnnotator::get_variant_overlaps_spliceregion(const vector<BED>& exons,
                                                      AnnotatedVariant& variant) {
    variant.score = "-1";
    variant.annotation = "non_splice_region";
    string transcript_strand = exons[0].strand;
    //check if variant inside transcript coords for positive strand
    if(transcript_strand == "+" &&
       exons[0].start - intronic_min_distance_ > variant.start &&
       exons[exons.size() - 1].end + intronic_min_distance_ < variant.start) {
        return;
    }
    //check if variant inside transcript coords for negative strand
    if(transcript_strand == "-" &&
       exons[exons.size() - 1].start - intronic_min_distance_ > variant.start &&
       exons[0].end + intronic_min_distance_ < variant.start) {
        return;
    }
    for(std::size_t i = 0; i < exons.size(); i++) {
        //the rest of the exons are outside the junction - ps
        if(transcript_strand == "+" && exons[i].start - intronic_min_distance_ > variant.end) {
            return;
        }
        //the rest of the exons are outside the junction - ns
        if(transcript_strand == "-" && exons[i].end + intronic_min_distance_ < variant.end) {
            return;
        }
        //exonic near start
        if(variant.end >= exons[i].start &&
                variant.end <= exons[i].start + exonic_min_distance_) {
            variant.score =  common::num_to_str(variant.end - exons[i].start);
            variant.annotation = "splicing_exonic";
            return;
        }
        //intronic near start
        if(variant.end < exons[i].start &&
                variant.end >= exons[i].start - intronic_min_distance_) {
            variant.score = common::num_to_str(exons[i].start - variant.end);
            variant.annotation = "splicing_intronic";
            return;
        }
        //exonic near end
        if(variant.end <= exons[i].end &&
                variant.end >= exons[i].end - exonic_min_distance_) {
            variant.score = common::num_to_str(exons[i].end - variant.end);
            variant.annotation = "splicing_exonic";
            return;
        }
        //intronic near end
        if(variant.end > exons[i].end &&
                variant.end <= exons[i].end + intronic_min_distance_) {
            variant.score = common::num_to_str(variant.end - exons[i].end);
            variant.annotation = "splicing_intronic";
            return;
        }
    }
    return;
}

//Annotate one line of a VCF
//The line to be annotated is in vcf_record_
AnnotatedVariant VariantsAnnotator::annotate_record_with_transcripts(bool write_output) {
    string overlapping_genes = "NA",
           overlapping_transcripts = "NA",
           overlapping_distances = "NA",
           annotations = "NA";
    set<string> unique_genes;
    string chr = std::string(bcf_hdr_id2name(vcf_header_in_, vcf_record_->rid));
    AnnotatedVariant variant(chr, vcf_record_->pos, (vcf_record_->pos) + 1);
    //While calculating BINs, incorporate intronic_distance since transcripts
    //which lie within that distance will be relevant.
    BIN start_bin = ((vcf_record_->pos - intronic_min_distance_) >> _binFirstShift);
    BIN end_bin = ((vcf_record_->pos + intronic_min_distance_ ) >> _binFirstShift);
    //Iterate over all BINs this variant could fall under
    for (BINLEVEL i = 0; i < _binLevels; ++i) {
        BIN offset = _binOffsetsExtended[i];
        for (BIN b = (start_bin + offset); b <= (end_bin + offset); ++b) {
            vector<string> transcripts = gtf_.transcripts_from_bin(chr.c_str(), b);
            for(std::size_t i = 0; i < transcripts.size(); i++) {
                const vector<BED> & exons =
                    gtf_.get_exons_from_transcript(transcripts[i]);
                if(!exons.size()) {
                    throw runtime_error("Unexpected error. No exons for transcript "
                            + transcripts[i]);
                }
                //Skip single exon genes
                if(skip_single_exon_genes_ && exons.size() == 1) {
                    continue;
                }
                //Use a AnnotatedVariant object to hold the result
                get_variant_overlaps_spliceregion(exons, variant);
                if(variant.annotation != "non_splice_region") {
                    string gene_id = gtf_.get_gene_from_transcript(transcripts[i]);
                    //Use sign to encode intronic/exonic
                    string annotation = variant.annotation;
                    string dist_str = variant.score;
                    //Add gene only once for multiple transcripts of the same gene.
                    if(overlapping_transcripts != "NA") {
                        //Check if this gene is new
                        if(unique_genes.find(gene_id) == unique_genes.end()) {
                            overlapping_genes += "," + gene_id;
                            unique_genes.insert(gene_id);
                        }
                        overlapping_distances += "," + dist_str;
                        overlapping_transcripts += "," + transcripts[i];
                        annotations += "," + annotation;
                    } else {
                        overlapping_genes = gene_id;
                        overlapping_distances = dist_str;
                        overlapping_transcripts = transcripts[i];
                        unique_genes.insert(gene_id);
                        annotations = annotation;
                    }
                }
            }
        }
        start_bin >>= _binNextShift;
        end_bin >>= _binNextShift;
    }
    if(write_output) {
        if(bcf_update_info_string(vcf_header_out_, vcf_record_,
                    "genes", overlapping_genes.c_str()) < 0 ||
                bcf_update_info_string(vcf_header_out_, vcf_record_,
                    "transcripts", overlapping_transcripts.c_str()) < 0 ||
                bcf_update_info_string(vcf_header_out_, vcf_record_,
                    "distances", overlapping_distances.c_str()) < 0 ||
                bcf_update_info_string(vcf_header_out_, vcf_record_,
                    "annotations", annotations.c_str()) < 0) {
            throw runtime_error("Unable to update info string");
        }
        bcf_write(vcf_fh_out_, vcf_header_out_, vcf_record_);
    }
    return variant;
}

//Read in next record
bool VariantsAnnotator::read_next_record() {
    return (bcf_read(vcf_fh_in_, vcf_header_in_, vcf_record_) == 0);
}

//Heavylifting happens here.
void VariantsAnnotator::annotate_vcf() {
    load_gtf();
    open_vcf_in();
    open_vcf_out();
    while(read_next_record()) {
        annotate_record_with_transcripts();
    }
}
