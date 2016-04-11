/*  cis_ase_identifier.cc -- 'cis-ase identify' methods

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

//Deal with the limits in htslib
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#include <stdexcept>
#include <sstream>
#include <cmath>
#include <cstring>
#include <htslib/sam.h>
#include <htslib/synced_bcf_reader.h>
#include "bam2bcf.h"
#include "bam_plcmd.h"
#include "common.h"
#include "cis_ase_identifier.h"
#include "sample.h"
#include "samtools.h"

using namespace std;

//Minimum posterior probability to be considered het
const double MIN_HET_PROB = 0.5;
const double MIN_HOM_PROB = 0.5;

//RR, RA1, A1A1, RA2, A1A2, A2A2, RA3, A1A3, A2A3, A3A3, RA4 ..
bool gt_het[15] = {0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0};

//Usage for this tool
void CisAseIdentifier::usage(ostream& out) {
    out << "\nUsage:\t\t"
        << "regtools cis-ase identify [options] somatic_variants.vcf polymorphism.vcf"
        << " tumor_dna_alignments.bam tumor_rna_alignments.bam ref.fa annotations.gtf";
    out << "\nOptions:";
    out << "\t"   << "-o STR Output file containing the aberrant splice junctions with annotations. [STDOUT]";
    out << "\n\t\t"   << "-d INT Minimum total read-depth for a somatic/ASE variant. [10]";
    out << "\n\t\t"   << "-w INT Window around a somatic variant to look for transcripts. ASE variants "
        << "will be in these transcripts[1000]";
    out << "\n";
}

//Parse command line options
void CisAseIdentifier::parse_options(int argc, char* argv[]) {
    optind = 1; //Reset before parsing again.
    char c;
    while((c = getopt(argc, argv, "d:o:w:h")) != -1) {
        switch(c) {
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'd':
                min_depth_ = atoi(optarg);
                break;
            case 'w':
                transcript_window_ = atoi(optarg);
                break;
            case 'h':
                usage(std::cerr);
                throw common::cmdline_help_exception("");
            default:
                usage(std::cerr);
                throw runtime_error("\nError parsing inputs!(1)");
        }
    }
    if(argc - optind >= 6) {
        somatic_vcf_ = string(argv[optind++]);
        poly_vcf_ = string(argv[optind++]);
        tumor_dna_ = string(argv[optind++]);
        tumor_rna_ = string(argv[optind++]);
        ref_ = string(argv[optind++]);
        gtf_ = string(argv[optind++]);
        gtf_parser_.set_gtffile(gtf_);
    }
    if(optind < argc ||
       somatic_vcf_ == "NA" ||
       tumor_dna_ == "NA" ||
       tumor_rna_ == "NA" ||
       ref_ == "NA" ||
       gtf_ == "NA"){
        usage(std::cerr);
        throw runtime_error("\nError parsing inputs!(2)\n");
    }
    cerr << "\nSomatic variants: " << somatic_vcf_;
    cerr << "\nPolymorphisms: " << poly_vcf_;
    cerr << "\nTumor DNA: " << tumor_dna_;
    cerr << "\nTumor RNA: " << tumor_rna_;
    cerr << "\nReference fasta file: " << ref_;
    cerr << "\nAnnotation file: " << gtf_;
    cerr << "\nMinimum read-depth for variants: " << min_depth_;
    cerr << endl;
}

//Open somatic VCF file
void CisAseIdentifier::open_somatic_vcf() {
    somatic_vcf_fh_ = bcf_open(somatic_vcf_.c_str(), "r");
    if(somatic_vcf_fh_ == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    somatic_vcf_header_ = bcf_hdr_read(somatic_vcf_fh_);
    if(somatic_vcf_header_ == NULL) {
        throw std::runtime_error("Unable to read header.");
    }
}

//Read in next record
bool CisAseIdentifier::read_somatic_record() {
    return (bcf_read(somatic_vcf_fh_, somatic_vcf_header_, somatic_vcf_record_) == 0);
}

//Set the region as the region-string
void CisAseIdentifier::set_mpileup_conf_region(mplp_conf_t &mplp_conf, string region) {
    mplp_conf.reg = strdup(region.c_str());
}

//Set the region as the somatic-vcf
void CisAseIdentifier::set_mpileup_conf_somatic_vcf(mplp_conf_t &mplp_conf) {
    mplp_conf.bed = bed_read(somatic_vcf_.c_str());
    if (!mplp_conf.bed) {
        throw runtime_error("Could not read file \"%s\"" + somatic_vcf_);
    }
}

//Init mpileup
bool CisAseIdentifier::mpileup_run(mplp_conf_t *conf,
        bool (CisAseIdentifier::*f)(bcf_hdr_t*, int, int, const bcf_call_t&, bcf1_t*),
                                    regtools_mpileup_conf& rmc1) {
    bool result = false;
    //set the iterator to the region amongst other things
    set_data_iter(conf, rmc1.file_names, rmc1.data, &rmc1.beg0, &rmc1.end0);
    // begin pileup
    rmc1.iter = bam_mplp_init(rmc1.n_samples, mplp_func, (void**)rmc1.data);
    bam_mplp_set_maxcnt(rmc1.iter, rmc1.max_depth);
    bam_mplp_init_overlaps(rmc1.iter);
    while (bam_mplp_auto(rmc1.iter, &rmc1.tid, &rmc1.pos, rmc1.n_plp, rmc1.plp) > 0) {
        if (conf->reg && (rmc1.pos < rmc1.beg0 || rmc1.pos >= rmc1.end0)) continue; // out of the region requested
        if(!rmc1.h->target_name) printf("\nNot defined target\n");
        if (conf->bed && rmc1.tid >= 0 &&
            !bed_overlap(conf->bed, rmc1.h->target_name[rmc1.tid], rmc1.pos, rmc1.pos+1)) {
            continue;
        }
        if(conf->reg)
            cerr << "Region within run_mpileup " << conf->reg << endl;
        mplp_get_ref(rmc1.data[0], rmc1.tid, &rmc1.ref, &rmc1.ref_len);
        if (conf->flag & MPLP_BCF) {
            int total_depth, _ref0, ref16;
            for (int i = total_depth = 0; i < rmc1.n_samples; ++i) total_depth += rmc1.n_plp[i];
            group_smpl(&rmc1.gplp, rmc1.sm, &rmc1.buf, rmc1.n_samples, rmc1.file_names,
                       rmc1.n_plp, rmc1.plp, conf->flag & MPLP_IGNORE_RG);
            _ref0 = (rmc1.ref && rmc1.pos < rmc1.ref_len)? rmc1.ref[rmc1.pos] : 'N';
            ref16 = seq_nt16_table[_ref0];
            bcf_callaux_clean(rmc1.bca, &rmc1.bc);
            for (int i = 0; i < rmc1.gplp.n; ++i)
                bcf_call_glfgen(rmc1.gplp.n_plp[i], rmc1.gplp.plp[i], ref16, rmc1.bca, rmc1.bcr + i);
            rmc1.bc.tid = rmc1.tid; rmc1.bc.pos = rmc1.pos;
            bcf_call_combine(rmc1.gplp.n, rmc1.bcr, rmc1.bca, ref16, &rmc1.bc);
            bcf_clear1(rmc1.bcf_rec);
            bcf_call2bcf(&rmc1.bc, rmc1.bcf_rec, rmc1.bcr, conf->fmt_flag, 0, 0);
            result = (this->*f)(rmc1.bcf_hdr, rmc1.tid, rmc1.pos, rmc1.bc, rmc1.bcf_rec);
        }
    }
    //Destroy the iterator for this region
    bam_mplp_destroy(rmc1.iter);
    for (int i = 0; i < rmc1.n_samples; ++i) {
        if (rmc1.data[i]->iter) hts_itr_destroy(rmc1.data[i]->iter);
    }
    return result;
}

//This function is ugly, pieces torn by hand from samtools
void CisAseIdentifier::mpileup_init(string bam, mplp_conf_t *conf, regtools_mpileup_conf& rmc1) {
    if ( conf->flag & MPLP_VCF )
        rmc1.mode = (conf->flag&MPLP_NO_COMP)? "wu" : "wz";   // uncompressed VCF or compressed VCF
    else
        rmc1.mode = (conf->flag&MPLP_NO_COMP)? "wub" : "wb";  // uncompressed BCF or compressed BCF
    rmc1.bcf_fp = bcf_open(conf->output_fname? conf->output_fname : "-", rmc1.mode);
    rmc1.bca = bcf_call_init(-1., conf->min_baseQ);
    rmc1.max_depth = conf->max_depth;
    mpileup_with_likelihoods(conf, rmc1.n_samples, rmc1.file_names, rmc1.data, rmc1.bca, rmc1.bcr,
            &rmc1.bc, &rmc1.gplp, rmc1.bcf_fp, rmc1.bcf_hdr, rmc1.sm, &rmc1.h, &rmc1.mp_ref);
    rmc1.is_initialized = true; //all pointers initialized
}

//Call genotypes using the posterior prob
genotype CisAseIdentifier::call_geno(const bcf_call_t& bc) {
    genotype geno;
    double sum_lik = 0, max_het_lik = 0;
    int n_gt = bc.n_alleles * (bc.n_alleles + 1) / 2;
    //disregard sites with more than 5 alleles in the VCF &&
    //Check for minimum coverage
    geno.n_reads = bc.depth;
    if(bc.n_alleles <= 5 && bc.depth >= min_depth_) {
        for (int i=0; i < n_gt; i++) {
            //convert back from phred
            double lik = pow(10.0, (-1.0 / 10.0 * bc.PL[i]));
            //printf(" PL %d lik %f", bc.PL[i], lik);
            sum_lik += lik;
            //True if GT is het
            if(gt_het[i]) {
                //perhaps switch from max to sum
                if (lik > max_het_lik) {
                    max_het_lik = lik;
                }
            }
        }
        geno.p_het = max_het_lik/sum_lik;
    }
    return geno;
}

//Callback for germline het
bool CisAseIdentifier::process_germline_het(bcf_hdr_t* bcf_hdr, int tid,
                                            int pos, const bcf_call_t& bc, bcf1_t* bcf_rec) {
    string region = common::create_region_string(bcf_hdr_id2name(bcf_hdr, bcf_rec->rid), pos + 1, pos + 1);
    germline_variants_[region].is_het_dna = false;
    genotype geno = call_geno(bc);
    germline_variants_[region].p_het_dna = geno.p_het;
    if(geno.is_het(min_depth_)) {
        germline_variants_[region].is_het_dna = true;
    } else {
        cerr << "Germline poly is hom" << endl;
    }
    cerr << "total, max " <<
        bcf_hdr_id2name(bcf_hdr, bcf_rec->rid) << " " <<
        pos + 1 << " " << geno.p_het << " " <<
        bcf_rec->d.als[0] << endl;
    return geno.is_het(min_depth_);
}

//Callback for hom in RNA(ASE)
bool CisAseIdentifier::process_rna_hom(bcf_hdr_t* bcf_hdr, int tid,
                                       int pos, const bcf_call_t& bc, bcf1_t* bcf_rec) {
    genotype geno = call_geno(bc);
    if(geno.is_hom(min_depth_)) {
        cerr << "RNA-hom" << endl;
    } else {
        cerr << "RNA variant is het" << endl;
    }
    cerr << "total, max " <<
        bcf_hdr_id2name(bcf_hdr, bcf_rec->rid) << " " <<
        pos + 1 << " " << geno.p_het << " " <<
        bcf_rec->d.als[0] << endl;
    return geno.is_hom(min_depth_);
}

//True if transcript within the variants transcript_window_
bool CisAseIdentifier::transcript_within_window(const vector<BED> &exons, uint32_t pos,
                                                string transcript_strand) {
    int n_exons = exons.size();
    if(transcript_strand == "+") {
        //variant inside transcript
        if(pos >= exons[0].start && pos <= exons[n_exons - 1].end) {
            return true;
        }
        //variant outside transcript
        if(exons[0].start - pos <= transcript_window_ &&
           exons[n_exons -1].start > pos) {
            return true;
        }
        //variant outside transcript
        if(pos - exons[n_exons - 1].end <= transcript_window_ &&
           exons[0].end < pos) {
            return true;
        }
    } else if(transcript_strand == "-") {
        //variant inside transcript
        if(pos >= exons[n_exons - 1].start && pos <= exons[0].end) {
            return true;
        }
        //variant outside transcript
        if(pos - exons[0].end <= transcript_window_ &&
           exons[n_exons -1].end < pos) {
            return true;
        }
        //variant outside transcript
        if(exons[n_exons - 1].start - pos <= transcript_window_ &&
           exons[0].start > pos) {
            return true;
        }
    } else {
        throw runtime_error("Unknown transcript strand.");
    }
    return false;
}

//Get the window pertinent to this variant
//Get the transcripts within a certain window
//Return the window that encompasses all these transcripts.
string CisAseIdentifier::get_relevant_window(const char* chr, int pos) {
    CHRPOS min_start = pos;
    CHRPOS max_end = pos;
    BIN start_bin = pos >> _binFirstShift;
    BIN end_bin = pos >> _binFirstShift;
    for (BINLEVEL i = 0; i < _binLevels; ++i) {
        BIN offset = _binOffsetsExtended[i];
        for (BIN b = (start_bin + offset); b <= (end_bin + offset); ++b) {
            vector<string> transcripts = gtf_parser_.transcripts_from_bin(chr, b);
            for(std::size_t i = 0; i < transcripts.size(); i++) {
                const vector<BED> & exons =
                    gtf_parser_.get_exons_from_transcript(transcripts[i]);
                //check if transcript within the window
                string transcript_strand = exons[0].strand;
                if(transcript_within_window(exons, pos, transcript_strand)) {
                    int n_exons = exons.size() - 1;
                    if(exons[0].start < min_start) {
                        min_start = exons[0].start;
                    }
                    if(exons[n_exons - 1].start < min_start) {
                        min_start = exons[n_exons - 1].start;
                    }
                    if(exons[exons.size() - 1].end > max_end) {
                        max_end = exons[exons.size() - 1].end;
                    }
                    if(exons[0].end > max_end) {
                        max_end = exons[0].end;
                    }
                }
            }
        }
        start_bin >>= _binNextShift;
        end_bin >>= _binNextShift;
    }
    return common::create_region_string(chr, min_start, max_end);
}

//Callback for somatic het
bool CisAseIdentifier::process_somatic_het(bcf_hdr_t* bcf_hdr, int tid,
                                           int pos, const bcf_call_t& bc, bcf1_t* bcf_rec) {
    genotype geno = call_geno(bc);
    if(geno.is_het(min_depth_)) {
        cerr << endl << "Somatic het. ";
        cerr << "total, max " <<
            bcf_hdr_id2name(bcf_hdr, bcf_rec->rid) << " " <<
            pos + 1 << " " << geno.p_het << " " <<
            bcf_rec->d.als[0];
        string window =
            get_relevant_window(bcf_hdr_id2name(bcf_hdr, bcf_rec->rid), pos);
        cerr << endl << "Window is ";
        cerr << window << endl;
        process_snps_in_window(window);
    } else {
        cerr << "Somatic variant is hom" << endl;
    }
    return geno.is_het(min_depth_);
}

//Open the polymorphism VCF file
void CisAseIdentifier::open_poly_vcf() {
    poly_vcf_fh_ = bcf_open(poly_vcf_.c_str(), "r");
    if(poly_vcf_fh_ == NULL) {
        throw std::runtime_error("Unable to open poly-vcf.");
    }
    poly_vcf_header_ = bcf_hdr_read(poly_vcf_fh_);
    if(poly_vcf_header_ == NULL) {
        throw std::runtime_error("Unable to read poly-vcf header.");
    }
}

//Get the information for SNPs within relevant window
void CisAseIdentifier::process_snps_in_window(string region) {
    std::cerr << "\ninside process_snps " << region << endl;
    if(!common::check_tabix_index(poly_vcf_)) {
        throw runtime_error("Tabix index does not exist for poly_vcf");
    }
    poly_sr_ = bcf_sr_init();
    bcf_sr_set_regions(poly_sr_, region.c_str(), 0);
    bcf_sr_add_reader(poly_sr_, poly_vcf_.c_str());
    while (bcf_sr_next_line(poly_sr_)) {
        bcf1_t *line = bcf_sr_get_line(poly_sr_, 0);
        string snp_region = common::create_region_string(bcf_hdr_id2name(poly_vcf_header_, line->rid),
                line->pos+1, line->pos+1);
        cerr << endl << "snp region is " << snp_region << endl;
        if(germline_variants_.count(snp_region)) {
            cerr << endl << "Variant exists in map ";
            cerr << germline_variants_[snp_region].is_het_dna;
            cerr << "\t";
            cerr << germline_variants_[snp_region].p_het_dna;
            cerr << endl;
            break;
        }
        cerr << "running rna mpileup" << endl;
        set_mpileup_conf_region(germline_conf_, snp_region);
        //Check if hom in RNA
        if(mpileup_run(&germline_conf_,
                    &CisAseIdentifier::process_rna_hom, germline_rna_rmc_)) {
            cerr << "rna is hom, now running RNA snp-mpileup" << endl;
            //Check if het in DNA
            if(mpileup_run(&germline_conf_,
                        &CisAseIdentifier::process_germline_het,
                        germline_dna_rmc_)) {
                cerr << "DNA is het. potential ASE " << snp_region << endl;
            } else {
                cerr << "DNA not het" << endl;
            }
        } else {
            cerr << "rna not hom" << endl;
        }
        free_mpileup_conf(germline_conf_);
        //bcf_destroy(line);
    }
    bcf_sr_destroy(poly_sr_);
}

//Free relevant pointers
void CisAseIdentifier::cleanup() {
    if(somatic_vcf_header_)
        bcf_hdr_destroy(somatic_vcf_header_);
    if(somatic_vcf_fh_)
        bcf_close(somatic_vcf_fh_);
    if(somatic_vcf_record_)
        bcf_destroy(somatic_vcf_record_);
    if(poly_vcf_header_)
        bcf_hdr_destroy(poly_vcf_header_);
    if(poly_vcf_fh_)
        bcf_close(poly_vcf_fh_);
    //Destroy pointer to reference
    if(ref_fai_)
        fai_destroy(ref_fai_);
}

//testing
void CisAseIdentifier::run2() {
    htsFile *test_bcf = NULL;
    bcf_hdr_t *test_header = NULL;
    bcf1_t *test_record = bcf_init();
    test_bcf = bcf_open(somatic_vcf_.c_str(), "r");
    if(test_bcf == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    test_header = bcf_hdr_read(test_bcf);
    if(test_header == NULL) {
        throw std::runtime_error("Unable to read header.");
    }
    while(bcf_read(test_bcf, test_header, test_record) == 0) {
        string somatic_region = common::create_region_string(bcf_hdr_id2name(test_header, test_record->rid),
                                                             test_record->pos+1, test_record->pos+1);
        cerr << endl << "somatic region is " << somatic_region << endl;
        set_mpileup_conf_region(somatic_conf_, somatic_region);
        mpileup_run(&somatic_conf_,
           &CisAseIdentifier::process_somatic_het,
           somatic_dna_rmc_);//The workhorse
        free_mpileup_conf(somatic_conf_);
    }
    bcf_hdr_destroy(test_header);
    bcf_destroy(test_record);
    bcf_close(test_bcf);
}

void CisAseIdentifier::run() {
    germline_dna_rmc_.init(tumor_dna_);
    somatic_dna_rmc_.init(tumor_dna_);
    germline_rna_rmc_.init(tumor_rna_);
    load_reference();
    //Load gene annotations
    gtf_parser_.load();
    somatic_conf_ = get_default_mpileup_conf();
    germline_conf_ = get_default_mpileup_conf();
    somatic_vcf_record_ = bcf_init();
    open_somatic_vcf();
    open_poly_vcf();
    //open file pointers
    mpileup_init(tumor_dna_, &germline_conf_, germline_dna_rmc_);
    mpileup_init(tumor_rna_, &germline_conf_, germline_rna_rmc_);
    mpileup_init(tumor_dna_, &somatic_conf_, somatic_dna_rmc_);
    run2();
    cleanup();
}
