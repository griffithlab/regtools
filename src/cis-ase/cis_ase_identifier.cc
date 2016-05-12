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
#include "beta_model.h"
#include "binomial_model.h"
#include "common.h"
#include "cis_ase_identifier.h"
#include "gtf_utils.h"
#include "sample.h"
#include "samtools.h"

using namespace std;

//Minimum posterior probability to be considered het
const double MIN_HET_PROB = 0.5;
const double MIN_HOM_PROB = 0.5;

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
    out << "\n\t\t"   << "-B Use the binomial model for modeling ASE in RNAseq(uses the beta model by default)";
    out << "\n";
}

//Parse command line options
void CisAseIdentifier::parse_options(int argc, char* argv[]) {
    optind = 1; //Reset before parsing again.
    char c;
    while((c = getopt(argc, argv, "Bd:o:w:h")) != -1) {
        switch(c) {
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'B':
                use_binomial_model_ = true;
                break;
            case 'd':
                min_depth_ = atoi(optarg);
                break;
            case 'w':
                transcript_variant_window_ = atoi(optarg);
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
    cerr << "\nWindow around somatic-variants to look for transcripts: " <<
            transcript_variant_window_;
    if(use_binomial_model_) {
        cerr << "\nUsing the binomial model for modeling RNAseq ASE";
    }
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
    somatic_vcf_record_ = bcf_init();
}

//Set the region as the region-string
void CisAseIdentifier::set_mpileup_conf_region(mplp_conf_t &mplp_conf, string region) {
    mplp_conf.reg = strdup(region.c_str());
}

//open BAMfile pointers
void CisAseIdentifier::mpileup_init_all() {
    somatic_conf_ = get_default_mpileup_conf(ref_, ref_fai_);
    germline_conf_ = get_default_mpileup_conf(ref_, ref_fai_);
    mpileup_init1(tumor_dna_, &germline_conf_, germline_dna_mmc_);
    mpileup_init1(tumor_rna_, &germline_conf_, germline_rna_mmc_);
    mpileup_init1(tumor_dna_, &somatic_conf_, somatic_dna_mmc_);
}

//init mmcs
void CisAseIdentifier::mmc_init_all() {
    germline_dna_mmc_.init(tumor_dna_);
    somatic_dna_mmc_.init(tumor_dna_);
    germline_rna_mmc_.init(tumor_rna_);
}

//This function is ugly, pieces torn by hand from samtools
void CisAseIdentifier::mpileup_init1(string bam, mplp_conf_t *conf, mpileup_conf_misc& mmc1) {
    if ( conf->flag & MPLP_VCF )
        mmc1.mode = (conf->flag&MPLP_NO_COMP)? "wu" : "wz";   // uncompressed VCF or compressed VCF
    else
        mmc1.mode = (conf->flag&MPLP_NO_COMP)? "wub" : "wb";  // uncompressed BCF or compressed BCF
    mmc1.bcf_fp = bcf_open(conf->output_fname? conf->output_fname : "-", mmc1.mode);
    mmc1.bca = bcf_call_init(-1., conf->min_baseQ);
    mmc1.max_depth = conf->max_depth;
    init_likelihoods(conf, mmc1.n_samples, mmc1.file_names, mmc1.data, mmc1.bca, mmc1.bcr,
            &mmc1.bc, &mmc1.gplp, mmc1.bcf_fp, mmc1.bcf_hdr, mmc1.sm, &mmc1.h, &mmc1.mp_ref);
    mmc1.is_initialized = true; //all pointers initialized
}

//Run mpileup
bool CisAseIdentifier::mpileup_run(mplp_conf_t *conf,
        bool (CisAseIdentifier::*f)(bcf_hdr_t*, int, int, const bcf_call_t&, bcf1_t*),
                                    mpileup_conf_misc& mmc1) {
    bool result = false;
    //set the iterator to the region amongst other things
    set_data_iter(conf, mmc1.file_names, mmc1.data, &mmc1.beg0, &mmc1.end0);
    // begin pileup
    mmc1.iter = bam_mplp_init(mmc1.n_samples, mplp_func, (void**)mmc1.data);
    bam_mplp_set_maxcnt(mmc1.iter, mmc1.max_depth);
    bam_mplp_init_overlaps(mmc1.iter);
    while (bam_mplp_auto(mmc1.iter, &mmc1.tid, &mmc1.pos, mmc1.n_plp, mmc1.plp) > 0) {
        if (conf->reg && (mmc1.pos < mmc1.beg0 || mmc1.pos >= mmc1.end0)) continue; // out of the region requested
        if(!mmc1.h->target_name) printf("\nNot defined target\n");
        if (conf->bed && mmc1.tid >= 0 &&
            !bed_overlap(conf->bed, mmc1.h->target_name[mmc1.tid], mmc1.pos, mmc1.pos+1)) {
            continue;
        }
        if(conf->reg)
            cerr << "Region within run_mpileup " << conf->reg << endl;
        mplp_get_ref(mmc1.data[0], mmc1.tid, &mmc1.ref, &mmc1.ref_len);
        if (conf->flag & MPLP_BCF) {
            int total_depth, _ref0, ref16;
            for (int i = total_depth = 0; i < mmc1.n_samples; ++i) total_depth += mmc1.n_plp[i];
            group_smpl(&mmc1.gplp, mmc1.sm, &mmc1.buf, mmc1.n_samples, mmc1.file_names,
                       mmc1.n_plp, mmc1.plp, conf->flag & MPLP_IGNORE_RG);
            _ref0 = (mmc1.ref && mmc1.pos < mmc1.ref_len)? mmc1.ref[mmc1.pos] : 'N';
            ref16 = seq_nt16_table[_ref0];
            bcf_callaux_clean(mmc1.bca, &mmc1.bc);
            for (int i = 0; i < mmc1.gplp.n; ++i)
                bcf_call_glfgen(mmc1.gplp.n_plp[i], mmc1.gplp.plp[i], ref16, mmc1.bca, mmc1.bcr + i);
            mmc1.bc.tid = mmc1.tid; mmc1.bc.pos = mmc1.pos;
            bcf_call_combine(mmc1.gplp.n, mmc1.bcr, mmc1.bca, ref16, &mmc1.bc);
            bcf_clear1(mmc1.bcf_rec);
            bcf_call2bcf(&mmc1.bc, mmc1.bcf_rec, mmc1.bcr, conf->fmt_flag, 0, 0);
            result = (this->*f)(mmc1.bcf_hdr, mmc1.tid, mmc1.pos, mmc1.bc, mmc1.bcf_rec);
        }
    }
    //Destroy the iterator for this region
    bam_mplp_destroy(mmc1.iter);
    for (int i = 0; i < mmc1.n_samples; ++i) {
        if (mmc1.data[i]->iter) hts_itr_destroy(mmc1.data[i]->iter);
    }
    return result;
}

//Call genotypes using the posterior prob under the beta/binomial model for RNA
genotype CisAseIdentifier::call_genotype_rna(const bcf_call_t& bc) {
    genotype geno;
    //disregard sites with more than 5 alleles in the VCF &&
    //Check for minimum coverage
    geno.n_reads = bc.depth;
    if(bc.n_alleles <= 5 && bc.depth >= min_depth_) {
        if(use_binomial_model_) {
            calculate_binomial_phet(bc, geno);
        } else {
            BetaModel bm(bc);
            bm.calculate_beta_phet(geno);
        }
    }
    return geno;
}

//Call genotypes using the posterior prob under the binomial model for DNA hets
genotype CisAseIdentifier::call_genotype_dna(const bcf_call_t& bc) {
    genotype geno;
    //disregard sites with more than 5 alleles in the VCF &&
    //Check for minimum coverage
    geno.n_reads = bc.depth;
    if(bc.n_alleles <= 5 && bc.depth >= min_depth_) {
        calculate_binomial_phet(bc, geno);
    }
    return geno;
}

//Callback for germline het
bool CisAseIdentifier::process_germline_het(bcf_hdr_t* bcf_hdr, int tid,
                                            int pos, const bcf_call_t& bc, bcf1_t* bcf_rec) {
    string region = common::create_region_string(bcf_hdr_id2name(bcf_hdr, bcf_rec->rid), pos + 1, pos + 1);
    genotype geno = call_genotype_dna(bc);
    dna_snps_[region].p_het_dna = geno.p_het;
    dna_snps_[region].is_het_dna = false;
    if(geno.is_het(min_depth_)) {
        dna_snps_[region].is_het_dna = true;
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
    string region = common::create_region_string(bcf_hdr_id2name(bcf_hdr, bcf_rec->rid), pos + 1, pos + 1);
    genotype geno = call_genotype_rna(bc);
    rna_snps_[region].p_het_dna = geno.p_het;
    rna_snps_[region].is_het_dna = true;
    if(geno.is_hom(min_depth_)) {
        rna_snps_[region].is_het_dna = false;
        cerr << "RNA-hom" << endl;
    } else {
        cerr << "RNA variant is het" << endl;
    }
    cerr << "total, max " <<
        bcf_hdr_id2name(bcf_hdr, bcf_rec->rid) << " " <<
        pos + 1 << " " << geno.p_het << " " << geno.het_type << " " <<
        bcf_rec->d.als[0] << endl;
    return geno.is_hom(min_depth_);
}

//Get the window pertinent to this variant
//Get the transcripts within a certain distance from the somatic variant
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
                if(is_variant_within_transcript_window(exons, pos, transcript_strand,
                                            transcript_variant_window_)) {
                    int last_exon = exons.size() - 1;
                    if(exons[0].start < min_start) {
                        min_start = exons[0].start;
                    }
                    if(exons[last_exon].start < min_start) {
                        min_start = exons[last_exon].start;
                    }
                    if(exons[last_exon].end > max_end) {
                        max_end = exons[last_exon].end;
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
    genotype geno = call_genotype_dna(bc);
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
        //Check if SNP analyzed in RNA before
        if(rna_snps_.count(snp_region)) {
            cerr << endl << "Variant in map - already analyzed";
            if(!rna_snps_[snp_region].is_het_dna) {
                cerr << "rna is hom, now running DNA snp-mpileup" << endl;
                if(dna_snps_.count(snp_region)) {
                    if(dna_snps_[snp_region].is_het_dna) {
                        cerr << "DNA is het. potential ASE " << snp_region << endl;
                    } else {
                        cerr << "DNA not het" << endl;
                    }
                }
            } else {
                cerr << "rna not hom" << endl;
            }
            break;
        }
        cerr << "running rna mpileup" << endl;
        set_mpileup_conf_region(germline_conf_, snp_region);
        //Check if hom in RNA
        if(mpileup_run(&germline_conf_,
                       &CisAseIdentifier::process_rna_hom,
                       germline_rna_mmc_)) {
            cerr << "rna is hom, now running DNA snp-mpileup" << endl;
            //Check if het in DNA
            if(mpileup_run(&germline_conf_,
                           &CisAseIdentifier::process_germline_het,
                           germline_dna_mmc_)) {
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

//ASE identification starts here
void CisAseIdentifier::identify_ase() {
    while(bcf_read(somatic_vcf_fh_,
                   somatic_vcf_header_, somatic_vcf_record_) == 0) {
        string somatic_region = common::create_region_string(bcf_hdr_id2name(somatic_vcf_header_, somatic_vcf_record_->rid),
                                                             somatic_vcf_record_->pos+1, somatic_vcf_record_->pos+1);
        cerr << endl << "somatic region is " << somatic_region << endl;
        set_mpileup_conf_region(somatic_conf_, somatic_region);
        mpileup_run(&somatic_conf_,
           &CisAseIdentifier::process_somatic_het,
           somatic_dna_mmc_);//The workhorse
        free_mpileup_conf(somatic_conf_);
    }
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

//The workflow starts here
void CisAseIdentifier::run() {
    mmc_init_all(); //load all the mmcs
    load_reference(); //load reference genome
    gtf_parser_.load(); //load gene annotations
    open_somatic_vcf();
    open_poly_vcf();
    mpileup_init_all();
    //Start running the pileups and looking at GTs
    identify_ase();
    //Cleanup file handles
    cleanup();
}
