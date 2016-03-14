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
#include "hts.h"
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
        << "regtools cis-splice-effects identify [options] somatic_variants.vcf polymorphism.vcf"
        << " tumor_dna_alignments.bam tumor_rna_alignments.bam ref.fa annotations.gtf";
    out << "\nOptions:";
    out << "\t"   << "-o STR Output file containing the aberrant splice junctions with annotations. [STDOUT]";
    out << "\n\t\t" << "-v STR Output file containing variants annotated as splice relevant (VCF format).";
    out << "\n\t\t" << "-w INT\tWindow size in b.p to identify splicing events in. "
        << "\n\t\t\t" << "The tool identifies events in variant.start +/- w basepairs."
        << "\n\t\t\t" << "Default behaviour is to look at the window between previous and next exons.";
    out << "\n\t\t" << "-j STR Output file containing the aberrant junctions in BED12 format.";
    out << "\n";
}

//Parse command line options
void CisAseIdentifier::parse_options(int argc, char* argv[]) {
    optind = 1; //Reset before parsing again.
    char c;
    while((c = getopt(argc, argv, "o:w:v:j:h")) != -1) {
        switch(c) {
            case 'o':
                output_file_ = string(optarg);
                break;
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



//This function is ugly, pieces torn by hand from samtools
bool CisAseIdentifier::run_mpileup(string bam, mplp_conf_t *conf, bool (CisAseIdentifier::*f)(bcf_hdr_t*, int, int, const bcf_call_t&, bcf1_t*)) {
    bool result = false;
    char* alignment1 = strdup(bam.c_str());
    char* file_names[] = { alignment1 };
    // init pileup
    int n_samples = 1;
    mplp_aux_t **data;
    int i, tid, pos, *n_plp, beg0 = 0, end0 = INT_MAX, ref_len, max_depth, max_indel_depth;
    const bam_pileup1_t **plp;
    bam_mplp_t iter;
    bam_hdr_t *h = NULL; /* header of first file in input list */
    char *ref;

    bcf_callaux_t *bca = NULL;
    bcf_callret1_t *bcr = NULL;
    bcf_call_t bc;
    htsFile *bcf_fp = NULL;
    const char *mode;
    if ( conf->flag & MPLP_VCF )
        mode = (conf->flag&MPLP_NO_COMP)? "wu" : "wz";   // uncompressed VCF or compressed VCF
    else
        mode = (conf->flag&MPLP_NO_COMP)? "wub" : "wb";  // uncompressed BCF or compressed BCF
    bcf_fp = bcf_open(conf->output_fname? conf->output_fname : "-", mode);
    bcf_hdr_t *bcf_hdr = NULL;
    // BCF header creation
    bcf_hdr = bcf_hdr_init("w");

    bam_sample_t *sm = NULL;
    kstring_t buf;
    mplp_pileup_t gplp;

    memset(&buf, 0, sizeof(kstring_t));
    memset(&bc, 0, sizeof(bcf_call_t));

    data = (mplp_aux_t**) calloc(n_samples, sizeof(mplp_aux_t*));
    plp = (const bam_pileup1_t**) calloc(n_samples, sizeof(bam_pileup1_t*));
    n_plp = (int *) calloc(n_samples, sizeof(int));
    sm = bam_smpl_init();
    memset(&gplp, 0, sizeof(mplp_pileup_t));
    mplp_ref_t mp_ref = MPLP_REF_INIT;
    bca = bcf_call_init(-1., conf->min_baseQ);
    mpileup_with_likelihoods(conf, 1, file_names, data, bca, bcr, &bc, &gplp, bcf_fp, bcf_hdr, sm, &h, &mp_ref, &beg0, &end0);

    iter = bam_mplp_init(n_samples, mplp_func, (void**)data);
    bcr = (bcf_callret1_t *) calloc(sm->n, sizeof(bcf_callret1_t));
    if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(iter);
    max_depth = conf->max_depth;
    max_indel_depth = conf->max_indel_depth * sm->n;
    bam_mplp_set_maxcnt(iter, max_depth);
    bcf1_t *bcf_rec = bcf_init1();
    int ret;

    if(conf->reg)
        cerr << "\nRegion outside loop within run_mpileup " << conf->reg;
    // begin pileup
    while ((ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
        if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
        if(!h->target_name) printf("\nNot defined target\n");
        if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, h->target_name[tid], pos, pos+1)) { continue; }
        if(conf->reg)
            cerr << "\nRegion within run_mpileup " << conf->reg;
        mplp_get_ref(data[0], tid, &ref, &ref_len);
        //printf("tid=%d len=%d ref=%p/%s\n", tid, ref_len, ref, ref);
        if (conf->flag & MPLP_BCF) {
            int total_depth, _ref0, ref16;
            for (i = total_depth = 0; i < n_samples; ++i) total_depth += n_plp[i];
            group_smpl(&gplp, sm, &buf, n_samples, file_names, n_plp, plp, conf->flag & MPLP_IGNORE_RG);
            _ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
            ref16 = seq_nt16_table[_ref0];
            bcf_callaux_clean(bca, &bc);
            for (i = 0; i < gplp.n; ++i)
                bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], ref16, bca, bcr + i);
            bc.tid = tid; bc.pos = pos;
            bcf_call_combine(gplp.n, bcr, bca, ref16, &bc);
            bcf_clear1(bcf_rec);
            bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, 0, 0);
            result = (this->*f)(bcf_hdr, tid, pos, bc, bcf_rec);
            //bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
        }
    }

    // clean up
    if(alignment1)
        free(alignment1);
    free(bc.tmp.s);
    bcf_destroy1(bcf_rec);
    if (bcf_fp) {
        hts_close(bcf_fp);
        bcf_hdr_destroy(bcf_hdr);
        bcf_call_destroy(bca);
        free(bc.PL);
        free(bc.DP4);
        free(bc.ADR);
        free(bc.ADF);
        free(bc.fmt_arr);
        free(bcr);
    }
    bam_smpl_destroy(sm); free(buf.s);
    for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
    free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);
    bam_mplp_destroy(iter);
    bam_hdr_destroy(h);
    for (i = 0; i < n_samples; ++i) {
        sam_close(data[i]->fp);
        if (data[i]->iter) hts_itr_destroy(data[i]->iter);
        free(data[i]);
    }
    free(data); free(plp); free(n_plp);
    free(mp_ref.ref[0]);
    free(mp_ref.ref[1]);
    return result;
}

//Call genotypes using the posterior prob
genotype CisAseIdentifier::call_geno(const bcf_call_t& bc) {
    genotype geno;
    double sum_lik = 0, max_het_lik = 0;
    int n_gt = bc.n_alleles * (bc.n_alleles + 1) / 2;
    //disregard sites with more than 5 alleles in the VCF
    if(bc.n_alleles <= 5) {
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
        geno.determine_het();
    }
    return geno;
}

//Callback for germline het
bool CisAseIdentifier::process_germline_het(bcf_hdr_t* bcf_hdr, int tid, int pos, const bcf_call_t& bc, bcf1_t* bcf_rec) {
    string region = common::create_region_string(bcf_hdr_id2name(bcf_hdr, bcf_rec->rid), pos + 1, pos + 1);
    germline_variants_[region].is_het_dna = false;
    genotype geno = call_geno(bc);
    germline_variants_[region].p_het_dna = geno.p_het;
    if(geno.is_het) {
        fprintf(stderr, "\ngermline-het-dna chr, position, total, max, als1, als2 %s %d %c %c",
                bcf_hdr_id2name(bcf_hdr, bcf_rec->rid),
                pos + 1, bcf_rec->d.als[0], bcf_rec->d.als[1]);
        germline_variants_[region].is_het_dna = true;
    }
    return geno.is_het;
}

//Callback for somatic het
bool CisAseIdentifier::process_somatic_het(bcf_hdr_t* bcf_hdr, int tid, int pos, const bcf_call_t& bc, bcf1_t* bcf_rec) {
    genotype geno = call_geno(bc);
    if(geno.is_het) {
        fprintf(stderr, "\nsomatic-var chr, position, "
                "total, max %s %d %f %c %c",
                bcf_hdr_id2name(bcf_hdr, bcf_rec->rid),
                pos + 1, geno.p_het,
                bcf_rec->d.als[0], bcf_rec->d.als[1]);
        cerr << endl << "Somatic het. Window is ";
        string window =
            common::create_region_string(bcf_hdr_id2name(bcf_hdr,
                        bcf_rec->rid), pos - 1000, pos + 1000);
        cerr << window << endl;
        process_snps_in_window(window);
    }
    return geno.is_het;
}

//Get the information for SNPs within relevant window
void CisAseIdentifier::process_snps_in_window(string region) {
    htsFile *bcf_fp = NULL;
    bcf_hdr_t *bcf_hdr = NULL;
    bcf_fp = bcf_open(poly_vcf_.c_str(), "r");
    if(bcf_fp == NULL) {
        throw std::runtime_error("Unable to open file.");
    }
    bcf_hdr = bcf_hdr_read(bcf_fp);
    if(bcf_hdr == NULL) {
        throw std::runtime_error("Unable to read header.");
    }
    bcf_srs_t *sr = bcf_sr_init();
    bcf_sr_set_regions(sr, region.c_str(), 0);
    bcf_sr_add_reader(sr, poly_vcf_.c_str());
    std::cerr << "\nchromosome\tposition\tnum_alleles" << std::endl;
    while (bcf_sr_next_line(sr)) {
        bcf1_t *line = bcf_sr_get_line(sr, 0);
        string snp_region = common::create_region_string(bcf_hdr_id2name(bcf_hdr, line->rid), line->pos+1, line->pos+1);
        cerr << endl << "snp region is " << snp_region << endl;
        if(germline_variants_.count(snp_region)) {
            cerr << endl << "Variant in map ";
            cerr << germline_variants_[snp_region].is_het_dna;
            cerr << "\t";
            cerr << germline_variants_[snp_region].p_het_dna;
            cerr << endl;
            return;
        }
        set_mpileup_conf_region(germline_conf_, snp_region);
        run_mpileup(tumor_dna_, &germline_conf_,
                    &CisAseIdentifier::process_germline_het);
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
    free_mpileup_conf(somatic_conf_);
    free_mpileup_conf(germline_conf_);
}

//The workhorse
void CisAseIdentifier::run() {
    load_reference();
    somatic_conf_ = get_default_mpileup_conf();
    germline_conf_ = get_default_mpileup_conf();
    somatic_vcf_record_ = bcf_init();
    open_somatic_vcf();
    set_mpileup_conf_somatic_vcf(somatic_conf_);
    run_mpileup(tumor_dna_, &somatic_conf_,
                &CisAseIdentifier::process_somatic_het);
    cleanup();
}

/*    cerr << endl << "In Run";
    for each_dna_somatic_variant {
        if variant_is_het_in_tumor_dna {
            for each_nearby_germline_polymorphism {
                if poly_is_het_tumor_dna and poly_is_not_het_in_tumor_rna {
                    print dna_somatic_variant "\t" dna_poly "\t" rna_poly "\t" evidence_level"
                    add_poly_to_list_of_ase_polys();
                }
            }
        }
    }
*/
