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
#include <cstring>
#include <htslib/sam.h>
#include "bam2bcf.h"
#include "bam_plcmd.h"
#include "cis_ase_identifier.h"
#include "hts.h"
#include "sample.h"
#include "samtools.h"

using namespace std;

//Usage for this tool
void CisAseIdentifier::usage(ostream& out) {
    out << "\nUsage:\t\t"
        << "regtools cis-splice-effects identify [options] somatic_variants.vcf"
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
    if(argc - optind >= 5) {
        somatic_vcf_ = string(argv[optind++]);
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

//Set the configuration for mpileup
void CisAseIdentifier::set_mpileup_conf() {
    memset(&mplp_conf_, 0, sizeof(mplp_conf_t));
    mplp_conf_.min_baseQ = 13;
    mplp_conf_.capQ_thres = 0;
    mplp_conf_.max_depth = 250;
    mplp_conf_.max_indel_depth = 250;
    mplp_conf_.openQ = 40;
    mplp_conf_.extQ = 20; mplp_conf_.tandemQ = 100;
    mplp_conf_.min_frac = 0.002; mplp_conf_.min_support = 1;
    mplp_conf_.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
    //uncompressed VCF
    mplp_conf_.flag |= MPLP_BCF | MPLP_VCF | MPLP_NO_COMP;
    mplp_conf_.bed = bed_read(somatic_vcf_.c_str());
    if (!mplp_conf_.bed) {
        throw runtime_error("Could not read file \"%s\"" + somatic_vcf_);
    }
    mplp_conf_.fai = fai_load(ref_.c_str());
    if (mplp_conf_.fai == NULL) throw runtime_error("Unable to open reference FASTA");
    mplp_conf_.fai_fname = (char*)ref_.c_str();
}

//Print somatic VCF record
//This function is ugly, pieces torn by hand from samtools
void CisAseIdentifier::run_mpileup() {
    char* alignment1 = strdup(tumor_dna_.c_str());
    char* file_names[] = { alignment1 };
    //mpileup(&mplp_conf_, 1, file_names);
    // init pileup
    int n_samples = 1;
    mplp_conf_t *conf = &mplp_conf_;
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
    mpileup_with_likelihoods(&mplp_conf_, 1, file_names, data, bca, bcr, &bc, &gplp, bcf_fp, bcf_hdr, sm, &h, &mp_ref);

    iter = bam_mplp_init(n_samples, mplp_func, (void**)data);
    bcr = (bcf_callret1_t *) calloc(sm->n, sizeof(bcf_callret1_t));
    if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(iter);
    max_depth = conf->max_depth;
    max_indel_depth = conf->max_indel_depth * sm->n;
    bam_mplp_set_maxcnt(iter, max_depth);
    bcf1_t *bcf_rec = bcf_init1();
    int ret;

    // begin pileup
    while ( (ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
        if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
        if(!h->target_name) printf("\nNot defined target\n");
        if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, h->target_name[tid], pos, pos+1)) { continue; }
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

            printf(" n_alleles %d ", bc.n_alleles);
            for (i=0; i< bc.n_alleles * (bc.n_alleles + 1) / 2; i++)
                printf(" PL %d ", bc.PL[i]);
            printf("\n");
            bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
        }
    }

    // clean up
    if(alignment1)
        free(alignment1);
    free(bc.tmp.s);
    bcf_destroy1(bcf_rec);
    if (bcf_fp)
    {
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
}

//Free relevant pointers
void CisAseIdentifier::cleanup() {
    if(somatic_vcf_header_)
        bcf_hdr_destroy(somatic_vcf_header_);
    if(somatic_vcf_fh_)
        bcf_close(somatic_vcf_fh_);
    if(somatic_vcf_record_)
        bcf_destroy(somatic_vcf_record_);
    if(mplp_conf_.fai)
        fai_destroy(mplp_conf_.fai);
    if (mplp_conf_.bed)
        bed_destroy(mplp_conf_.bed);
}

//The workhorse
void CisAseIdentifier::run() {
    open_somatic_vcf();
    set_mpileup_conf();
    run_mpileup();
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
