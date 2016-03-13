/*  bam_plcmd_regtools.c -- based on bam_plcmd.c, modified
 *  by Avinash Ramu.

    Copyright (C) 2008-2015 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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


#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include "bam_plcmd.h"
#include "sam_header.h"
#include "samtools.h"
#include "sam_opts.h"
#include "bam2bcf.h"
#include "sample.h"

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
                int n, char *const*fn, int *n_plp,
                const bam_pileup1_t **plp, int ignore_rg);
int mplp_func(void *data, bam1_t *b);
int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len);

/*
 * Performs pileup
 * @param conf configuration for this pileup
 * @param n number of files specified in fn
 * @param fn filenames
 */
int mpileup_with_likelihoods(mplp_conf_t *conf, int n, char **fn,
                             mplp_aux_t **data,
                             bcf_callaux_t *bca,
                             bcf_callret1_t *bcr,
                             bcf_call_t *bc,
                             mplp_pileup_t *gplp,
                             htsFile *bcf_fp,
                             bcf_hdr_t *bcf_hdr,
                             bam_sample_t *sm,
                             bam_hdr_t **h,
                             mplp_ref_t *mp_ref,
                             int *beg0,
                             int *end0
                             )
{
    extern void *bcf_call_add_rg(void *rghash, const char *hdtext, const char *list);
    extern void bcf_call_del_rghash(void *rghash);
    int i, max_depth, max_indel_depth;
    void *rghash = NULL;

    if (n == 0) {
        fprintf(stderr,"[%s] no input file/data given\n", __func__);
        exit(EXIT_FAILURE);
    }

    // read the header of each file in the list and initialize data
    for (i = 0; i < n; ++i) {
        bam_hdr_t *h_tmp;
        data[i] = calloc(1, sizeof(mplp_aux_t));
        data[i]->fp = sam_open_format(fn[i], "rb", &conf->ga.in);
        if ( !data[i]->fp )
        {
            fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, fn[i], strerror(errno));
            exit(EXIT_FAILURE);
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            exit(EXIT_FAILURE);
        }
        if (conf->fai_fname && hts_set_fai_filename(data[i]->fp, conf->fai_fname) != 0) {
            fprintf(stderr, "[%s] failed to process %s: %s\n",
                    __func__, conf->fai_fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        data[i]->conf = conf;
        data[i]->ref = mp_ref;
        h_tmp = sam_hdr_read(data[i]->fp);
        if ( !h_tmp ) {
            fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
            exit(EXIT_FAILURE);
        }
        bam_smpl_add(sm, fn[i], (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text);
        // Collect read group IDs with PL (platform) listed in pl_list (note: fragile, strstr search)
        rghash = bcf_call_add_rg(rghash, h_tmp->text, conf->pl_list);
        if (conf->reg) {
            hts_idx_t *idx = sam_index_load(data[i]->fp, fn[i]);
            if (idx == NULL) {
                fprintf(stderr, "[%s] fail to load index for %s\n", __func__, fn[i]);
                exit(EXIT_FAILURE);
            }
            if ( (data[i]->iter=sam_itr_querys(idx, h_tmp, conf->reg)) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s' with %s\n", __func__, conf->reg, fn[i]);
                exit(EXIT_FAILURE);
            }
            if (i == 0) *beg0 = data[i]->iter->beg, *end0 = data[i]->iter->end;
            hts_idx_destroy(idx);
        }
        else
            data[i]->iter = NULL;

        if (i == 0) *h = data[i]->h = h_tmp; // save the header of the first file
        else {
            // FIXME: check consistency between h and h_tmp
            bam_hdr_destroy(h_tmp);

            // we store only the first file's header; it's (alleged to be)
            // compatible with the i-th file's target_name lookup needs
            data[i]->h = *h;
        }
    }
    // allocate data storage proportionate to number of samples being studied sm->n
    gplp->n = sm->n;
    gplp->n_plp = calloc(sm->n, sizeof(int));
    gplp->m_plp = calloc(sm->n, sizeof(int));
    gplp->plp = calloc(sm->n, sizeof(bam_pileup1_t*));

    fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, sm->n, n);
    fprintf(stderr, "[%s]\n", conf->reg);
    // write the VCF header
    if (conf->flag & MPLP_BCF)
    {

        if (bcf_fp == NULL) {
            fprintf(stderr, "[%s] failed to write to %s: %s\n", __func__, conf->output_fname? conf->output_fname : "standard output", strerror(errno));
            exit(EXIT_FAILURE);
        }

        kstring_t str = {0,0,NULL};

        // Translate BAM @SQ tags to BCF ##contig tags
        // todo: use/write new BAM header manipulation routines, fill also UR, M5
        for (i=0; i<(*h)->n_targets; i++)
        {
            str.l = 0;
            ksprintf(&str, "##contig=<ID=%s,length=%d>", (*h)->target_name[i], (*h)->target_len[i]);
            bcf_hdr_append(bcf_hdr, str.s);
        }
        free(str.s);
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
        for (i=0; i<sm->n; i++)
            bcf_hdr_add_sample(bcf_hdr, sm->smpl[i]);
        bcf_hdr_add_sample(bcf_hdr, NULL);
        //This needed to be extracted out from the vcf_hdr_write
        if ( bcf_hdr->dirty ) bcf_hdr_sync(bcf_hdr);
        //bcf_hdr_write(bcf_fp, bcf_hdr);

        // Initialise the calling algorithm
        bca->rghash = rghash;
        bca->openQ = conf->openQ, bca->extQ = conf->extQ, bca->tandemQ = conf->tandemQ;
        bca->min_frac = conf->min_frac;
        bca->min_support = conf->min_support;
        bca->per_sample_flt = conf->flag & MPLP_PER_SAMPLE;

        bc->bcf_hdr = bcf_hdr;
        bc->n = sm->n;
        bc->PL = malloc(15 * sm->n * sizeof(*bc->PL));
        if (conf->fmt_flag)
        {
            assert( sizeof(float)==sizeof(int32_t) );
            bc->DP4 = malloc(sm->n * sizeof(int32_t) * 4);
            bc->fmt_arr = malloc(sm->n * sizeof(float)); // all fmt_flag fields
            if ( conf->fmt_flag&(B2B_INFO_DPR|B2B_FMT_DPR|B2B_INFO_AD|B2B_INFO_ADF|B2B_INFO_ADR|B2B_FMT_AD|B2B_FMT_ADF|B2B_FMT_ADR) )
            {
                // first B2B_MAX_ALLELES fields for total numbers, the rest per-sample
                bc->ADR = (int32_t*) malloc((sm->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
                bc->ADF = (int32_t*) malloc((sm->n+1)*B2B_MAX_ALLELES*sizeof(int32_t));
                for (i=0; i<sm->n; i++)
                {
                    bcr[i].ADR = bc->ADR + (i+1)*B2B_MAX_ALLELES;
                    bcr[i].ADF = bc->ADF + (i+1)*B2B_MAX_ALLELES;
                }
            }
        }
    }
    bcf_call_del_rghash(rghash);
    return 0;
}
