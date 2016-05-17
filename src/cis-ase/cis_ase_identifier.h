/*  cis_ase_identifier.h -- 'cis-ase identify'

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

#ifndef CIS_ASE_IDENTIFIER_
#define CIS_ASE_IDENTIFIER_

#include <iostream>
#include <fstream>
#include <map>
#include "gtf_parser.h"
#include "htslib/sam.h"
#include "htslib/synced_bcf_reader.h"
#include "bam2bcf.h"
#include "bam_plcmd.h"
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "sample.h"
#include "variants_annotator.h"

using namespace std;

extern const double MIN_HET_PROB;
extern const double MIN_HOM_PROB;

extern "C" {
    void *bed_read(const char *fn);
    void bed_destroy(void *_h);
    int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len);
    void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
            int n, char *const*fn, int *n_plp,
            const bam_pileup1_t **plp, int ignore_rg);
    int bed_overlap(const void *_h, const char *chr, int beg, int end);
    int init_likelihoods(mplp_conf_t *conf, int n, char **fn, mplp_aux_t **data, bcf_callaux_t *bca, bcf_callret1_t *bcr, bcf_call_t *bc, mplp_pileup_t *gplp, htsFile *bcf_fp, bcf_hdr_t *bcf_hdr, bam_sample_t *sm, bam_hdr_t **h, mplp_ref_t *mp_ref);
    void set_data_iter(mplp_conf_t *conf, char** fn, mplp_aux_t **data, int *beg0, int *end0);
}

//Results of genotype call
struct genotype {
    //prob of het
    double p_het;
    //Read depth at the sites
    int n_reads;
    //Metadata describing the het type
    string het_type;
    genotype() {
        p_het = -1.0;
        n_reads = -1;
        het_type = "NA";
    }
    bool is_het(int min_depth) {
        if(n_reads == -1.0) {
            throw runtime_error("Uninitialized genotype!");
        }
        if(p_het >= MIN_HET_PROB && n_reads >= min_depth) {
            return true;
        }
        return false;
    }
    bool is_hom(int min_depth) {
        if(n_reads == -1.0) {
            throw runtime_error("Uninitialized genotype!");
        }
        if(1 - p_het >= MIN_HOM_PROB && n_reads >= min_depth) {
            return true;
        }
        return false;
    }
};

//results of mpileup for a variant
struct locus_info {
    //Flag for hom variants in RNA.
    bool is_hom_rna;
    //probability of hom genotype in RNA
    double p_hom_rna;
    //read-depth at the site - RNA
    int n_reads_rna;
    //Flag for het variants in DNA.
    bool is_het_dna;
    //probability of het variant in DNA.
    double p_het_dna;
    //read-depth at the site - DNA
    int n_reads_dna;
    locus_info() {
        p_het_dna = -1;
        p_hom_rna = -1;
        n_reads_dna = -1;
        n_reads_rna = -1;
        is_hom_rna = false;
        is_het_dna = false;
    }
};

//shared config variables for calling mpileup
//Doing this allows us to initialize these
//only once for a file and run on multiple
//positions. More efficient than the -r option
//which seems to iterate through the entire file
struct mpileup_conf_misc {
    bool result;
    char* file_names[1];
    int n_samples;
    mplp_aux_t **data;
    int i, tid, pos, *n_plp, beg0, end0, ref_len, max_depth;
    const bam_pileup1_t **plp;
    bam_mplp_t iter;
    bam_hdr_t *h; /* header of first file in input list */
    char *ref;
    bcf_callaux_t *bca;
    bcf_callret1_t *bcr;
    bcf_call_t bc;
    htsFile *bcf_fp;
    const char *mode;
    bcf_hdr_t *bcf_hdr;
    kstring_t buf;
    mplp_pileup_t gplp;
    bam_sample_t *sm;
    mplp_ref_t mp_ref;
    bcf1_t *bcf_rec;
    //set to true once initialized
    bool is_initialized;
    mpileup_conf_misc() {
        result = false;
        n_samples = 1;
        beg0 = 0;
        end0 = INT_MAX;
        bca = NULL;
        bcr = NULL;
        bcf_fp = NULL;
        bcf_hdr = NULL;
        file_names[0] = NULL;
        is_initialized = false;
        bcf_hdr = bcf_hdr_init("w");
        sm = bam_smpl_init();
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
        for (int i=0; i<sm->n; i++) {
            cerr << "\nAdding sample  " << sm->smpl[i];
            bcf_hdr_add_sample(bcf_hdr, sm->smpl[i]);
        }
        bcf_hdr_add_sample(bcf_hdr, NULL);
        //This needed to be extracted out from the vcf_hdr_write
        if ( bcf_hdr->dirty ) bcf_hdr_sync(bcf_hdr);
        mplp_ref_t mp_ref1 = MPLP_REF_INIT;
        mp_ref = mp_ref1;
        memset(&buf, 0, sizeof(kstring_t));
        plp = (const bam_pileup1_t**) calloc(n_samples, sizeof(bam_pileup1_t*));
        n_plp = (int *) calloc(n_samples, sizeof(int));
        data = (mplp_aux_t**) calloc(n_samples, sizeof(mplp_aux_t*));
        for (int i = 0; i < n_samples; ++i) {
            data[i] = (mplp_aux_t *) calloc(1, sizeof(mplp_aux_t));
        }
        bcr = (bcf_callret1_t *) calloc(n_samples, sizeof(bcf_callret1_t));
        bcf_rec = bcf_init1();
        memset(&gplp, 0, sizeof(mplp_pileup_t));
        // allocate data storage proportionate to number of samples being studied n_samples
        gplp.n = n_samples;
        gplp.n_plp = (int *) calloc(n_samples, sizeof(int));
        gplp.m_plp = (int *) calloc(n_samples, sizeof(int));
        gplp.plp = (bam_pileup1_t **) calloc(n_samples, sizeof(bam_pileup1_t*));
        memset(&bc, 0, sizeof(bcf_call_t));
        bc.PL = (int32_t *) malloc(15 * n_samples * sizeof(bc.PL));
        bc.DP4 = (int32_t *) malloc(n_samples * sizeof(int32_t) * 4);
        // all fmt_flag fields
        bc.fmt_arr = (uint8_t *) malloc(n_samples * sizeof(float));
        bc.ADR = (int32_t*) malloc((n_samples+1)*B2B_MAX_ALLELES*sizeof(int32_t));
        bc.ADF = (int32_t*) malloc((n_samples+1)*B2B_MAX_ALLELES*sizeof(int32_t));
        bc.bcf_hdr = bcf_hdr;
        bc.n = sm->n;
    }
    void init(string bam) {
        file_names[0] = strdup(bam.c_str());
        bam_smpl_add(sm, file_names[0], 0);
    }
    ~mpileup_conf_misc() {
        free(bc.tmp.s);
        bcf_destroy1(bcf_rec);
        bcf_hdr_destroy(bcf_hdr);
        free(bc.PL);
        free(bc.DP4);
        free(bc.ADR);
        free(bc.ADF);
        free(bc.fmt_arr);
        free(bcr);
        bam_smpl_destroy(sm); free(buf.s);
        for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
        free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);
        if (bcf_fp) {
            hts_close(bcf_fp);
            bcf_call_destroy(bca);
        }
        //These are allocated within the mpileup_with_likelihoods
        if(is_initialized) {
            bam_hdr_destroy(h);
            for (i = 0; i < n_samples; ++i) {
                sam_close(data[i]->fp);
            }
        }
        for (i = 0; i < n_samples; ++i) {
            free(data[i]);
        }
        free(data); free(plp); free(n_plp);
        free(mp_ref.ref[0]);
        free(mp_ref.ref[1]);
        if(file_names[0]) {
            free(file_names[0]);
        }
        is_initialized = false;
    }
};

//Set the configuration for mpileup
inline mplp_conf_t get_default_mpileup_conf(string ref, faidx_t *ref_fai) {
    mplp_conf_t mplp_conf;
    memset(&mplp_conf, 0, sizeof(mplp_conf_t));
    mplp_conf.min_baseQ = 13;
    mplp_conf.capQ_thres = 0;
    mplp_conf.max_depth = 250;
    mplp_conf.max_indel_depth = 250;
    mplp_conf.openQ = 40;
    mplp_conf.extQ = 20; mplp_conf.tandemQ = 100;
    mplp_conf.min_frac = 0.002; mplp_conf.min_support = 1;
    mplp_conf.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
    //uncompressed VCF
    mplp_conf.flag |= MPLP_BCF | MPLP_VCF | MPLP_NO_COMP;
    mplp_conf.fai = ref_fai;
    mplp_conf.fai_fname = (char*)ref.c_str();
    return mplp_conf;
}

//Free the mplp_conf_t object members
inline void free_mpileup_conf(mplp_conf_t conf) {
    if (conf.bed)
        bed_destroy(conf.bed);
    if (conf.reg)
        free(conf.reg);
}

//Workhorse for "cis-ase identify"
class CisAseIdentifier {
    private:
        //GTF parser object - holds the GTF in memory
        GtfParser gtf_parser_;
        //Minimum depth to consider somatic/ASE
        uint32_t min_depth_;
        //Window around somatic variants to look for transcripts
        uint32_t transcript_variant_window_;
        //VCF file with somatic variants
        string somatic_vcf_;
        //VCF file with polymorphisms
        string poly_vcf_;
        //Tumor RNAseq alignments
        string tumor_rna_;
        //Tumor DNA alignments
        string tumor_dna_;
        //Reference sequence FASTA
        string ref_;
        //GTF file with annotations
        string gtf_;
        //File to write output to
        string output_file_;
        //Which polymorphisms to look at
        string relevant_poly_annot_;
        //output stream to output annotated junctions file
        ofstream ofs_;
        //Somatic VCF file handle
        htsFile *somatic_vcf_fh_;
        //Somatic VCF Header
        bcf_hdr_t *somatic_vcf_header_;
        //Reference FASTA object
        faidx_t *ref_fai_;
        //mpileup conf for the somatic dna BAM
        mpileup_conf_misc somatic_dna_mmc_;
        //mpileup conf for the germline dna BAM
        mpileup_conf_misc germline_dna_mmc_;
        //mpileup conf for the germline rna BAM
        mpileup_conf_misc germline_rna_mmc_;
        //Configuration for somatic mpileup
        mplp_conf_t somatic_conf_;
        //Configuration for germline mpileup
        mplp_conf_t germline_conf_;
        //Somatic VCF record
        bcf1_t *somatic_vcf_record_;
        //Polymorphism VCF file handle
        htsFile *poly_vcf_fh_;
        //Polymorphism VCF Header
        bcf_hdr_t *poly_vcf_header_;
        //processed DNA snps - helps to pileup only once
        map<string, locus_info> dna_snps_;
        //processed RNA snps - helps to pileup only once
        map<string, locus_info> rna_snps_;
        //synced-reader for polymorphism vcf
        bcf_srs_t *poly_sr_;
        //Use binomial model for modeling ase?
        bool use_binomial_model_;
        //list of exonic variants indexed by "chr:BIN"
        map<string, vector<AnnotatedVariant> > bin_to_exonic_variants_;
    public:
        //Constructor
        CisAseIdentifier() : min_depth_(10),
                             transcript_variant_window_(1000),
                             somatic_vcf_("NA"),
                             tumor_rna_("NA"),
                             tumor_dna_("NA"), ref_("NA"), gtf_("NA"),
                             output_file_("NA"),
                             relevant_poly_annot_("exonic"),
                             somatic_vcf_fh_(NULL),
                             somatic_vcf_header_(NULL),
                             somatic_vcf_record_(NULL),
                             poly_vcf_fh_(NULL),
                             poly_vcf_header_(NULL),
                             use_binomial_model_(false) {
        }
        //Destructor
        ~CisAseIdentifier() {
            if(ofs_.is_open()) {
                ofs_.close();
            }
        }
        //Parse command line arguments
        void parse_options(int argc, char* argv[]);
        //Usage for this tool
        void usage(ostream &out);
        //The workhorse
        void run();
        //ASE identification starts here
        void identify_ase();
        //Open somatic VCF file
        void open_somatic_vcf();
        //init mpileup misc conf
        void mpileup_init1(string bam, mplp_conf_t *conf, mpileup_conf_misc& mmc1);
        //Run mpileup and get the genotype likelihoods
        bool mpileup_run(mplp_conf_t *conf, bool (CisAseIdentifier::*f)(bcf_hdr_t*, int, int, const bcf_call_t&, bcf1_t*), mpileup_conf_misc& mmc1);
        //Call genotypes using the binomial model - DNA
        genotype call_genotype_dna(const bcf_call_t& bc);
        //Call genotypes using the beta/binomial model - RNA
        genotype call_genotype_rna(const bcf_call_t& bc);
        //Get the SNPs within relevant window
        void process_snps_in_window(BED window);
        //Process homs in RNA(ASE)
        bool process_rna_hom(bcf_hdr_t* bcf_hdr, int tid, int pos, const bcf_call_t& bc, bcf1_t* bcf_rec);
        //Process somatic variants
        bool process_somatic_het(bcf_hdr_t* bcf_hdr, int tid, int pos, const bcf_call_t& bc, bcf1_t* bcf_rec);
        //Process hets in DNA
        bool process_germline_het(bcf_hdr_t* bcf_hdr, int tid, int pos, const bcf_call_t& bc, bcf1_t* bcf_rec);
        //Set the region as the region-string
        void set_mpileup_conf_region(mplp_conf_t & mplp_conf, string region);
        //Open the polymorphism VCF file
        void open_poly_vcf();
        //Free relevant pointers
        void cleanup();
        //Load pointer to reference
        void load_reference() {
            ref_fai_ = fai_load(ref_.c_str());
            if (ref_fai_ == NULL) throw runtime_error("Unable to open reference FASTA");
        }
        //Get the region of interest for a somatic variant
        //This depends on the transcripts in the region
        BED get_relevant_window(const char* chr, int pos);
        //Open BAM file-pointers
        void mpileup_init_all();
        //init mmcs
        void mmc_init_all();
        //create the map, where list of exonic variants are
        //indexed by "chr:bin"
        void annotate_exonic_polymorphisms();
};

#endif //CIS_ASE_IDENTIFIER_
