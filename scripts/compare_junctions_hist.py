import csv
from doctest import master
from itertools import groupby
import pandas as pd
from dfply import *
import numpy as np
import glob
import os

tag = 'E'
os.chdir('/Users/kcotto/Desktop/CHOL/')
splicing_variants_inputfile = '/Users/kcotto/Desktop/CHOL/all_splicing_variants_E.bed'
samples_inputfile = '/Users/kcotto/Desktop/CHOL/dir_names.tsv'

# read in all splicing variants
all_splicing_variants = pd.read_csv(splicing_variants_inputfile, delimiter='\t', header=0)

# create key to match regtools variant_info column and key2 that is the same as key but with sample name added
def createkey(row):
    key = row[0] + ':' + str(row[1]) + '-' + str(row[2])
    return key
all_splicing_variants['key'] = all_splicing_variants.apply(lambda row: createkey(row), axis=1)

def createkey(row):
    key = row[0] + ':' + str(row[1]) + '-' + str(row[2]) + '_' + row[3]
    return key
all_splicing_variants['key2'] = all_splicing_variants.apply(lambda row: createkey(row), axis=1)

# read in the sample names
all_samples = []
with open(samples_inputfile, 'r') as samples:
    reader = csv.reader(samples, delimiter='\t')
    for line in reader:
        all_samples.append(line[0])

### read in all of the regtools cse output for this cohort ###
# create list to hold each sample's df
dfs = []

# read each sample's output file into a df and subset columns, split variants into multirows,
# and require that variant is in all_splicing_variants
for sample in all_samples:
    path = f'samples/{sample}/output/cse_identify_filtered_compare_{tag}.tsv'
    df = f'df_{sample}'
    print(f'Reading in {sample}')
    df = pd.read_csv(path, delimiter='\t', header=0)
    df['sample'] = sample
    df = df[['sample', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor', 'score', 'name', 'genes']]
    df = df.dropna(subset=['variant_info'])
    df = df.set_index(['sample', 'chrom', 'start', 'end', 'strand', 'anchor', 'score', 'name', 'genes']).apply(lambda x: x.str.split(',').explode()).reset_index()
    df = df.loc[df['variant_info'].isin(all_splicing_variants['key'])]
    dfs.append(df)

# concat all individual dfs into one df
print("Concatenating each sample's df together")
master_df = pd.concat(dfs, axis=0, ignore_index=True)
del dfs

# create various keys
def createkey(row):
    key = row[1] + '_' + str(row[2]) + '_' + str(row[3]) + '_' + row[5] + '_' + row[9]
    return key
master_df['info'] = master_df.apply(lambda row: createkey(row), axis=1)

def createkey(row):
    key = row[9] + '_' + row[0]
    return key
master_df['key'] = master_df.apply(lambda row: createkey(row), axis=1)

def createkey(row):
    key = row[1] + '_' + str(row[2]) + '_' + str(row[3])
    return key
master_df['junction'] = master_df.apply(lambda row: createkey(row), axis=1)

# subset data to work on samples with splicing variant of interest
samples_w_variant_df = master_df.loc[master_df['key'].isin(all_splicing_variants['key2'])]
# print(samples_w_variant_df.info(verbose=True))

# start performing the calculations for this subset of data
print('Calculating for samples with variants of interest')
mode = 'blah'
print(samples_w_variant_df.head(10))
if mode == 'group':
    samples_w_variant_df = (samples_w_variant_df >>
                        group_by(X.key) >>
                        summarize(score_tmp = X.score.sum()) >>
                        outer_join(samples_w_variant_df, by='key')
                        )
    samples_w_variant_df['norm_score'] = samples_w_variant_df['score']/samples_w_variant_df['score_tmp']
    samples_w_variant_df = (samples_w_variant_df >>
                        group_by('junction') >>
                        summarize(mean_norm_score_variant=X.norm_score.mean(), sd_norm_score_variant=X.norm_score.std(), total_score_variant=X.score.sum()) >>
                        outer_join(samples_w_variant_df, by='junction')
                        )
    tmp_df = samples_w_variant_df.groupby('junction')[['norm_score', 'score', 'variant_info', 'sample', 'name', 'info']].aggregate(lambda x: x.tolist()).reset_index()
    samples_w_variant_df = pd.merge(samples_w_variant_df, tmp_df, on='junction')
    samples_w_variant_df = samples_w_variant_df[['sample_y', 'variant_info_y', 'chrom', 'start', 'end', 'strand', 'anchor',
                                             'info_y', 'genes', 'name_y', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                             'norm_score_y', 'score_y', 'junction', 'total_score_variant']]
    samples_w_variant_df.columns = ['samples', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                             'info', 'genes', 'names', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                             'norm_scores_variant', 'scores', 'junction', 'total_score_variant']
    samples_w_variant_df = samples_w_variant_df.astype(str).drop_duplicates()
else:
    samples_w_variant_df = (samples_w_variant_df >>
                            group_by(X.key) >>
                            summarize(score_tmp = X.score.sum()) >>
                            outer_join(samples_w_variant_df, by='key')
                            )
    samples_w_variant_df['norm_score'] = samples_w_variant_df['score']/samples_w_variant_df['score_tmp']
    samples_w_variant_df = (samples_w_variant_df >>
                        group_by('info') >>
                        summarize(mean_norm_score_variant=X.norm_score.mean(), sd_norm_score_variant=X.norm_score.std(), total_score_variant=X.score.sum()) >>
                        outer_join(samples_w_variant_df, by='info')
                        )
    tmp_df = samples_w_variant_df.groupby('info')[['norm_score', 'score', 'sd_norm_score_variant', 'mean_norm_score_variant', 'sample', 'name']].aggregate(lambda x: x.tolist()).reset_index()
    samples_w_variant_df = pd.merge(samples_w_variant_df, tmp_df, on='info')
    samples_w_variant_df = samples_w_variant_df[['sample_y', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                             'info', 'genes', 'name_y', 'mean_norm_score_variant_y', 'sd_norm_score_variant_y',
                                             'norm_score_y', 'score_y', 'junction', 'total_score_variant']]
    samples_w_variant_df.columns = ['samples', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                             'info', 'genes', 'names', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                             'norm_scores_variant', 'scores', 'junction', 'total_score_variant']

# work on samples that don't have the variant of interest

samples_wout_variant_df = master_df[-master_df['key'].isin(all_splicing_variants['key2'])]
samples_wout_variant_df = (samples_wout_variant_df >>
                        group_by(X.key) >>
                        summarize(score_tmp = X.score.sum()) >>
                        outer_join(samples_wout_variant_df, by='key')
                        )
samples_wout_variant_df['norm_score'] = samples_wout_variant_df['score']/samples_wout_variant_df['score_tmp']

mode = 'strict' #others include 'exclude' and 'group'

# if mode == 'strict':
    


