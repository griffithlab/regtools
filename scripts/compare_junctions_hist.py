import csv
from itertools import groupby
import pandas as pd
from dfply import *
import numpy as np
from scipy import stats
import os
import argparse

# input_parser = argparse.ArgumentParser(
#     description="Run RegTools stats script",
# )
# input_parser.add_argument(
#     '-t',
#     '--tag',
#     help="Variant tag parameter used to run RegTools.",
# )
# input_parser.add_argument(
#     '-i',
#     '--variants_file',
#     help="File containing variants to be considered as splicing relevant."
# )
# input_parser.add_argument(
#     '-d',
#     '--dir_names',
#     help="File containing directory names corresponding to each sample that is to be processed."
# )
# input_parser.add_argument(
#     '-v',
#     '--variant-grouping',
#     help="",
#     choices=['strict', 'exclude', 'include']
# )

# args = input_parser.parse_args()

# tag = args.tag
# splicing_variants_inputfile = args.variants_file
# samples_inputfile = args.dir_names
# variant_grouping_mode = args.variant_grouping

tag = 'default'
splicing_variants_inputfile = '/Users/kcotto/Desktop/MET_samples/MET_splicing_variants.bed'
samples_inputfile = '/Users/kcotto/Desktop/MET_samples/samples.txt'
variant_grouping_mode = 'strict'
os.chdir('/Users/kcotto/Desktop/MET_samples/')

# read in all splicing variants
all_splicing_variants = pd.read_csv(
    splicing_variants_inputfile, delimiter='\t', header=None)

# create key to match regtools variant_info column and key2 that is the same as key but with sample name added


def createkey(row):
    key = row[0] + ':' + str(row[1]) + '-' + str(row[2])
    return key

all_splicing_variants['key'] = all_splicing_variants.apply(
    lambda row: createkey(row), axis=1)


def createkey(row):
    key = row[0] + ':' + str(row[1]) + '-' + str(row[2]) + '_' + row[3]
    return key


all_splicing_variants['key2'] = all_splicing_variants.apply(
    lambda row: createkey(row), axis=1)

# read in the sample names
all_samples = []
with open(samples_inputfile, 'r') as samples:
    reader = csv.reader(samples, delimiter='\t')
    for line in reader:
        all_samples.append(line[0])

num_of_samples = len(all_samples)

### read in all of the regtools cse output for this cohort ###
# create list to hold each sample's df
dfs = []

# read each sample's output file into a df and subset columns, split variants into multirows,
# and require that variant is in all_splicing_variants
for sample in all_samples:
    path = f'{sample}/output/cse_identify_filtered_compare_{tag}.tsv'
    df = f'df_{sample}'
    print(f'Reading in {sample}')
    df = pd.read_csv(path, delimiter='\t', header=0)
    df['sample'] = sample
    df = df[['sample', 'variant_info', 'chrom', 'start',
             'end', 'strand', 'anchor', 'score', 'name', 'genes']]
    df = df.dropna(subset=['variant_info'])
    df = df.set_index(['sample', 'chrom', 'start', 'end', 'strand', 'anchor', 'score',
                      'name', 'genes']).apply(lambda x: x.str.split(',').explode()).reset_index()
    df = df.loc[df['variant_info'].isin(all_splicing_variants['key'])]
    dfs.append(df)

# concat all individual dfs into one df
print("Concatenating each sample's df together")
master_df = pd.concat(dfs, axis=0, ignore_index=True)
del dfs

# create various keys


def createkey(row):
    key = row[1] + '_' + str(row[2]) + '_' + \
        str(row[3]) + '_' + row[5] + '_' + row[9]
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
samples_w_variant_df = master_df.loc[master_df['key'].isin(
    all_splicing_variants['key2'])]
# print(samples_w_variant_df.info(verbose=True))

def add_zeros_variant(row):
    norm_scores = row[1]
    if norm_scores == 0:
        norm_scores = [0]
    samples_wout_variant = row[2]
    samples_w_variant = row[3]    
    num_of_zeros_toadd = num_of_samples - samples_wout_variant - samples_w_variant
    zeros = np.repeat(0, num_of_zeros_toadd).tolist()
    norm_scores = norm_scores + zeros
    norm_scores.sort(reverse=True)
    new_norm_score_value = (',').join(map(str, norm_scores))
    return new_norm_score_value

# tmp_df['new_norm_scores'] = tmp_df.apply(lambda row: add_zeros_nonvariant(row), axis=1)

# start performing the calculations for this subset of data
print('Calculating normalized scores for samples with variants of interest')
if variant_grouping_mode == 'include':
    samples_w_variant_df = (samples_w_variant_df >>
                            group_by(X.key) >>
                            summarize(score_tmp=X.score.sum()) >>
                            outer_join(samples_w_variant_df, by='key')
                            )
    samples_w_variant_df['norm_score'] = samples_w_variant_df['score'] / \
        samples_w_variant_df['score_tmp']
    samples_w_variant_df = (samples_w_variant_df >>
                            group_by('junction') >>
                            summarize(mean_norm_score_variant=X.norm_score.mean(), sd_norm_score_variant=X.norm_score.std(), total_score_variant=X.score.sum()) >>
                            outer_join(samples_w_variant_df, by='junction')
                            )
    tmp_df = samples_w_variant_df.groupby('junction')[
        ['norm_score', 'score', 'variant_info', 'sample', 'name']].aggregate(lambda x: x.tolist()).reset_index()
    samples_w_variant_df = pd.merge(
        samples_w_variant_df, tmp_df, on='junction')
    samples_w_variant_df = samples_w_variant_df[['sample_y', 'variant_info_y', 'chrom', 'start', 'end', 'strand', 'anchor',
                                                 'info', 'genes', 'name_y', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                                 'norm_score_y', 'score_y', 'junction', 'total_score_variant']]
    samples_w_variant_df.columns = ['junction_samples', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                    'info', 'genes', 'names', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                    'norm_scores_variant', 'scores', 'junction', 'total_score_variant']
    samples_w_variant_df['variant_samples'] = samples_w_variant_df['junction_samples']
    samples_w_variant_df = samples_w_variant_df[~samples_w_variant_df.astype(
        str).duplicated()]
else:
    samples_w_variant_df = (samples_w_variant_df >>
                            group_by(X.key) >>
                            summarize(score_tmp=X.score.sum()) >>
                            outer_join(samples_w_variant_df, by='key')
                            )
    samples_w_variant_df['norm_score'] = samples_w_variant_df['score'] / \
        samples_w_variant_df['score_tmp']
    samples_w_variant_df = (samples_w_variant_df >>
                            group_by('info') >>
                            summarize(mean_norm_score_variant=X.norm_score.mean(), sd_norm_score_variant=X.norm_score.std(), total_score_variant=X.score.sum()) >>
                            outer_join(samples_w_variant_df, by='info')
                            )
    tmp_df = samples_w_variant_df.groupby('info')[['norm_score', 'score', 'sd_norm_score_variant',
                                                   'mean_norm_score_variant', 'sample', 'name']].aggregate(lambda x: x.tolist()).reset_index()
    samples_w_variant_df = pd.merge(samples_w_variant_df, tmp_df, on='info')
    samples_w_variant_df = samples_w_variant_df[['sample_x', 'sample_y', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                                 'info', 'genes', 'name_y', 'mean_norm_score_variant_y', 'sd_norm_score_variant_y',
                                                 'norm_score_y', 'score_y', 'junction', 'total_score_variant']]
    samples_w_variant_df.columns = ['sample_x', 'sample_y', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                    'info', 'genes', 'names', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                    'norm_scores_variant', 'scores', 'junction', 'total_score_variant']
    tmp_df = samples_w_variant_df.groupby('variant_info')[['sample_x']].aggregate(lambda x: list(set(x.tolist()))).reset_index()
    samples_w_variant_df = pd.merge(samples_w_variant_df, tmp_df, on='variant_info')
    samples_w_variant_df = samples_w_variant_df[['sample_y', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                                 'info', 'genes', 'names', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                                 'norm_scores_variant', 'scores', 'junction', 'total_score_variant', 'sample_x_y']]
    samples_w_variant_df.columns = ['junction_samples', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                    'info', 'genes', 'names', 'mean_norm_score_variant', 'sd_norm_score_variant',
                                    'norm_scores_variant', 'scores', 'junction', 'total_score_variant', 'variant_samples']
    samples_w_variant_df = samples_w_variant_df[~samples_w_variant_df.astype(
        str).duplicated()]

# work on samples that don't have the variant of interest
print('Calculating normalized scores for samples without variants of interest')
samples_wout_variant_df = master_df[~master_df['key'].isin(
    all_splicing_variants['key2'])]
del (master_df)

# mode = 'strict' #others include 'include' and 'exclude'
# if mode == 'strict':
samples_wout_variant_df = (samples_wout_variant_df >>
                           group_by(X.key) >>
                           summarize(score_tmp=X.score.sum()) >>
                           outer_join(samples_wout_variant_df, by='key')
                           )
samples_wout_variant_df['norm_score'] = samples_wout_variant_df['score'] / \
    samples_wout_variant_df['score_tmp']
samples_wout_variant_df = samples_wout_variant_df.loc[samples_wout_variant_df['variant_info'].isin(
    all_splicing_variants['key'])]
tmp_df = samples_wout_variant_df.groupby(
    'info')[['norm_score', 'sample']].aggregate(lambda x: x.tolist()).reset_index()
samples_wout_variant_df = pd.merge(samples_wout_variant_df, tmp_df, on='info')
samples_wout_variant_df['samples_wout_variant_count'] = samples_wout_variant_df['norm_score_y'].astype(
    str).str.count(',') + 1
if variant_grouping_mode == 'include' or variant_grouping_mode == 'exclude':
    samples_wout_variant_df = samples_wout_variant_df[~samples_wout_variant_df['junction'].isin(
        samples_w_variant_df['junction'])]
    tmp_df = samples_wout_variant_df.groupby(
        'info')[['norm_score_x']].aggregate(lambda x: x.tolist()).reset_index()
    samples_wout_variant_df = pd.merge(
        samples_wout_variant_df, tmp_df, on='info')
    samples_wout_variant_df = (samples_wout_variant_df >>
                           group_by('info') >>
                           summarize(total_score_non=X.score.sum()) >>
                           outer_join(samples_wout_variant_df, by='info')
                           )
    print(samples_wout_variant_df.info())
    samples_wout_variant_df = samples_wout_variant_df[['sample_y', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                                       'info', 'genes', 'norm_score_x_y', 'junction', 'total_score_non', 'samples_wout_variant_count']]
else:
    samples_wout_variant_df = (samples_wout_variant_df >>
                           group_by('info') >>
                           summarize(total_score_non=X.score.sum()) >>
                           outer_join(samples_wout_variant_df, by='info')
                           )
    samples_wout_variant_df = samples_wout_variant_df[['sample_y', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                                       'info', 'genes', 'norm_score_y', 'junction', 'total_score_non', 'samples_wout_variant_count']]
samples_wout_variant_df.columns = ['sample', 'variant_info', 'chrom', 'start', 'end', 'strand', 'anchor',
                                'info', 'genes', 'norm_scores_non', 'junction', 'total_score_non', 'samples_wout_variant_count']

print('Merging dataframes')
master_df = pd.merge(samples_w_variant_df, samples_wout_variant_df, how='left' ,on='info')
master_df = master_df[-master_df.astype(
        str).duplicated()]
del(samples_wout_variant_df)
del(samples_w_variant_df)

master_df['samples_w_variant_count'] = master_df['variant_samples'].astype(
    str).str.count(',') + 1

tmp_df = master_df[['info', 'norm_scores_non', 'samples_wout_variant_count', 'samples_w_variant_count']]
tmp_df = tmp_df.fillna(0)

def add_zeros_nonvariant(row):
    norm_scores = row[1]
    if norm_scores == 0:
        norm_scores = [0]
    samples_wout_variant = row[2]
    samples_w_variant = row[3]    
    num_of_zeros_toadd = num_of_samples - samples_wout_variant - samples_w_variant
    zeros = np.repeat(0, num_of_zeros_toadd).tolist()
    norm_scores = norm_scores + zeros
    norm_scores.sort(reverse=True)
    new_norm_score_value = (',').join(map(str, norm_scores))
    return new_norm_score_value

tmp_df['new_norm_scores'] = tmp_df.apply(lambda row: add_zeros_nonvariant(row), axis=1)
master_df = pd.merge(master_df, tmp_df, how='left' ,on='info')
del(tmp_df)

def get_mean(row):
    values = row[33].split(',')
    values = [float(i) for i in values]
    mean = np.mean(values)
    return mean

master_df['mean_norm_score_non'] = master_df.apply(lambda row: get_mean(row), axis=1)

def get_sd(row):
    values = row[33].split(',')
    values = [float(i) for i in values]
    std = np.std(values)
    return std

master_df['sd_norm_score_non'] = master_df.apply(lambda row: get_sd(row), axis=1)

print('getting p-values for associations')

def get_min(row):
    values = row[12]
    values = [float(i) for i in values]
    minimum = min(values)
    return(minimum)

master_df['min_norm_score_variant'] = master_df.apply(lambda row: get_min(row), axis=1)

def get_pvalue_mean(row):
    values = row[33].split(',')
    values = [float(i) for i in values]
    mean_value = row[10]
    pvalue = stats.percentileofscore(values, mean_value)
    pvalue = 1 - (pvalue/100.0)
    pvalue = re.sub('[\[\]]', '', np.array_str(pvalue))
    return pvalue

master_df['p_value_mean'] = master_df.apply(lambda row: get_pvalue_mean(row), axis=1)

def get_pvalue_min(row):
    values = row[33].split(',')
    values = [float(i) for i in values]
    min_value = row[36]
    pvalue = stats.percentileofscore(values, min_value)
    pvalue = 1 - (pvalue/100.0)
    pvalue = re.sub('[\[\]]', '', np.array_str(pvalue))
    return pvalue
    
master_df['p_value_min'] = master_df.apply(lambda row: get_pvalue_mean(row), axis=1)

master_df = master_df[['variant_samples', 'variant_info_x', 'genes_x', 'junction_samples', 
                       'chrom_x', 'start_x', 'end_x', 'strand_x', 'anchor_x', 'info', 
                       'names', 'mean_norm_score_variant', 'sd_norm_score_variant', 
                       'norm_scores_variant', 'total_score_variant', 'mean_norm_score_non',
                       'sd_norm_score_non', 'new_norm_scores', 'total_score_non', 'p_value_mean','p_value_min']]
master_df.columns = ['variant_samples', 'variant_info', 'genes', 'junction_samples', 
                       'chrom', 'start', 'end', 'strand', 'anchor', 'variant_junction_info', 
                       'names', 'mean_norm_score_variant', 'sd_norm_score_variant', 
                       'norm_scores_variant', 'total_score_variant', 'mean_norm_score_non',
                       'sd_norm_score_non', 'norm_scores_non', 'total_score_non', 'p_value_mean','p_value_min']

master_df = master_df.applymap(lambda x: x[0] if isinstance(x, list) else x)
master_df = master_df.fillna(0)

master_df.to_csv(f'junction_pvalues_{tag}_{variant_grouping_mode}.tsv', sep='\t', index=False)
print(master_df.info())