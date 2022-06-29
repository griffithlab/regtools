# compare_junctions_hist.R
# Rscript --vanilla compare_junctions_hist.R <tag>


# load libraries
library(data.table)
library(plyr)
library(tidyverse)

# get options tag
# args = commandArgs(trailingOnly = TRUE)
# wd = args[1]
# samples = args[2]
# splice_variants_file = args[3]

setwd('wd')
# setwd(wd)

sample_names = 'all_samples.txt'
input_file = 'all_tumor_splice_variants.bed'
tumor_sample_barcodes = Sys.glob('starting_barcode_files/*tumor*.txt')
normal_sample_barcodes = Sys.glob('starting_barcode_files/*normal*.txt')

# All splicing relevant variants (union of rows from variants.bed files; add column with comma-separated list of sample names)
all_splicing_variants = unique(data.table::fread(input_file), sep = '\t', header = T, stringsAsFactors = FALSE)
colnames(all_splicing_variants) <- c("chrom", "start", "end", "samples")

# all_splicing_variants = as.data.table(aggregate(samples ~ chrom + start + end, paste, data=all_splicing_variants)) # I don't think this is doing anything?
all_splicing_variants$key <- paste0(all_splicing_variants$chrom, ":", all_splicing_variants$start, "-", all_splicing_variants$end) #this key is just a 1bp-long chrom:start-stop designed to match the regtools output variant_info column

## Get all of the samples
samples_w_data = strsplit(scan(sample_names, what="", sep="\n"), "[[:space:]]+")

## Get metrics
tumor_samples <- 
  tumor_sample_barcodes %>%
  map_df(~fread(.,header=F))
normal_samples <- 
  normal_sample_barcodes %>%
  map_df(~fread(.,header=F))
number_of_starting_tumor_cells = dim(tumor_samples)[1]
number_of_starting_normal_cells = dim(normal_samples)[1]
number_of_tumor_cells_w_data = sum(str_count(samples_w_data, pattern = 'tumor'))
number_of_normal_cells_w_data = sum(str_count(samples_w_data, pattern = 'normal'))

print(paste("Number of starting tumor cells:", number_of_starting_tumor_cells))
print(paste("Number of tumor cells with data:", number_of_tumor_cells_w_data))
print(paste("Number of starting normal cells:", number_of_starting_normal_cells))
print(paste("Number of normal cells with data:", number_of_normal_cells_w_data))


################################################################################
##### Helper functions #########################################################
################################################################################

# basically combines all the cse_identify_filtered_compare files for this cohort
get_sample_data <- function(sample, all_splicing_variants){
  file_name = paste(sample, ".tsv", sep = "")
  cse_identify_data = data.table::fread(file_name, sep = '\t', header = T, stringsAsFactors = FALSE, strip.white = TRUE)
  cse_identify_data$sample <- sample
  cse_identify_data <- cse_identify_data[,.(sample,variant_info,chrom,start,end,strand,anchor,score,name,gene_names)]
  cse_identify_data <- cse_identify_data[!is.na(variant_info)]
  cse_identify_data <- cse_identify_data[, variant_info := as.list(strsplit(variant_info, ",", fixed=TRUE))]
  cse_identify_data <- cse_identify_data[rep(cse_identify_data[,.I], lengths(variant_info))][, variant_info := unlist(cse_identify_data$variant_info)][]
  cse_identify_data <- cse_identify_data[variant_info %chin% all_splicing_variants$key]
  return(cse_identify_data)
}

################################################################################
##### Combine cse_identify_filtered_compare.tsv files for each sample ##########
################################################################################

# this is regtools compare output across samples (so it contains variant-junction lines even from samples without variant)
dt <- rbindlist(lapply(samples_w_data, get_sample_data, all_splicing_variants=all_splicing_variants))

dt[,info := paste(chrom, start, end, anchor, variant_info, sep="_")]

################ subset out entries where we can for now #######################

# make a sample/variant_info key
dt[,key := paste0(variant_info, "_", sample)]

print("zl")
cse_identify_v1 <- dt
rm(dt)
print("kc")
cse_identify_v2 <- cse_identify_v1
print("yy")



######### aggregate scores for variant/samples identified in cse_identify ######

# key files with variant and sample keys
all_splicing_variants$key2 <- paste0(all_splicing_variants$key, "_", all_splicing_variants$samples)

# work on variants which have a sample to go with them first
cse_identify_v1 <- cse_identify_v1[key %chin% all_splicing_variants$key2]

# perform the required calculations for this subset
cse_identify_v1[,score.tmp := sum(score), by=.(sample, variant_info)]
cse_identify_v1[,score_norm := score/score.tmp, by=.(sample, variant_info)]
cse_identify_v1[,mean_norm_score_variant := mean(score_norm), by=.(sample,variant_info,chrom,start,end,strand,anchor,info)]
cse_identify_v1[,sd_norm_score_variant := sd(score_norm), .(variant_info,chrom,start,end,strand,anchor,info)]
cse_identify_v1[,total_score_variant := sum(score), .(variant_info,chrom,start,end,strand,anchor,info)]
cse_identify_v1 <- cse_identify_v1
cse_identify_v1_tmp = ddply(cse_identify_v1, .(variant_info, chrom, start, end, strand, anchor, info), summarise, mean_norm_score_variant=mean(mean_norm_score_variant), name=paste(name,collapse = ','), sample=paste(sample,collapse = ','), norm_scores_variant=paste(score_norm, collapse = ','))
cse_identify_v1 = merge(cse_identify_v1, cse_identify_v1_tmp, by=c('variant_info', 'chrom', "start", "end", "strand", "anchor", 'info'))
rm(cse_identify_v1_tmp)

print("test4")

# subset and rename columns to match the original output
cse_identify_v1 <- cse_identify_v1[,c("sample.y", "variant_info", "chrom", "start", "end", "strand", "anchor",
                                      "variant_info", "info", "gene_names","name.y", "mean_norm_score_variant.y",
                                      "sd_norm_score_variant", "norm_scores_variant", "total_score_variant")]
colnames(cse_identify_v1) <- c("sample", "key", "chrom", "start", "end", "strand", "anchor", "variant_info",
                               "info", "gene_names", "names", "mean_norm_score_variant", "sd_norm_score_variant",
                               "norm_scores_variant", "total_score_variant")

################ agrregate samples with no variant of interest ############################

# second, we just want entries where the variant is not in the sample we care about
cse_identify_v2 <- cse_identify_v2[!key %chin% all_splicing_variants$key2]

print("test5")

cse_identify_v2[,score.tmp := sum(score), by=.(sample, variant_info)]
cse_identify_v2[,norm_score := score/score.tmp, by=.(sample, variant_info)]

# work on variants which have a sample to go with them first
cse_identify_v2 <- cse_identify_v2[variant_info %chin% all_splicing_variants$key]


# this will work for now but can be improved I think
a <- function(x){
  x <- x[,list(mean_norm_score_non=mean(norm_score),
               sd_norm_score_non=sd(norm_score),
               norm_scores_non=list(norm_score),
               total_score_non=sum(score)),
         by=list(chrom,start,end,strand,anchor,variant_info,info)]
  return(x)
}
cse_identify_v2 <- split(cse_identify_v2, cse_identify_v2$variant_info)
cse_identify_v2 <- lapply(cse_identify_v2, a)
cse_identify_v2 <- rbindlist(cse_identify_v2)

print("test6")

################ Merge the two DT's we've been working on together #############

regtools_data <- merge(cse_identify_v1, cse_identify_v2, by="info", all.x = TRUE)
rm(cse_identify_v1)
rm(cse_identify_v2)

tmp <- all_splicing_variants[order(key)]
tmp <- ddply(tmp, .(chrom, start, end, key), summarise, samples=paste(samples,collapse = ','))
regtools_data = merge(regtools_data, tmp, by.x=c('variant_info.x'), by.y=c('key'))
columns_to_keep = c('samples', 'variant_info.x', 'sample', "info", "chrom.x", "start.x", "end.x", 'strand.x', 'anchor.x', 'gene_names',
                    'names', 'mean_norm_score_variant', 'sd_norm_score_variant', 'norm_scores_variant',
                    'total_score_variant', 'mean_norm_score_non', 'sd_norm_score_non', 'norm_scores_non',
                    'total_score_non')
regtools_data = subset(regtools_data, select=columns_to_keep)

print('this works')

# zeroes need to be added in for some samples
a <- function(x, y){
  toAdd <- y - length(x)
  # browser()
  toAdd <- rep(0.0000000, toAdd)
  x <- c(x, toAdd)
  return(x)
}

#this needs to be thought about
x <- mapply(a, regtools_data$norm_scores_non, number_of_starting_normal_cells)

get_mean <- function(x){
  x <- mean(as.numeric(x))
  return(x)
}

get_sd <- function(x){
  x <- sd(as.numeric(x))
  return(x)
}

x <- mapply(get_mean, regtools_data$norm_scores_non) 
regtools_data$mean_norm_score_non <- x

x <- mapply(get_sd, regtools_data$norm_scores_non) 
regtools_data$sd_norm_score_non <- x

a <- function(x, y){
  # if(y == "TCGA-ZH-A8Y2-01A,TCGA-ZH-A8Y5-01A"){
  #    browser()
  # }
  x <- unlist(strsplit(x, ","))
  toAdd <- y - length(x) 
  # browser()
  if (toAdd > 0) {
    toAdd <- rep(0.0000000, toAdd)
    x <- c(x, toAdd)
  }
  x <- list(x)
  return(x)
}

x <- mapply(a, regtools_data$norm_scores_variant, number_of_starting_tumor_cells)
regtools_data$norm_scores_variant = x

x <- mapply(get_mean, regtools_data$norm_scores_variant) 
regtools_data$mean_norm_score_variant <- x

x <- mapply(get_sd, regtools_data$norm_scores_variant) 
regtools_data$sd_norm_score_variant <- x

get_min <- function(x){
  x <- min(as.numeric(x))
  return(x)
}

x <- mapply(get_min, regtools_data$norm_scores_variant) 
regtools_data$min_norm_score_variant <- x

print("test7")

################ calculate p-values ############################################

# calculate using mean
a <- function(x){
  variant_norm_score = mean(as.numeric(unlist(strsplit(x[['norm_scores_variant']], ',', fixed=TRUE))))
  if(length(x[['norm_scores_non']]) <= 1){
    return(0)
  }
  
  all_norm_scores = c(x$norm_scores_non, variant_norm_score)
  countable = rank(all_norm_scores)
  num_samples = str_count(x$norm_scores_variant, ',') + 1
  non_variant_norm_scores_ranked = head(countable, -1)
  variant_norm_score_ranked = tail(countable, 1)
  histinfo = hist(non_variant_norm_scores_ranked, 
                  breaks = seq(0.5, max(non_variant_norm_scores_ranked)+1.5, by=1), plot=F)
  mids = histinfo$mids
  cd = cumsum(histinfo$density)
  underestimate = max(which(mids <= variant_norm_score_ranked))
  pvalue = 1-cd[underestimate]
  return(pvalue)
}

# calculate using min
b <- function(x){
  variant_norm_score = min(as.numeric(unlist(strsplit(x[['norm_scores_variant']], ',', fixed=TRUE))))
  if(variant_norm_score == 0){
    return(1)
  }
  if(length(x[['norm_scores_non']]) <= 1){
    return(0)
  }
  
  all_norm_scores = c(x$norm_scores_non, variant_norm_score)
  countable = rank(all_norm_scores)
  num_samples = str_count(x$norm_scores_variant, ',') + 1
  non_variant_norm_scores_ranked = head(countable, -1)
  variant_norm_score_ranked = tail(countable, 1)
  histinfo = hist(non_variant_norm_scores_ranked, 
                  breaks = seq(0.5, max(non_variant_norm_scores_ranked)+1.5, by=1), plot=F)
  mids = histinfo$mids
  cd = cumsum(histinfo$density)
  underestimate = max(which(mids <= variant_norm_score_ranked))
  pvalue = 1-cd[underestimate]
  return(pvalue)
}

regtools_data$p_value_mean <- apply(regtools_data, 1, a)
regtools_data$p_value_min <- apply(regtools_data, 1, b)
print("Number of rows in data.table")
print(length(regtools_data$samples))

paste_commas <- function(v){
  return(paste(v,collapse = ","))
}
regtools_data$norm_scores_variant <- unlist(lapply(regtools_data$norm_scores_variant,paste_commas))
regtools_data$norm_scores_non <- unlist(lapply(regtools_data$norm_scores_non,paste_commas))
columns_to_keep = c('samples', 'variant_info.x', 'gene_names', 'sample', "chrom.x", "start.x", "end.x", 'strand.x', 'anchor.x', 'info',
                    'names', 'mean_norm_score_variant', 'sd_norm_score_variant', 'norm_scores_variant',
                    'total_score_variant', 'mean_norm_score_non', 'sd_norm_score_non', 'norm_scores_non',
                    'total_score_non', 'p_value_mean','p_value_min')
regtools_data = subset(regtools_data, select=columns_to_keep)
colnames(regtools_data) <- c('variant_samples', 'variant_info', 'gene_names', 'junction_samples', "chrom", "start", "end", 'strand', 'anchor', 'variant_junction_info',
                             'names', 'mean_norm_score_variant', 'sd_norm_score_variant', 'norm_scores_variant',
                             'total_score_variant', 'mean_norm_score_non', 'sd_norm_score_non', 'norm_scores_non',
                             'total_score_non', 'p_value_mean','p_value_min')
regtools_data$sd_norm_score_variant[is.na(regtools_data$sd_norm_score_variant)] = 0
regtools_data$mean_norm_score_non[is.na(regtools_data$mean_norm_score_non)] = 0
regtools_data$sd_norm_score_non[is.na(regtools_data$sd_norm_score_non)] = 0
regtools_data$total_score_variant[is.na(regtools_data$total_score_variant)] = 0
regtools_data$total_score_non[is.na(regtools_data$total_score_non)] = 0
all_splicing_variants <- as.data.table(all_splicing_variants)
regtools_data = regtools_data %>% distinct()


write.table(regtools_data, file='output.tsv', quote = FALSE, sep = '\t', row.names = F)