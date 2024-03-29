# compare_junctions_hist.R
# Rscript --vanilla compare_junctions_hist.R <tag>


# load libraries
library(data.table)
library(plyr)
library(tidyverse)

debug = F

# system.time({
# if (debug){
#   tag = paste("_", "default", sep="")
# } else {
#   # get options tag
#   args = commandArgs(trailingOnly = TRUE)
#   tag = args[1]
#   input_file = args[2]
#   if ( substr(tag, 2, 3) == "--"){
#     stop("Please specify an option tag (e.g. \"default\", \"i20e5\")")
#   }
# }

setwd('~/Desktop/CHOL')
tag = 'E'
input_file = '/Users/kcotto/Desktop/CHOL/all_splicing_variants_E.bed'

# All splicing relevant variants (union of rows from variants.bed files; add column with comma-separated list of sample names)
all_splicing_variants = unique(data.table::fread(input_file), sep = '\t', header = T, stringsAsFactors = FALSE)
colnames(all_splicing_variants) <- c("chrom", "start", "end", "samples")

# all_splicing_variants = as.data.table(aggregate(samples ~ chrom + start + end, paste, data=all_splicing_variants)) # I don't think this is doing anything?
all_splicing_variants$key <- paste0(all_splicing_variants$chrom, ":", all_splicing_variants$start, "-", all_splicing_variants$end) #this key is just a 1bp-long chrom:start-stop designed to match the regtools output variant_info column

## Get all of the samples
all_samples = strsplit(scan("/Users/kcotto/Desktop/CHOL/dir_names.tsv", what="", sep="\n"), "[[:space:]]+")

################################################################################
##### Helper functions #########################################################
################################################################################

# basically combines all the cse_identify_filtered_compare files for this cohort
get_sample_data <- function(sample, all_splicing_variants){
  file_name = paste("samples/", sample, "/output/cse_identify_filtered_compare_", tag,".tsv", sep = "")
  cse_identify_data = data.table::fread(file_name, sep = '\t', header = T, stringsAsFactors = FALSE, strip.white = TRUE)
  cse_identify_data$sample <- sample
  cse_identify_data <- cse_identify_data[,.(sample,variant_info,chrom,start,end,strand,anchor,score,name,genes)]
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
dt <- rbindlist(lapply(all_samples, get_sample_data, all_splicing_variants=all_splicing_variants))

dt[,info := paste(chrom, start, end, anchor, variant_info, sep="_")]

################ subset out entries where we can for now #######################

# make a sample/variant_info key
dt[,key := paste0(variant_info, "_", sample)]
dt[,junction := paste0(chrom, ":", start, "-", end, "_", sample)]

print("zl")
cse_identify_v1 <- dt
rm(dt)
print("kc")
cse_identify_v2 <- cse_identify_v1
print("yy")



######### aggregate scores for variant/samples identified in cse_identify ######

# key files with variant and sample keys
all_splicing_variants$key2 <- paste0(all_splicing_variants$key, "_",all_splicing_variants$samples)

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
                                      "variant_info", "info", "genes","name.y", "mean_norm_score_variant.y",
                                      "sd_norm_score_variant", "norm_scores_variant", "total_score_variant", "junction")]
colnames(cse_identify_v1) <- c("sample", "key", "chrom", "start", "end", "strand", "anchor", "variant_info",
                               "info", "genes", "names", "mean_norm_score_variant", "sd_norm_score_variant",
                               "norm_scores_variant", "total_score_variant", "junction_key")

################ aggregate variants with no sample ############################

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
#this is where non-variant samples are aggregated
cse_identify_v2 <- rbindlist(cse_identify_v2)

print("test6")

################ Merge the two DT's we've been working on together #############

regtools_data <- merge(cse_identify_v1, cse_identify_v2, by="info", all.x = TRUE)
rm(cse_identify_v1)
rm(cse_identify_v2)

tmp <- all_splicing_variants[order(key)]
tmp <- ddply(tmp, .(chrom, start, end, key), summarise, samples=paste(samples,collapse = ','))
regtools_data = merge(regtools_data, tmp, by.x=c('variant_info.x'), by.y=c('key'))
columns_to_keep = c('samples', 'variant_info.x', 'sample', "info", "chrom.x", "start.x", "end.x", 'strand.x', 'anchor.x', 'genes',
                    'names', 'mean_norm_score_variant', 'sd_norm_score_variant', 'norm_scores_variant',
                    'total_score_variant', 'mean_norm_score_non', 'sd_norm_score_non', 'norm_scores_non',
                    'total_score_non')
regtools_data = subset(regtools_data, select=columns_to_keep)


# zeroes need to be added in for some samples
a <- function(x, y){
  toAdd <- y - length(x) - 1
  # browser()
  toAdd <- rep(0.0000000, toAdd)
  x <- c(x, toAdd)
  return(x)
}
x <- mapply(a, regtools_data$norm_scores_non, length(all_samples))

get_num_zeros_to_rm <- function(z){
  num_zeroes_to_rm = str_count(z, ',') 
  return(num_zeroes_to_rm)
}

num_zeroes_to_rm <- mapply(get_num_zeros_to_rm, regtools_data$samples)

x = split(x, rep(1:ncol(x), each = nrow(x)))
regtools_data$norm_scores_non = x
regtools_data$zeroes_to_rm = num_zeroes_to_rm

rm_zeroes <- function(x,y){
  new_length <- length(x) - y
  x <- sort(x,decreasing = TRUE)
  x <- x[1:new_length]
  return(x)
}

if (max(num_zeroes_to_rm > 0)) {
x <- mapply(rm_zeroes, regtools_data$norm_scores_non, regtools_data$zeroes_to_rm)
regtools_data$norm_scores_non = x
}

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
  toAdd <- (str_count(y, ',') + 1) - (str_count(x, ',') + 1) 
  # browser()
  if (toAdd > 0) {
    toAdd <- rep(0.0000000, toAdd)
    x <- c(x, toAdd)
  } else {
    x <- unlist(strsplit(x, ","))
  }
  x <- list(x)
  return(x)
}
x <- mapply(a, regtools_data$norm_scores_variant, regtools_data$samples)
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
columns_to_keep = c('samples', 'variant_info.x', 'genes', 'sample', "chrom.x", "start.x", "end.x", 'strand.x', 'anchor.x', 'info',
                    'names', 'mean_norm_score_variant', 'sd_norm_score_variant', 'norm_scores_variant',
                    'total_score_variant', 'mean_norm_score_non', 'sd_norm_score_non', 'norm_scores_non',
                    'total_score_non', 'p_value_mean','p_value_min')
regtools_data = subset(regtools_data, select=columns_to_keep)
colnames(regtools_data) <- c('variant_samples', 'variant_info', 'genes', 'junction_samples', "chrom", "start", "end", 'strand', 'anchor', 'variant_junction_info',
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


write.table(regtools_data, file=paste(input_file, "_out.tsv", sep=""), quote=FALSE, sep='\t', row.names = F)

# })
