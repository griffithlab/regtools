# compare_junctions_hist.R
# Rscript --vanilla compare_junctions_hist.R <tag>

# set dir
setwd("~/Desktop/regtools/CHOL")

# load libraries
library(data.table)
library(profvis)



debug = F

# system.time({
# if (debug){
#   tag = paste("_", "default", sep="")
# } else {
#   # get options tag
#   argc = length(commandArgs())
#   tag = paste("_", commandArgs(trailingOnly = F)[argc], sep="")
# 
#   if ( substr(tag, 2, 3) == "--"){
#     stop("Please specify an option tag (e.g. \"default\", \"i20e5\")")
#   }
# }

tag <- "default"
sample <- "TCGA-3X-AAV9-01A"

## All splicing relevant variants (union of rows from variants.bed files; add column with comma-separated list of sample names)
all_splicing_variants = unique(data.table::fread(paste("all_splicing_variants_", tag, ".bed",sep=""), sep = '\t', header = T, stringsAsFactors = FALSE))
colnames(all_splicing_variants) <- c("chrom", "start", "end", "samples")

# all_splicing_variants = as.data.table(aggregate(samples ~ chrom + start + end, paste, data=all_splicing_variants)) # I don't think this is doing anything?
all_splicing_variants$key <- paste0(all_splicing_variants$chrom, ":", all_splicing_variants$start, "-", all_splicing_variants$end) #this key is just a 1bp-long chrom:start-stop designed to match the regtools output variant_info column

## Get all of the samples
all_samples = strsplit(scan("dir_names.tsv", what="", sep="\n"), "[[:space:]]+")

################################################################################
##### Helper functions #########################################################
################################################################################

# basically combines all the cse_identify_filtered_compare files for this cohort
get_sample_data <- function(sample, all_splicing_variants){
  file_name = paste("samples/", sample, "/output/cse_identify_filtered_compare_", tag,".tsv", sep = "")
  cse_identify_data = data.table::fread(file_name, sep = '\t', header = T, stringsAsFactors = FALSE, strip.white = TRUE)
  cse_identify_data$sample <- sample
  cse_identify_data <- cse_identify_data[,.(sample,variant_info,chrom,start,end,strand,anchor,score,name)]
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
cse_identify_v1[,score_norm := score/score.tmp]
cse_identify_v1[,mean_norm_score_variant := mean(score_norm), by=.(sample,variant_info,chrom,start,end,strand,anchor,info)]
cse_identify_v1[,sd_norm_score_variant := sd(score_norm), .(sample,variant_info,chrom,start,end,strand,anchor,info)]
cse_identify_v1[,total_score_variant := sum(score), .(sample,variant_info,chrom,start,end,strand,anchor,info)]
cse_identify_v1 <- cse_identify_v1

print("test4")

# subset and rename columns to match the original output
cse_identify_v1 <- cse_identify_v1[,c("sample", "variant_info", "chrom", "start", "end", "strand", "anchor",
                                      "variant_info", "info", "name", "mean_norm_score_variant",
                                      "sd_norm_score_variant", "score_norm", "total_score_variant")]
colnames(cse_identify_v1) <- c("samples", "key", "chrom", "start", "end", "strand", "anchor", "variant_info",
                               "info", "names", "mean_norm_score_variant", "sd_norm_score_variant",
                               "norm_scores_variant", "total_score_variant")

################ aggrregate variants with no sample ############################

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

# zeroes need to be added in for some samples
a <- function(x, y){
  toAdd <- y - length(x)
  toAdd <- rep(0.0000000, toAdd)
  x <- c(x, toAdd)
  return(x)
}
regtools_data$norm_scores_non <- lapply(regtools_data$norm_scores_non, a, length(all_samples))
print("test7")
################ calculate p-values ############################################

a <- function(x){
  variant_norm_score = as.numeric(x$norm_scores_variant)
  if(length(x$norm_scores_non) <= 1){
    return(0)
  } 
  all_norm_scores = c(x$norm_scores_non, variant_norm_score)
  countable = rank(all_norm_scores)
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

regtools_data$p_value <- apply(regtools_data, 1, a)
print("test8")


paste_commas <- function(v){
   return(paste(v,collapse = ","))
}

regtools_data$norm_scores_variant <- unlist(lapply(regtools_data$norm_scores_variant,paste_commas))
regtools_data$norm_scores_non <- unlist(lapply(regtools_data$norm_scores_non,paste_commas))
write.table(regtools_data, file=paste("compare_junctions/hist/", "junction_pvalues", tag, ".tsv", sep=""), quote=FALSE, sep='\t', row.names = F)


