# compare_junctions_hist.R
# Rscript --vanilla compare_junctions_hist.R <tag>

# load libraries
library(data.table)
library(graphics)
library(plyr)
library(foreach)
library(doParallel)
registerDoParallel(cores=32)

system.time({
# get options tag
#argc = length(commandArgs())
tag = paste("_", "default", sep="")

if ( substr(tag, 2, 3) == "--"){
  stop("Please specify an option tag (e.g. \"default\", \"i20e5\")")
}

## All splicing relevant variants (union of rows from variants.bed files; add column with comma-separated list of sample names)
all_splicing_variants = unique(data.table::fread(paste("all_splicing_variants", tag, ".bed",sep=""), sep = '\t', header = T, stringsAsFactors = FALSE))
colnames(all_splicing_variants) <- c("chrom", "start", "end", "samples")

all_splicing_variants = as.data.table(aggregate(samples ~ chrom + start + end, paste, data=all_splicing_variants))
all_splicing_variants$key <- paste0(all_splicing_variants$chrom, ":", all_splicing_variants$start, "-", all_splicing_variants$end)

## Get all of the samples 
all_samples = strsplit(scan("dir_names.tsv", what="", sep="\n"), "[[:space:]]+")

################################################################################
##### Helper functions #########################################################
################################################################################

get_sample_data <- function(sample){
  
  file_name = paste("samples/", sample, "/output/cse_identify_filtered_compare", tag,".tsv", sep = "")
  cse_identify_data = data.table::fread(file_name, sep = '\t', header = T, stringsAsFactors = FALSE, strip.white = TRUE)
  cse_identify_data$sample <- sample
  s <- strsplit(cse_identify_data$variant_info, split = ",")
  cse_identify_data <- data.table::data.table(chrom=rep(cse_identify_data$chrom, sapply(s, length)),
                                              start=rep(cse_identify_data$start, sapply(s, length)),
                                              end=rep(cse_identify_data$end, sapply(s, length)),
                                              name=rep(cse_identify_data$name, sapply(s, length)),
                                              score=rep(cse_identify_data$score, sapply(s, length)),
                                              strand=rep(cse_identify_data$strand, sapply(s, length)),
                                              splice_site=rep(cse_identify_data$splice_site, sapply(s, length)),
                                              acceptors_skipped=rep(cse_identify_data$acceptors_skipped, sapply(s, length)),
                                              exons_skipped=rep(cse_identify_data$exons_skipped, sapply(s, length)),
                                              donors_skipped=rep(cse_identify_data$donors_skipped, sapply(s, length)),
                                              anchor=rep(cse_identify_data$anchor, sapply(s, length)),
                                              known_donor=rep(cse_identify_data$known_donor, sapply(s, length)),
                                              known_acceptor=rep(cse_identify_data$known_acceptor, sapply(s, length)),
                                              known_junction=rep(cse_identify_data$known_junction, sapply(s, length)),
                                              genes=rep(cse_identify_data$genes, sapply(s, length)),
                                              transcripts=rep(cse_identify_data$transcripts, sapply(s, length)),
                                              sample=rep(cse_identify_data$sample, sapply(s, length)),
                                              variant_info=unlist(s))
  return(cse_identify_data)
}

################################################################################
##### Combine cse_identify_filtered_compare.tsv files for each sample ##########
################################################################################

# this is regtools output per sample
df <- rbindlist(lapply(all_samples, get_sample_data))
df$info <- paste(df$chrom, df$start, 
                 df$end, df$anchor, 
                 df$variant_info, sep = "_")

################################################################################
##### Function to create reviewable regtools output files ######################
################################################################################

a <- function(x, all_cse_identify_data){
  
  ## Get data from samples with variant
  variant_samples <- unlist(x$samples)
  
  ##############################################################################
  ##### Get junction and variant information for variant samples ###############
  ##############################################################################
  
  ## Get the data from sample with the variant and the significant junctions
  variant_junctions_data <- all_cse_identify_data[all_cse_identify_data$variant_info==x$key & all_cse_identify_data$sample %in% variant_samples,]
  
  ## Make dataset with the sample and the junction/variant information
  sample_data <- variant_junctions_data[,c("sample", "info")]
  
  # Normalize junction scores in this variant region for each sample
  summed_variant_scores = variant_junctions_data[, list(score=sum(score)), by=list(sample)]
  colnames(summed_variant_scores) <- c("sample", "score.tmp")
  variant_junctions_data <- merge(variant_junctions_data, summed_variant_scores, by="sample")
  variant_junctions_data$norm_score <- variant_junctions_data$score / variant_junctions_data$score.tmp
  
  # Aggregate data across junction-samples
  if (nrow(variant_junctions_data)[[1]] > 0){
    variant_junctions_aggr = variant_junctions_data[, list(average_norm_score=mean(norm_score),total_score=sum(score)), 
                                                    by=list(chrom,start,end,name,strand,anchor,variant_info,info)]
    
  } else {
    return(data.table())
  }
  
  ## Identify mutated samples
  a <- function(x, tempData){
    return(paste(as.character(tempData[info == x, "sample"]$sample), collapse="|"))
  }
  variant_junctions_aggr$sample <- sapply(variant_junctions_aggr$info, a, tempData = sample_data)
  
  ##############################################################################
  ##### Get junction and variant information for non-variant samples###############
  ##############################################################################
  
  ## Samples without variant
  non_samples = setdiff(all_samples, variant_samples)
  non_variant_junctions_data <- all_cse_identify_data[all_cse_identify_data$variant_info==x$key & all_cse_identify_data$sample %in% non_samples]
  
  ## Get variant info
  non_variant_junctions_data$info <- paste(non_variant_junctions_data$chrom, non_variant_junctions_data$start, 
                                           non_variant_junctions_data$end, non_variant_junctions_data$anchor, 
                                           non_variant_junctions_data$variant_info, sep = "_")
  
  ## Make dataset with the sample and the junction/variant information
  non_variant_sample_data <- non_variant_junctions_data[,c("sample", "info")]
  
  # Normalize junction scores in this variant region for each sample
  summed_non_variant_scores = non_variant_junctions_data[, list(score=sum(score)), by=list(sample)]
  colnames(summed_non_variant_scores) <- c("sample", "score.tmp")
  non_variant_junctions_data <- merge(non_variant_junctions_data, summed_non_variant_scores, by="sample")
  non_variant_junctions_data$norm_score <- non_variant_junctions_data$score / non_variant_junctions_data$score.tmp
  
  

  
  # Generate histogram of normalized junction scores
  collapsed_norm_scores = non_variant_junctions_data[, list(norm_scores=list(norm_score)), by=list(chrom,start,end,strand)]
  
  # split out what p-values can be calculated on vs not
  collapsed_norm_scores$key <- paste(collapsed_norm_scores$chrom, collapsed_norm_scores$start, collapsed_norm_scores$end, collapsed_norm_scores$strand, sep='.')
  variant_junctions_aggr$key <- paste(variant_junctions_aggr$chrom, variant_junctions_aggr$start, variant_junctions_aggr$end, variant_junctions_aggr$strand, sep='.')
  variant_junctions_aggr_run_pvalues <- variant_junctions_aggr[variant_junctions_aggr$key %in% collapsed_norm_scores$key,]
  variant_junctions_aggr_norun_pvalues <- variant_junctions_aggr[!variant_junctions_aggr$key %in% collapsed_norm_scores$key,]
  variant_junctions_aggr_norun_pvalues$pvalue <- 0
  
  calculate_pvalues <- function(junction_data){ #### pull out
    chrom = trimws(junction_data[["chrom"]])
    start = trimws(junction_data[["start"]])
    end = trimws(junction_data[["end"]])
    strand = trimws(junction_data[["strand"]])
    variant_norm_score = junction_data[["average_norm_score"]]
    ### pass this variable
    i = which(collapsed_norm_scores$chrom == chrom & collapsed_norm_scores$start == start & collapsed_norm_scores$end == end & collapsed_norm_scores$strand == strand)  
    if (length(i)){
      non_variant_norm_scores = collapsed_norm_scores[i]$norm_scores[[1]]
      missing = length(non_samples) - length(non_variant_norm_scores)
      missing_zeros = integer(missing)
      non_variant_norm_scores = c(non_variant_norm_scores, missing_zeros)
      
      variant_norm_score = as.numeric(variant_norm_score)
      all_norm_scores = c(non_variant_norm_scores, variant_norm_score)
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
    } else {
      return(0)
    }
  }
  if(nrow(variant_junctions_aggr_run_pvalues) > 0){
    # apply above function to each aggregated variant junction
    pvalues = apply(variant_junctions_aggr_run_pvalues, 1, FUN = calculate_pvalues)
    variant_junctions_aggr_run_pvalues$pvalue = pvalues
    
    variant_junctions_aggr <- rbind(variant_junctions_aggr_norun_pvalues, variant_junctions_aggr_run_pvalues)
  } else {
    variant_junctions_aggr <- variant_junctions_aggr_norun_pvalues
  }
  
  return(variant_junctions_aggr)
  
}

regtools_data <- adply(all_splicing_variants, 1, a, all_cse_identify_data=df, .parallel=TRUE)

#regtools_data <- rbindlist(regtools_data, fill=TRUE)

regtools_data <- regtools_data[order(anchor, sample, pvalue)]
regtools_data$samples <- paste(regtools_data$samples)
})
write.table(regtools_data, file=paste("compare_junctions/hist/", "junction_pvalues", tag, ".tsv", sep=""), quote=FALSE, sep='\t', row.names = F)

