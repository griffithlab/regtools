# compare_junctions_hist.R
# Rscript --vanilla compare_junctions_hist.R <tag>

# load libraries
library(data.table)
library(graphics)
library(plyr)
library(foreach)
library(doParallel)
registerDoParallel(cores=32)

debug = F

system.time({
if (debug){
  tag = paste("_", "default", sep="")
} else {
  # get options tag
  argc = length(commandArgs())
  tag = paste("_", commandArgs(trailingOnly = F)[argc], sep="")

  if ( substr(tag, 2, 3) == "--"){
    stop("Please specify an option tag (e.g. \"default\", \"i20e5\")")
  }
}

## All splicing relevant variants (union of rows from variants.bed files; add column with comma-separated list of sample names)
all_splicing_variants = unique(data.table::fread(paste("all_splicing_variants", tag, ".bed",sep=""), sep = '\t', header = T, stringsAsFactors = FALSE))
colnames(all_splicing_variants) <- c("chrom", "start", "end", "samples")

all_splicing_variants = as.data.table(aggregate(samples ~ chrom + start + end, paste, data=all_splicing_variants))
all_splicing_variants$key <- paste0(all_splicing_variants$chrom, ":", all_splicing_variants$start, "-", all_splicing_variants$end) #this key is just a 1bp-long chrom:start-stop designed to match the regtools output variant_info column

## Get all of the samples
all_samples = strsplit(scan("dir_names.tsv", what="", sep="\n"), "[[:space:]]+")

################################################################################
##### Helper functions #########################################################
################################################################################

# basically combines all the cse_identify_filtered_compare files for this cohort
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

# this is regtools compare output across samples (so it contains variant-junction lines even from samples without variant)
df <- rbindlist(lapply(all_samples, get_sample_data)) # eventually becomes all_cse_identify_data in function a
df$info <- paste(df$chrom, df$start, 
                 df$end, df$anchor, 
                 df$variant_info, sep = "_")

################################################################################
##### Function to create reviewable regtools output files ######################
################################################################################

a <- function(x, all_cse_identify_data){
  ## x comes from all_splicing_variants, so this just gives the samples with the variant
  variant_samples <- unlist(x$samples) 
  
  ##############################################################################
  ##### Get junction and variant information for variant samples ###############
  ##############################################################################
  
  ## Get junctions from concatenated regtools output associated with each variant, only if the data came from a sample actually containing the variant
  variant_junctions_data <- all_cse_identify_data[all_cse_identify_data$variant_info==x$key & all_cse_identify_data$sample %in% variant_samples,]
  
  ## Make dataset with the sample and the junction/variant information
  sample_data <- variant_junctions_data[,c("sample", "info")]
  
  # Normalize junction scores in this variant region for each sample
  summed_variant_scores = variant_junctions_data[, list(score=sum(score)), by=list(sample)]
  colnames(summed_variant_scores) <- c("sample", "score.tmp")
  variant_junctions_data <- merge(variant_junctions_data, summed_variant_scores, by="sample")
  variant_junctions_data$norm_score <- variant_junctions_data$score / variant_junctions_data$score.tmp
  
  # Aggregate data across junction-samples
  if (nrow(variant_junctions_data)[[1]] > 0){ # names column refers to names that regtools gave each junction in the variant output (only will be multiple if variant and junction was found in multiple samples)
    variant_junctions_aggr = variant_junctions_data[, list(names=list(name),mean_norm_score_variant=mean(norm_score),sd_norm_score_variant=sd(norm_score),norm_scores_variant=list(norm_score),total_score_variant=sum(score)), 
                                                    by=list(chrom,start,end,strand,anchor,variant_info,info)]
    
  } else {
    return(data.table())
  }
  
  ## Identify mutated samples
  a <- function(x, tempData){
    return(paste(as.character(tempData[info == x, "sample"]$sample), collapse=",")) # for super-junction aggregated across samples with the variant, get list of samples with variant (and junction)
  }
  variant_junctions_aggr$sample <- sapply(variant_junctions_aggr$info, a, tempData = sample_data)
  
  ##############################################################################
  ##### Get junction and variant information for non-variant samples ###########
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

  # Aggregate data across junction-samples
  if (nrow(variant_junctions_data)[[1]] > 0){
    non_variant_junctions_aggr = non_variant_junctions_data[, list(mean_norm_score_non=mean(norm_score),sd_norm_score_non=sd(norm_score),norm_scores_non=sd(norm_score),total_score_non=sum(score)), 
                                                    by=list(chrom,start,end,strand,anchor,variant_info,info)]
  } else {
    return(data.table())
  }
  
  # Generate histogram of normalized junction scores
  collapsed_norm_scores = non_variant_junctions_data[, list(norm_scores=list(norm_score)), by=list(chrom,start,end,strand)]
  
  # Add non-variant sample aggregate stats to output
  variant_junctions_aggr = merge(x=variant_junctions_aggr, y=non_variant_junctions_aggr,by= 
                c("chrom","start","end","strand","anchor","variant_info","info"),all.x=T)
  
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
    variant_norm_score = junction_data[["mean_norm_score_variant"]]
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

# for each iteration (i.e. line in all_splicing_variants),
# return a dataframe of regtools info and stats, basically what you end up seeing in regtools_data
# so the samples column comes from all_splicing_variants, since adply uses rbind.fill, which fills in missing columns
regtools_data <- adply(all_splicing_variants, 1, a, all_cse_identify_data=df, .parallel=!debug)

paste_commas <- function(v){
  return(paste(v,collapse = ","))
}

regtools_data <- regtools_data[order(anchor, sample, pvalue)]
regtools_data$norm_scores_variant <- unlist(lapply(regtools_data$norm_scores_variant,paste_commas))
regtools_data$norm_scores_non <- unlist(lapply(regtools_data$norm_scores_non,paste_commas))
regtools_data$samples <- unlist(lapply(regtools_data$samples,paste_commas))
regtools_data$names <- unlist(lapply(regtools_data$names,paste_commas))
regtools_data$sd_norm_score_variant[is.na(regtools_data$sd_norm_score_variant)] = 0
regtools_data$mean_norm_score_non[is.na(regtools_data$mean_norm_score_non)] = 0
regtools_data$sd_norm_score_non[is.na(regtools_data$sd_norm_score_non)] = 0
regtools_data$total_score_variant[is.na(regtools_data$total_score_variant)] = 0
regtools_data$total_score_non[is.na(regtools_data$total_score_non)] = 0

regtools_data = subset(regtools_data, select = !(names(regtools_data) %in% c("key")))
colnames(regtools_data)[colnames(regtools_data)=="sample"] <- "variant_junction_samples"
colnames(regtools_data)[colnames(regtools_data)=="samples"] <- "variant_samples"
colnames(regtools_data)[colnames(regtools_data)=="info"] <- "variant_junction_info"
})

write.table(regtools_data, file=paste("compare_junctions/hist/", "junction_pvalues", tag, ".tsv", sep=""), quote=FALSE, sep='\t', row.names = F)
