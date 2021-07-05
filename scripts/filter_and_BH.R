# filter_and_BH.R
library(data.table)
library(stats)

debug = F

if (debug){
  tag = "_default"
} else {
  # get options tag
  argc = length(commandArgs())
  tag = paste("_", commandArgs(trailingOnly = F)[argc], sep="")
  
  if ( substr(tag, 2, 3) == "--"){
    stop("Please specify an option tag (e.g. \"default\", \"i20e5\")")
  }
}


read_file=paste("compare_junctions/hist/", "junction_pvalues", tag, ".tsv", sep="")
regtools_data = unique(data.table::fread(file=read_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE))
regtools_data_filtered = regtools_data[(regtools_data$total_score_variant > 5 & 
											regtools_data$p_value_mean >= 0 & 
											(regtools_data$anchor == "DA"))]

p = regtools_data_filtered$p_value_mean
adjusted_p = p.adjust(p, method = "BH")
regtools_data_filtered$adjusted_p = adjusted_p
regtools_data_filtered_sorted = regtools_data_filtered[order(adjusted_p)]

write_file = paste("compare_junctions/hist/", "junction_pvalues_filtered_BH_DA_junctions", tag, ".tsv", sep="")
write.table(regtools_data_filtered_sorted, file=write_file, quote=FALSE, sep='\t', row.names = FALSE)

threshold = 0.05
is_significant = regtools_data_filtered_sorted$adjusted_p < threshold
regtools_data_significant_filtered_sorted = regtools_data_filtered_sorted[is_significant] 

write_file = paste("compare_junctions/hist/", "junction_pvalues_significant_",threshold,"_filtered_BH_DA_junctions", tag, ".tsv", sep="")
write.table(regtools_data_significant_filtered_sorted, file=write_file, quote=FALSE, sep='\t', row.names = FALSE)
