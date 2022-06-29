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


read_file='input_file'
regtools_data = unique(data.table::fread(file=read_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE))
regtools_data_filtered = regtools_data[(regtools_data$p_value_mean >= 0 & 
											(regtools_data$anchor == "D" | 
												regtools_data$anchor == "A" | 
												regtools_data$anchor == "NDA"))]

p = regtools_data_filtered$p_value_mean
adjusted_p = p.adjust(p, method = "BH")
regtools_data_filtered$adjusted_p = adjusted_p
regtools_data_filtered_sorted = regtools_data_filtered[order(adjusted_p)]

write_file = '/Volumes/mgriffit/Active/regtools/scrna/all_cells_3rd_variant_samples_rm_w_cellcounts_zeros_rm_filtered_D_A_NDA.tsv'
write.table(regtools_data_filtered_sorted, file=write_file, quote=FALSE, sep='\t', row.names = FALSE)

threshold = 0.05
is_significant = regtools_data_filtered_sorted$adjusted_p < threshold
regtools_data_significant_filtered_sorted = regtools_data_filtered_sorted[is_significant] 

write_file = '/Volumes/mgriffit/Active/regtools/scrna/all_cells_3rd_variant_samples_rm_w_cellcounts_zeros_rm_filtered_sig_D_A_NDA.tsv'
write.table(regtools_data_significant_filtered_sorted, file=write_file, quote=FALSE, sep='\t', row.names = FALSE)
