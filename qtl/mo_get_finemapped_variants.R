#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_finemapped_variants.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(mdfiver)
library(qvalue)
library(data.table)

####################
# Functions        #
####################

merge_finemapped_results <- function(finemap_loc, finemap_prepend, finemap_append, out_loc) {
  # list all the files
  finemap_files <- list.files(finemap_loc)
  # filter on the ones that match the pattern
  finemap_files <- finemap_files[grepl(paste(finemap_prepend, '\\d+', finemap_append, sep = ''), finemap_files)]
  # now read each one into a list
  finemap_files_list <- list()
  for (finemap_file in finemap_files) {
    # get the full location
    full_file_loc <- paste(finemap_loc, '/', finemap_file, sep = '')
    # read the file
    finemapping_chrom <- fread(full_file_loc, header = T, sep = '\t')
    # put into the list
    finemap_files_list[[finemap_file]] <- finemapping_chrom
  }
  # merge them all
  finemap_files_merged <- do.call('rbind', finemap_files_list)
  # store the gz connection if we need it
  full_output_loc_wzip <- out_loc
  # gzip it if the extention ends on gz
  if (grepl('.gz$', out_loc)) {
    full_output_loc_wzip <- gzfile(out_loc)
  }
  # write result
  write.table(finemap_files_merged, full_output_loc_wzip, sep = '\t', row.names= F, col.names = T)
  # create md5
  mdfiver::create_md5_for_file(out_loc)
}


get_finemapped_variants_per_feature_hardcutoff <- function(finemapped_file_loc, feature_column='feature_id', pip_column='SusieRss_pip', pip_cutoff=0.8) {
  # read the file
  finemapping <- fread(finemapped_file_loc, header = T, sep = '\t')
  # we'll store the rows first
  finemapped_rows <- list()
  # check each gene
  for (feature in unique(finemapping[[feature_column]])) {
    # get the rows that have a pip of higher than what was supplied
    finemapped_rows_feature <- finemapping[
      !is.na(finemapping[[feature_column]]) &
        finemapping[[feature_column]] == feature &
        !is.na(finemapping[[pip_column]]) &
        finemapping[[pip_column]] > pip_cutoff
    ]
    # put in a list if there is at least one finemapped variant
    if (nrow(finemapped_rows_feature) > 0) {
      finemapped_rows[[feature]] <- finemapped_rows_feature
    }
  }
  # merge all the finemapped rows
  finemapped_all <- do.call('rbind', finemapped_rows)
  return(finemapped_all)
}

####################
# Main Code        #
####################

# merge the finemapped results
merge_finemapped_results(
  finemap_loc = '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_14/Finemapping/', 
  finemap_prepend = '062024_susie_finemap_Mono_', 
  finemap_append = '_formatted.txt', 
  out_loc = '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_14/Finemapping/finemapped_monocyte.tsv.gz'
)
# get the finemapped variants
finemapped_eqtls_monocyte <- get_finemapped_variants_per_feature_hardcutoff('/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_14/Finemapping/finemapped_monocyte.tsv.gz')
# remove where there is no credible set defined
finemapped_eqtls_monocyte <- finemapped_eqtls_monocyte[!is.na(finemapped_eqtls_monocyte[['SusieRss_CS']])]
# write the result
write.table(data.frame(coQTL = paste(finemapped_eqtls_monocyte[['snp_id']], finemapped_eqtls_monocyte[['feature_id']], sep = ';')), gzfile('/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_14/Finemapping/finemapped_monocyte_enrichment_input.tsv.gz'), row.names = F, col.names = T, quote = F, sep = '\t')
