#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_caqtl_feature_filters.R
# Function: create feature filter files for caQTL mapping, based on percentage of cells having a peak as described in beds
############################################################################################################################


####################
# libraries        #
####################

# use this to load large tables
library(data.table)
# make md5 checksums
library(mdfiver)


####################
# Functions        #
####################

#' Convert Peak Files to Feature Lists
#'
#' This function reads peak files for specified cell types and conditions, merges the conditions for a cell type, filters the data based on a given column and value, and writes unique features to output files.
#'
#' @param peak_input_dir Character. Directory containing the peak files in bed format.
#' @param output_dir Character. Directory to save the output feature lists.
#' @param cell_types Character vector. Cell types to process. Default is c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK').
#' @param conditions Character vector. Conditions to process. Default is c('UT', '24hCA').
#' @param file_prepend Character. Prefix for the input peak file names. Default is 'mo_peaks_lane1to80_'.
#' @param file_append Character. Suffix for the input peak file names. Default is '.bed'.
#' @param output_prepend Character. Prefix for the output feature file names. Default is 'test_features_'.
#' @param output_append Character. Suffix for the output feature file names. Default is '.txt.gz'.
#' @param filter_column Character. Column name to apply the filter on. Default is 'pct_exp'.
#' @param filter_value Numeric. Minimum value for filtering the data. Default is 0.001.
#' @param feature_column Character. Column name to extract unique features. Default is 'name'.
#' @return Integer. Returns 0 upon successful completion.
#' @examples
#' \dontrun{
#' peak_files_to_feature_lists("path/to/peak_input", "path/to/output")
#' }
peak_files_to_feature_lists <- function(peak_input_dir, output_dir, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', '24hCA'), file_prepend='mo_peaks_lane1to80_', file_append='.bed', output_prepend='test_features_', output_append='.txt.gz', filter_column='pct_exp', filter_value=0.001, feature_column='name') {
  # check cell type
  for (cell_type in cell_types) {
    # we'll store each condition
    results_conditions_list <- list()
    # now check each condition
    for (condition in conditions) {
      # paste together the full path
      full_path_condition_celltype <- paste(peak_input_dir, '/', file_prepend, condition, '_', cell_type, file_append, sep = '')
      # check if this file exists
      if (file.exists(full_path_condition_celltype)) {
        # read the file
        condition_celltype <- fread(full_path_condition_celltype, header = T, sep = '\t', check.names = F)
        # put into the list
        results_conditions_list[[condition]] <- condition_celltype
      }
      else {
        warning(paste('cell type', cell_type, 'and condition', condition, 'does not have a file at', full_path_condition_celltype, 'and is skipped'))
      }
    }
    # merge conditions
    results_conditions <- do.call('rbind', results_conditions_list)
    # filter on minimum
    results_conditions <- results_conditions[results_conditions[[filter_column]] >= filter_value, ]
    # now get those features
    features_to_write <- unique(results_conditions[[feature_column]])
    # prepare a file to write
    features_to_write_loc <- paste(output_dir, '/', output_prepend, cell_type, output_append, sep = '')
    # gz file ends with .gz
    features_to_write_loc_gz <- features_to_write_loc
    if (grepl('.gz$', features_to_write_loc)) {
      features_to_write_loc_gz <- gzfile(features_to_write_loc)
    }
    # write result
    write.table(data.frame(x = features_to_write), features_to_write_loc_gz, row.names = F, col.names = F, quote = F)
    # and make an md5
    mdfiver::create_md5_for_file(features_to_write_loc)
  }
  return(0)
}



####################
# Settings        #
####################


####################
# Main Code        #
####################

# where the peaks are
peak_pct_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/'
# where the output should be placed
peak_confinement_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/confinements/pc01/'

peak_files_to_feature_lists(
  peak_pct_dir, 
  peak_confinement_dir
)
