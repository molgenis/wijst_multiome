#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_filter_down_significant_results.R
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

#' get the number eGenes per cell type from QTL output
#' 
#' @param input_dir base location of the folder containing files to filter
#' @param input_file which output file to read for the results
#' @param output_dir base location of the folder to put filtered output
#' @param output_file_prepend prepend of filtered file to create
#' @param output_file_append append of filtered file to create
#' @param significance_column column denoting significance
#' @param significance_cutoff cutoff for which to set significance
#' @param split_column which column to split the output on
#' @param add_mtc add multiple testing before filtering down
#' @param mtc_column the column of values to apply multiple testing on
#' @param feature_mtc_column the column that has the feature group to perform the multiple testing on
#' @param mtc_column_to_add the name of the column that has the mtc-corrected values
#' @param filter_alpha remove entries that have an abhorrant alpha param
#' @param verbose print progress
#' @returns 0 if success
#' 
split_output_by_column <- function(input_dir, input_file='qtl_results_all.txt.gz', output_dir=NULL, output_file_prepend='qtl_results_all_qval_', output_file_append='.txt.gz', split_column='feature_chromosome', add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', verbose=T, filter_alpha=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(input_dir, full.names = F, recursive = F)
  # we will store the results in a list for now
  numbers_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(input_dir, '/', cell_type, '/', input_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- fread(full_cell_type_path, sep = '\t', header = T)
    
    # filter on alpha if requested
    if (filter_alpha) {
      cell_type_output <- cell_type_output[!(cell_type_output[['alpha_param']] > 5 | cell_type_output[['alpha_param']] < .2), ]
    }
    
    # get the features and the emperical p value
    if (add_mtc) {
      if (verbose) {
        print(paste('adding MTC to', full_cell_type_path))
      }
      # get just the two columns we care about
      cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column), with = F]
      # order by significance
      cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
      # remove where the feature is smaller than zero
      cell_type_output_features <- cell_type_output_features[!(cell_type_output_features[[mtc_column]] < 0), ]
      # keep only the first entry
      cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
      # set the values that are larger than 1, to be 1, problem with precision
      cell_type_output_features[cell_type_output_features[[mtc_column]] > 1, mtc_column] <- 1
      print(min(cell_type_output_features[[mtc_column]]))
      print(max(cell_type_output_features[[mtc_column]]))
      # add multiple testing correction
      cell_type_output_features[['qvalue']] <- qvalue(cell_type_output_features[[mtc_column]])$qvalues
      # now add back to the original table
      cell_type_output[[mtc_column_to_add]] <- cell_type_output_features[match(cell_type_output[[feature_mtc_column]], cell_type_output_features[[feature_mtc_column]]), 'qvalue']
    }
    
    # print progress if requested
    if (verbose) {
      print(paste('variant+phenotype number', nrow(cell_type_output)))
    }
    
    # get the output location
    output_directory <- input_dir
    # if supplied, set the output directory
    if (!is.null(output_dir)) {
      output_directory <- output_dir
    }
    
    # check each split value
    for (split_value in unique(cell_type_output[[split_column]])) {
      if (verbose) {
        print(paste('splitting', full_cell_type_path, 'and writing split', split_value))
      }
      # subset to that value
      cell_type_output_split <- cell_type_output[!is.na(cell_type_output[[split_column]]) & cell_type_output[[split_column]] == split_value, ]
      # make the full path
      full_output_loc <- paste(output_directory, '/', cell_type, '/', output_file_prepend, split_value, output_file_append, sep = '')
      # store the gz connection if we need it
      full_output_loc_wzip <- full_output_loc
      # gzip it if the extention ends on gz
      if (grepl('.gz$', full_output_loc)) {
        full_output_loc_wzip <- gzfile(full_output_loc)
      }
      # write result
      write.table(cell_type_output_split, full_output_loc_wzip, sep = '\t', row.names= F, col.names = T)
      # create md5
      mdfiver::create_md5_for_file(full_output_loc)
    }
    if (verbose) {
      print(paste('finished', full_cell_type_path))
    }
  }
  return(0)
}


#' get the number eGenes per cell type from QTL output
#' 
#' @param unfiltered_loc base location of the folder containing files to filter
#' @param unfiltered_file which output file to read for the results
#' @param filtered_loc base location of the folder to put filtered output
#' @param filtered_file name of filtered file to create
#' @param significance_column column denoting significance
#' @param significance_cutoff cutoff for which to set significance
#' @param verbose print progress
#' @param add_mtc add multiple testing before filtering down
#' @param mtc_column the column of values to apply multiple testing on
#' @param feature_mtc_column the column that has the feature group to perform the multiple testing on
#' @param mtc_column_to_add the name of the column that has the mtc-corrected values
#' @returns 0 if success
#' 
filter_output_by_significance <- function(unfiltered_loc, unfiltered_file='qtl_results_all.txt.gz', filtered_loc=NULL, filtered_file=NULL, significance_column='p_value', significance_cutoff=0.05, verbose=T, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value') {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(unfiltered_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  numbers_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(unfiltered_loc, '/', cell_type, '/', unfiltered_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- fread(full_cell_type_path, sep = '\t', header = T)
    
    # get the features and the emperical p value
    if (add_mtc) {
      # get just the two columns we care about
      cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column)]
      # order by significance
      cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
      # keep only the first entry
      cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
      # set the values that are larger than 1, to be 1, problem with precision
      cell_type_output_features[cell_type_output_features[[mtc_column]] > 1, mtc_column] <- 1
      # add multiple testing correction
      cell_type_output_features[['qvalue']] <- qvalue(cell_type_output_features[[mtc_column]])$qvalues
      # now add back to the original table
      cell_type_output[[mtc_column_to_add]] <- cell_type_output_features[match(cell_type_output[[feature_mtc_column]], cell_type_output_features[[feature_mtc_column]]), 'qvalue'][['qvalue']]
    }
    
    # print progress if requested
    if (verbose) {
      print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
    }
    # filter
    cell_type_output <- cell_type_output[
      !is.na(cell_type_output[[significance_column]]) &
        cell_type_output[[significance_column]] < significance_cutoff, 
    ]
    if (verbose) {
      print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
    }
    # get the output location
    output_dir <- unfiltered_loc
    # if supplied, set the output directory
    if (!is.null(filtered_loc)) {
      output_dir <- filtered_loc
    }
    # set the output file location
    output_file <- NULL
    # if supplied set that
    if (!is.null(filtered_file)) {
      output_file <- filtered_file
    }
    else {
      output_file <- paste('filtered', '.', unfiltered_file, sep = '')
    }
    # make the full path
    full_output_loc <- paste(output_dir, '/', cell_type, '/', output_file, sep = '')
    # store the gz connection if we need it
    full_output_loc_wzip <- full_output_loc
    # gzip it if the extention ends on gz
    if (grepl('.gz$', full_output_loc)) {
      full_output_loc_wzip <- gzfile(full_output_loc)
    }
    # write result
    write.table(cell_type_output, full_output_loc_wzip, sep = '\t', row.names= F, col.names = T)
    # create md5
    mdfiver::create_md5_for_file(full_output_loc)
  }
  # return 0 upon success
  return(0)
}


####################
# Main Code        #
####################

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/UT/'
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output_cis_50kb_updated_features/L1/UT/'
# perform splitting
split_output_by_column(
  input_dir=eqtl_output_loc,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T
)
# perform filtering
filter_output_by_significance(
  unfiltered_loc=eqtl_output_loc, 
  unfiltered_file='qtl_results_all.txt.gz', 
  filtered_loc=NULL, 
  filtered_file='qtl_results_all_nominally_significant.txt.gz', 
  significance_column='p_value', 
  significance_cutoff=0.05, 
  verbose=T
)
filter_output_by_significance(
  unfiltered_loc=caqtl_output_loc, 
  unfiltered_file='qtl_results_all.txt.gz', 
  filtered_loc=NULL, 
  filtered_file='qtl_results_all_nominally_significant.txt.gz', 
  significance_column='p_value', 
  significance_cutoff=0.05, 
  verbose=T
)
filter_output_by_significance(
  unfiltered_loc=eqtl_output_loc, 
  unfiltered_file='qtl_results_all_nominally_significant.txt.gz', 
  filtered_loc=NULL, 
  filtered_file='qtl_results_all_fdr01_significant.txt.gz', 
  significance_column='feature_q_value', 
  significance_cutoff=0.1, 
  verbose=T, 
  add_mtc = F
)
filter_output_by_significance(
  unfiltered_loc=caqtl_output_loc, 
  unfiltered_file='qtl_results_all_nominally_significant.txt.gz', 
  filtered_loc=NULL, 
  filtered_file='qtl_results_all_fdr01_significant.txt.gz', 
  significance_column='feature_q_value', 
  significance_cutoff=0.1, 
  verbose=T, 
  add_mtc = F
)
