#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_filter_down_significant_results.R
# Function: add multiple testing correction per feature to qtl outputs, and split them per chromosome
############################################################################################################################

####################
# libraries        #
####################

library(mdfiver)
library(qvalue)
library(data.table)
library(stats)


####################
# Functions        #
####################


#' Calculate Nominal Thresholds
#'
#' This function calculates nominal thresholds for p-values based on a given false discovery rate (FDR).
#'
#' @param res_df A data frame containing the results with p-values and other relevant columns.
#' @param fdr A numeric value specifying the false discovery rate threshold. Default is 0.05.
#' @param pval_col A character string specifying the name of the column with p-values. Default is 'p_value'.
#' @param nominal_threshold_column A character string specifying the name of the column to store the nominal thresholds. Default is 'pval_nominal_threshold'.
#' @param cutoff_column A character string specifying the name of the column with feature q-values. Default is 'feature_q_value'.
#' @param alpha_column A character string specifying the name of the column with alpha parameters for the beta distribution. Default is 'alpha_param'.
#' @param beta_column A character string specifying the name of the column with beta parameters for the beta distribution. Default is 'beta_param'.
#'
#' @return A data frame with an additional column for nominal thresholds.
#' @export
#'
#' @examples
#' \dontrun{
#'   res_df <- data.frame(
#'     p_value = runif(100),
#'     feature_q_value = runif(100),
#'     alpha_param = rep(1, 100),
#'     beta_param = rep(1, 100)
#'   )
#'   calculate_nominal_thresholds(res_df)
#' }
calculate_nominal_thresholds <- function(res_df, fdr=0.05, pval_col='p_value', nominal_threshold_column='pval_nominal_threshold', cutoff_column='feature_q_value', alpha_column='alpha_param', beta_column='beta_param') {
  # get the lowerbound p values, so the ones that are smaller than the FDR
  indices_lb <- res_df[[cutoff_column]] < fdr
  lb <- as.vector(res_df[indices_lb, ][[pval_col]])
  # put then in ascending order
  lb <- lb[order(lb)]
  # get the upperbound p values, so the ones that are bigger than the FDR
  indices_ub <- res_df[[cutoff_column]] > fdr
  ub <- as.vector(res_df[indices_ub, ][[pval_col]])
  # and order them
  ub <- ub[order(ub)]
  # if we have any significant effects, we can get a cutoff
  if (length(lb) > 0) {
    # get the highest (p) significant value
    highest_in_lb <- tail(lb, 1)
    # if there are any non significant effects
    if (length(ub) > 0) {
      # get the lowest (p) non-significant value
      lowest_in_ub <- head(ub, 1)
      # and calculate the threshold
      pthreshold <- (highest_in_lb + lowest_in_ub) / 2
    } else {
      # otherwise the highest effect will just be the cutoff
      pthreshold <- highest_in_lb
    }
    # ge the threshold, based on the shapes of the beta distribution and the significance threshold
    res_df[[nominal_threshold_column]] <- stats::qbeta(pthreshold, as.vector(res_df[[alpha_column]]), as.vector(res_df[[beta_column]]))
  }
  else {
    # otherwise it would have to be zero
    res_df[[nominal_threshold_column]] <- 0
  }
  return(res_df)
}


#' Filter File by Significance
#'
#' This function filters a file based on significance levels and optionally adds multiple testing correction (MTC) and global nominal thresholds.
#'
#' @param input_loc A character string specifying the location of the input file.
#' @param output_loc A character string specifying the location to save the filtered output file.
#' @param significance_column A character string specifying the name of the column with significance values. Default is 'p_value'.
#' @param significance_cutoff A numeric value specifying the significance cutoff threshold. Default is 0.05.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#' @param add_mtc A logical value indicating whether to add multiple testing correction. Default is TRUE.
#' @param mtc_column A character string specifying the name of the column with empirical feature p-values. Default is 'empirical_feature_p_value'.
#' @param feature_mtc_column A character string specifying the name of the column with feature IDs for MTC. Default is 'feature_id'.
#' @param mtc_column_to_add A character string specifying the name of the column to store the MTC values. Default is 'feature_q_value'.
#' @param add_global_nominal_threshold A logical value indicating whether to add global nominal thresholds. Default is FALSE.
#' @param global_nominal_threshold_column_to_add A character string specifying the name of the column to store global nominal thresholds. Default is 'pval_nominal_threshold_global'.
#' @param alpha_column A character string specifying the name of the column with alpha parameters for the beta distribution. Default is 'alpha_param'.
#' @param beta_column A character string specifying the name of the column with beta parameters for the beta distribution. Default is 'beta_param'.
#'
#' @return None. The function writes the filtered data to the specified output location.
#' @export
#'
#' @examples
#' \dontrun{
#'   filter_file_by_significance(
#'     input_loc = "path/to/input/file.txt",
#'     output_loc = "path/to/output/file.txt",
#'     significance_column = 'p_value',
#'     significance_cutoff = 0.05,
#'     verbose = TRUE,
#'     add_mtc = TRUE,
#'     mtc_column = 'empirical_feature_p_value',
#'     feature_mtc_column = 'feature_id',
#'     mtc_column_to_add = 'feature_q_value',
#'     add_global_nominal_threshold = FALSE,
#'     global_nominal_threshold_column_to_add = 'pval_nominal_threshold_global',
#'     alpha_column = 'alpha_param',
#'     beta_column = 'beta_param'
#'   )
#' }
filter_file_by_significance <- function(input_loc, output_loc, significance_column='p_value', significance_cutoff=0.05, verbose=T, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', add_global_nominal_threshold=F, global_nominal_threshold_column_to_add='pval_nominal_threshold_global', alpha_column='alpha_param', beta_column='beta_param') {
  # paste together the full path
  full_cell_type_path <- input_loc
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
    # based on this MTC column, we can now also add a cuttoff
    if (add_global_nominal_threshold) {
      cell_type_output_global_threshold <- calculate_nominal_thresholds(cell_type_output_features, fdr=significance_cutoff, pval_col=significance_column, nominal_threshold_column='nomthres', cutoff_column=mtc_column_to_add, alpha_column = alpha_column, beta_column = beta_column)
      print(head(cell_type_output_global_threshold))
      # now add the nominal threshold to the full table
      cell_type_output[[global_nominal_threshold_column_to_add]] <- cell_type_output_global_threshold[match(cell_type_output[[feature_mtc_column]], cell_type_output_global_threshold[[feature_mtc_column]]), 'nomthres'][['nomthres']]
    }
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
  # store the gz connection if we need it
  full_output_loc_wzip <- output_loc
  # gzip it if the extention ends on gz
  if (grepl('.gz$', output_loc)) {
    full_output_loc_wzip <- gzfile(output_loc)
  }
  # write result
  write.table(cell_type_output, full_output_loc_wzip, sep = '\t', row.names= F, col.names = T)
  # create md5
  mdfiver::create_md5_for_file(output_loc)
}


#' get the number eGenes per cell type from QTL output
#' 
#' @param input_dir base location of the folder containing files to filter
#' @param input_file which output file to read for the results
#' @param output_dir base location of the folder to put filtered output
#' @param output_file_prepend prepend of filtered file to create
#' @param output_file_append append of filtered file to create
#' @param significance_cutoff cutoff for which to set significance
#' @param split_column which column to split the output on
#' @param add_mtc add multiple testing before filtering down
#' @param mtc_column the column of values to apply multiple testing on
#' @param feature_mtc_column the column that has the feature group to perform the multiple testing on
#' @param mtc_column_to_add the name of the column that has the mtc-corrected values
#' @param filter_alpha remove entries that have an abhorrant alpha param
#' @param verbose print progress
#' @param folders vector of folders to consider. optional, if not supplied, all folders will be considered
#' @param sep value separator in table
#' @returns 0 if success
#' 
split_output_by_column <- function(input_dir, input_file='qtl_results_all.txt.gz', output_dir=NULL, output_file_prepend='qtl_results_all_qval_', output_file_append='.txt.gz', significance_cutoff=0.05, split_column='feature_chromosome', add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', add_global_nominal_threshold=T, global_nominal_threshold_column_to_add='pval_nominal_threshold_global', alpha_column='alpha_param', beta_column='beta_param', nominal_p_column='p_value', verbose=T, filter_alpha=T, folders=NULL, sep='\t') {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(input_dir, full.names = F, recursive = F)
  # subset if a set of folders was supplied
  if (!is.null(folders)) {
    cell_types <- intersect(cell_types, folders)
  }
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
    cell_type_output <- fread(full_cell_type_path, sep = sep, header = T)
    
    # filter on alpha if requested
    if (filter_alpha) {
      cell_type_output <- cell_type_output[!(cell_type_output[['alpha_param']] > 5 | cell_type_output[['alpha_param']] < .2), ]
    }
    
    # get the features and the emperical p value
    if (add_mtc) {
      if (verbose) {
        print(paste('adding MTC to', full_cell_type_path))
      }
      # if the feature mtc column was multiple, combine them
      if (length(feature_mtc_column) > 1) {
        # by taking the first column
        cell_type_output[[paste(feature_mtc_column, collapse = '_')]] <- cell_type_output[[feature_mtc_column[1]]]
        # and adding the other columns
        for (i in 1 : length(feature_mtc_column)) {
          cell_type_output[[paste(feature_mtc_column, collapse = '_')]] <- paste(cell_type_output[[paste(feature_mtc_column, collapse = '_')]], cell_type_output[[feature_mtc_column[i]]])
        }
        # and setting the new column as the feature column
        feature_mtc_column <- paste(feature_mtc_column, collapse = '_')
      }
      # subset to what we need
      cell_type_output_features <- NULL
      # which is a bit if we care about the nominal threshold
      if (add_global_nominal_threshold) {
        cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column, nominal_p_column, alpha_column, beta_column), with = F]
      }
      # even less if we don't try to get the nominal threshold as well
      else {
        cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column), with = F]
      }
      # order by significance
      cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
      # remove where the feature is smaller than zero
      cell_type_output_features <- cell_type_output_features[!(cell_type_output_features[[mtc_column]] < 0), ]
      # keep only the first entry
      cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
      # set the values that are larger than 1, to be 1, problem with precision
      cell_type_output_features[cell_type_output_features[[mtc_column]] > 1, mtc_column] <- 1
      # add multiple testing correction
      cell_type_output_features[['qvalue']] <- qvalue(cell_type_output_features[[mtc_column]])$qvalues
      # now add back to the original table
      cell_type_output[[mtc_column_to_add]] <- cell_type_output_features[match(cell_type_output[[feature_mtc_column]], cell_type_output_features[[feature_mtc_column]]), 'qvalue'][['qvalue']]
      # based on this MTC column, we can now also add a cuttoff
      if (add_global_nominal_threshold) {
        if (verbose) {
          print(paste('adding global nominal p-value cutoff to', full_cell_type_path))
        }
        # get the nominal p value cutoff based on the p values and the beta distribution
        cell_type_output_global_threshold <- calculate_nominal_thresholds(cell_type_output_features, fdr=significance_cutoff, pval_col=nominal_p_column, nominal_threshold_column='nomthres', cutoff_column='qvalue', alpha_column = alpha_column, beta_column = beta_column)
        # now add the nominal threshold to the full table
        cell_type_output[[global_nominal_threshold_column_to_add]] <- cell_type_output_global_threshold[match(cell_type_output[[feature_mtc_column]], cell_type_output_global_threshold[[feature_mtc_column]]), ][['nomthres']]
      }
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
#' @param folders optional vector of folders to look at specifically
#' @returns 0 if success
#' 
filter_output_by_significance <- function(unfiltered_loc, unfiltered_file='qtl_results_all.txt.gz', filtered_loc=NULL, filtered_file=NULL, significance_column='p_value', significance_cutoff=0.05, verbose=T, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', folders=NULL) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(unfiltered_loc, full.names = F, recursive = F)
  # check if overlaps with the folders that we want to take a look at
  if (!is.null(folders)) {
    cell_types <- intersect(cell_types, folders)
  }
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


#' Merge Chromosome Output
#'
#' This function merges chromosome-specific output files from different cell types into a single file.
#'
#' @param input_dir A character string specifying the directory containing the input files.
#' @param input_prepend A character string specifying the prefix of the input files. Default is 'qtl_results_all_qval_'.
#' @param input_append A character string specifying the suffix of the input files. Default is '_fdr01_significant.txt.gz'.
#' @param output_dir A character string specifying the directory to save the merged output file. Default is NULL.
#' @param output_file A character string specifying the name of the merged output file. Default is 'qtl_results_all_qval_allchroms_fdr01_significant.txt.gz'.
#' @param folders A character vector specifying the folders to include. Default is NULL.
#'
#' @return An integer value indicating the success of the operation. The function writes the merged data to the specified output location.
#' @export
#'
#' @examples
#' \dontrun{
#'   merge_chromosome_output(
#'     input_dir = "path/to/input/dir",
#'     input_prepend = 'qtl_results_all_qval_',
#'     input_append = '_fdr01_significant.txt.gz',
#'     output_dir = "path/to/output/dir",
#'     output_file = 'qtl_results_all_qval_allchroms_fdr01_significant.txt.gz',
#'     folders = c("folder1", "folder2")
#'   )
#' }
merge_chromosome_output <- function(input_dir, input_prepend='qtl_results_all_qval_', input_append='_fdr01_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr01_significant.txt.gz', folders=NULL) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(input_dir, full.names = F, recursive = F)
  # subset to folders we are interested in, if supplied with that option
  if (!is.null(folders)) {
    cell_types <- intersect(cell_types, folders)
  }
  # we will store the results in a list for now
  numbers_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the cell type folder
    input_celltype_dir <- paste(input_dir, '/', cell_type, '/', sep = '')
    # list the files
    input_files_celltype <- list.files(input_celltype_dir)
    # now filter only for the ones we want
    input_files_celltype <- input_files_celltype[grep(paste(input_prepend, '\\d+', input_append, '$', sep  = ''), input_files_celltype)]
    # we'll store each file in a list
    input_files_celltype_list <- list()
    # and go through each file
    for (input_file_celltype in input_files_celltype) {
      # read the file
      input_celltype_chrom <- fread(paste(input_celltype_dir, '/', input_file_celltype, sep = ''), header = T, sep = '\t')
      # put in the list
      if (nrow(input_celltype_chrom) > 0) {
        input_files_celltype_list[[input_file_celltype]] <- input_celltype_chrom
      }
      else {
        warning(paste('skipping', paste(input_celltype_dir, '/', input_file_celltype, sep = ''), 'because it has no rows'))
      }
    }
    # now merge all of them
    input_celltypes_all <- do.call('rbind', input_files_celltype_list)
    # get the output location
    output_location <- input_dir
    # if supplied, set the output directory
    if (!is.null(output_dir)) {
      output_location <- output_dir
    }
    # make full output location
    full_output_loc <- paste(output_location, '/', cell_type, '/', output_file, sep = '')
    # store the gz connection if we need it
    full_output_loc_wzip <- full_output_loc
    # gzip it if the extention ends on gz
    if (grepl('.gz$', full_output_loc)) {
      full_output_loc_wzip <- gzfile(full_output_loc)
    }
    # write result
    write.table(input_celltypes_all, full_output_loc_wzip, sep = '\t', row.names= F, col.names = T)
    # create md5
    mdfiver::create_md5_for_file(full_output_loc)
  }
  return(0)
}


perform_qvalue_correction <- function(cell_type_output, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value') {
  # subset to only the important columns
  if (is.data.table(cell_type_output)){
    cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column), with = F]
  }
  else {
    cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column)]
  }
  # order by significance
  cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
  # remove where the feature is smaller than zero
  cell_type_output_features <- cell_type_output_features[!(cell_type_output_features[[mtc_column]] < 0), ]
  # keep only the first entry
  cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
  # set the values that are larger than 1, to be 1, problem with precision
  cell_type_output_features[cell_type_output_features[[mtc_column]] > 1, mtc_column] <- 1
  # add multiple testing correction
  cell_type_output_features[['qvalue']] <- qvalue::qvalue(cell_type_output_features[[mtc_column]])$qvalues
  # now add back to the original table
  cell_type_output[[mtc_column_to_add]] <- cell_type_output_features[match(cell_type_output[[feature_mtc_column]], cell_type_output_features[[feature_mtc_column]]), ][['qvalue']]
  return(cell_type_output)
}


filter_interactions_by_qtls <- function(input_dir_interactions, input_dir_qtls, output_dir_interactions=NULL, output_file_interactions='/inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', input_interactions_filename='/inflammation_final/iqtl_results_all_eigenmt.tsv.gz', input_qtls_prepend='qtl_results_all_qval_', input_qtls_append='_fdr005_significant.txt.gz', significance_cutoffs_qlts=list('feature_q_value'=0.05), feature_column_qtls='feature_id', feature_column_interactions='feature', add_qvalue_interactions=T, qvalue_column_interactions='feature_q_value', column_to_mtc_interactions='feature_bf_eigen', add_total_eigen=T, total_eigen_column_interactions='total_bf_eigen', nominal_p_column_interactions='p_value', n_tests_column_interaction='n_tests_feature') {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(input_dir_interactions, full.names = F, recursive = F)
  # get the output location
  output_location <- input_dir_interactions
  # if supplied, set the output directory
  if (!is.null(output_dir_interactions)) {
    output_location <- output_dir_interactions
  }
  # we will also store the results
  results_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the path to the interactions
    input_file_interactions <- paste(input_dir_interactions, '/', cell_type, '/', input_interactions_filename, sep = '')
    # paste together the cell type folder for the qtls
    input_celltype_dir_qtls <- paste(input_dir_qtls, '/', cell_type, '/', sep = '')
    # check if the file and directory exist
    if (file.exists(input_file_interactions) & dir.exists(input_celltype_dir_qtls)) {
      # list the files
      input_files_celltype_qtls <- list.files(input_celltype_dir_qtls, full.names = F, recursive = F)
      # now filter only for the ones we want
      input_files_celltype_qtls <- input_files_celltype_qtls[grep(paste(input_qtls_prepend, '\\d+', input_qtls_append, '$', sep  = ''), input_files_celltype_qtls)]
      # we'll save the entire table in parts first
      qtl_files_list <- list()
      # and go through each file
      for (input_file_celltype in input_files_celltype_qtls) {
        # read the file
        input_single_file <- read.table(paste(input_celltype_dir_qtls, input_file_celltype, sep = '/'), header = T, sep = '\t')
        # subset using subsets
        if (!is.null(significance_cutoffs_qlts)) {
          # where each key is the column, and the cutoff the value
          for (cutoff_column in names(significance_cutoffs_qlts)) {
            input_single_file <- input_single_file[!is.na(input_single_file[[cutoff_column]]) & input_single_file[[cutoff_column]] < significance_cutoffs_qlts[[cutoff_column]], ]
          }
          # add to the list
          qtl_files_list[[input_file_celltype]] <- input_single_file
        }
      }
      # merge all
      qtl_files_all <- do.call('rbind', qtl_files_list)
      # get the features
      qtl_files_all_features <- qtl_files_all[[feature_column_qtls]]
      # read the interactions
      input_interactions <- read.table(input_file_interactions, header = T, sep = '\t')
      # subset to the features that were significant before
      input_interactions <- input_interactions[input_interactions[[feature_column_interactions]] %in% qtl_files_all_features, ]
      # (re-) add qvalue correction
      if (add_qvalue_interactions) {
        input_interactions <- perform_qvalue_correction(input_interactions, mtc_column = column_to_mtc_interactions, feature_mtc_column = feature_column_interactions, mtc_column_to_add = qvalue_column_interactions)
      }
      # add the total eigenMT
      if (add_total_eigen) {
        # get the sum of tests
        n_tests <- sum(input_interactions[!duplicated(input_interactions[[feature_column_interactions]]), n_tests_column_interaction])
        # bonferroni
        input_interactions[[total_eigen_column_interactions]] <- input_interactions[[nominal_p_column_interactions]] * n_tests
        # but of course no more than 1
        input_interactions[input_interactions[[total_eigen_column_interactions]] > 1, total_eigen_column_interactions] <- 1
      }
      # paste together the full output location
      full_output_loc <- paste(output_location, cell_type, output_file_interactions, sep = '/')
      # store the gz connection if we need it
      full_output_loc_wzip <- full_output_loc
      # gzip it if the extention ends on gz
      if (grepl('.gz$', full_output_loc)) {
        full_output_loc_wzip <- gzfile(full_output_loc)
      }
      # write the file
      write.table(input_interactions, full_output_loc_wzip, row.names = F, col.names = T, sep = '\t')
      # make checksum
      mdfiver::create_md5_for_file(full_output_loc)
      # put in the list as well
      results_per_celltype[[cell_type]] <- input_interactions
    }
    else {
      warning(paste('skipping', cell_type, 'because either the interaction file or qtl directory does not exist:', input_dir_interactions, input_dir_qtls))
    }
  }
  return(results_per_celltype)
}





####################
# Main Code        #
####################


###################
# mo eQTLs        #
###################

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# for UT and 24hCA as well
eqtl_output_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/UT/'
eqtl_output_24hca_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/24hCA/'

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
  verbose=T, 
  folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)
split_output_by_column(
  input_dir=eqtl_output_ut_loc,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)
split_output_by_column(
  input_dir=eqtl_output_24hca_loc,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)

# check each chromosome
for (chrom in 1:22) {
  # in location
  in_file <- paste('qtl_results_all_qval_', chrom, '.txt.gz', sep = '')
  # now filter on FDR as well
  fdr_file <- paste('qtl_results_all_qval_', chrom, '_fdr005_significant.txt.gz', sep = '')
  filter_output_by_significance(
    unfiltered_loc=eqtl_output_loc, 
    unfiltered_file=in_file, 
    filtered_loc=NULL, 
    filtered_file=fdr_file, 
    significance_column='feature_q_value', 
    significance_cutoff=0.05, 
    verbose=T, 
    add_mtc = F, 
    folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  )
  filter_output_by_significance(
    unfiltered_loc=eqtl_output_ut_loc, 
    unfiltered_file=in_file, 
    filtered_loc=NULL, 
    filtered_file=fdr_file, 
    significance_column='feature_q_value', 
    significance_cutoff=0.05, 
    verbose=T, 
    add_mtc = F,  
    folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  )
  filter_output_by_significance(
    unfiltered_loc=eqtl_output_24hca_loc, 
    unfiltered_file=in_file, 
    filtered_loc=NULL, 
    filtered_file=fdr_file, 
    significance_column='feature_q_value', 
    significance_cutoff=0.05, 
    verbose=T, 
    add_mtc = F,  
    folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  )
}

# merge chromosome outputs
merge_chromosome_output(eqtl_output_ut_loc, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')
merge_chromosome_output(eqtl_output_24hca_loc, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')
merge_chromosome_output(eqtl_output_loc, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')


###################
# oneK1K eQTLs    #
###################

# do onek1k as well
eqtl_output_onek1k_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/output/L1/'

# perform splitting
split_output_by_column(
  input_dir=eqtl_output_onek1k_loc,
  #input_file='qtl_results_all.txt.gz',
  input_file='qtl_results_all.txt',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  #folders = c('monocyte', 'NK', 'DC'), 
  folders = c('CD8T')
)


#################
# LCL caQTLs    #
#################

# and the LCL data
caqtl_lcl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/replication/'

# perform splitting
split_output_by_column(
  input_dir=caqtl_lcl_output_loc,
  input_file='qtl_results_all.txt',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='p_value',
  feature_mtc_column=c('snp_id', 'feature_id'),
  mtc_column_to_add='qtl_q_value',
  verbose=T, 
  filter_alpha = F
)


#################
# mo caQTLs     #
#################

# location of the caQTL
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'
# for UT and 24hCA as well
caqtl_output_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/combined_output_50kb/L1/UT/'
caqtl_output_24hca_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/combined_output_50kb/L1/24hCA/'

# for caQTL as well
split_output_by_column(
  input_dir=caqtl_output_loc,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  folders = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)
split_output_by_column(
  input_dir=caqtl_output_ut_loc,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  sep = ',', 
  folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)
split_output_by_column(
  input_dir=caqtl_output_24hca_loc,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  sep = ',', 
  folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)

# check each chromosome
for (chrom in 1:22) {
  # in location
  in_file <- paste('qtl_results_all_qval_', chrom, '.txt.gz', sep = '')
  # now filter on FDR as well
  fdr_file <- paste('qtl_results_all_qval_', chrom, '_fdr005_significant.txt.gz', sep = '')
  filter_output_by_significance(
    unfiltered_loc=caqtl_output_ut_loc,
    unfiltered_file=in_file,
    filtered_loc=NULL,
    filtered_file=fdr_file,
    significance_column='feature_q_value',
    significance_cutoff=0.05,
    verbose=T,
    add_mtc = F,
    folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  )
  filter_output_by_significance(
    unfiltered_loc=caqtl_output_24hca_loc,
    unfiltered_file=in_file,
    filtered_loc=NULL,
    filtered_file=fdr_file,
    significance_column='feature_q_value',
    significance_cutoff=0.05,
    verbose=T,
    add_mtc = F,
    folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  )
  filter_output_by_significance(
    unfiltered_loc=caqtl_output_loc, 
    unfiltered_file=in_file, 
    filtered_loc=NULL, 
    filtered_file=fdr_file, 
    significance_column='feature_q_value', 
    significance_cutoff=0.05, 
    verbose=T, 
    add_mtc = F, 
    folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  )
}

# merge chromosome outputs
merge_chromosome_output(caqtl_output_ut_loc, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')
merge_chromosome_output(caqtl_output_24hca_loc, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')
merge_chromosome_output(caqtl_output_loc, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')

# location of the interaction-eQTL outputs
icaqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/output/combined_significant/L1/'
# read the icaqtl output, and filter by ones that are FDR significant in the combined mapping
filter_interactions_by_qtls(icaqtl_output_loc, caqtl_output_loc)

########################
# sc-eQTLgen eQTLs     #
########################

# now do the eQTL output of sc-eQTLgen
sceqtlgen_base_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_20250212_freeze1/'
# for caQTL as well
split_output_by_column(
  input_dir=sceqtlgen_base_loc,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  folders = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)
# check each chromosome
for (chrom in 1:22) {
  # in location
  in_file <- paste('qtl_results_all_qval_', chrom, '.txt.gz', sep = '')
  # now filter on FDR as well
  fdr_file <- paste('qtl_results_all_qval_', chrom, '_fdr005_significant.txt.gz', sep = '')
  filter_output_by_significance(
    unfiltered_loc=sceqtlgen_base_loc,
    unfiltered_file=in_file,
    filtered_loc=NULL,
    filtered_file=fdr_file,
    significance_column='feature_q_value',
    significance_cutoff=0.05,
    verbose=T,
    add_mtc = F,
    folders=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  )
}
# merge significant results
merge_chromosome_output(sceqtlgen_base_loc, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')


###########################
# test directory QTLs     #
###########################

# any directory with QTL files
qtl_test_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined_20perm/'
# and which cell types are there
qtl_test_celltypes <- 'monocyte'
# perform splitting
split_output_by_column(
  input_dir=qtl_test_dir,
  input_file='qtl_results_all.txt.gz',
  output_dir=NULL,
  output_file_prepend='qtl_results_all_qval_',
  output_file_append='.txt.gz',
  split_column='feature_chromosome',
  add_mtc=T,
  mtc_column='empirical_feature_p_value',
  feature_mtc_column='feature_id',
  mtc_column_to_add='feature_q_value',
  verbose=T, 
  folders=qtl_test_celltypes
)
# check each chromosome
for (chrom in 1:22) {
  # in location
  in_file <- paste('qtl_results_all_qval_', chrom, '.txt.gz', sep = '')
  # now filter on FDR as well
  fdr_file <- paste('qtl_results_all_qval_', chrom, '_fdr005_significant.txt.gz', sep = '')
  # check if the file exists
  if (file.exists(in_file)) {
    filter_output_by_significance(
      unfiltered_loc=qtl_test_dir,
      unfiltered_file=in_file,
      filtered_loc=NULL,
      filtered_file=fdr_file,
      significance_column='feature_q_value',
      significance_cutoff=0.05,
      verbose=T,
      add_mtc = F,
      folders=qtl_test_celltypes
    )
  } else {
    # otherwise warn
    warning(paste('missing expected input file', in_file, '!'))
  }
  
}
# merge significant results
merge_chromosome_output(qtl_test_dir, input_append='_fdr005_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz', folders = qtl_test_celltypes)
