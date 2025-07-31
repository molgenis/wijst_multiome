#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_merge_hybrid_cres.R
# Function: merged chunked CRE outputs
# Example: Rscript mo_merge_hybrid_cres.R \
# --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/B/ \
# --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/B/qtl_results_all.txt.gz
#
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(optparse)
library(mdfiver)
library(qvalue)


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



#' Merge QTL result chunks from a directory
# processed.
#'
#' @param chunk_directory Character string. Path to the directory containing chunk subdirectories.
#' @param filename_output Character string. Name of the file to read from each chunk subdirectory. Default is `'qtl_results_all.txt.gz'`.
#' @param chunk_pattern Character string. Regular expression pattern to filter relevant chunk subdirectories. Default is `'chr\\d+\\-\\d+\\-\\d+'`.
#'
#' @return A data.table containing the merged contents of all valid chunk files, with an additional column `'chunk'` indicating the source subdirectory.
#'
#' @import data.table
#' @export
#'
#' @examples
#' \dontrun{
#' merged_data <- merge_chunks_in_directory("results/chunks")
#' head(merged_data)
#' }
merge_chunks_in_directory <- function(chunk_directory, filename_output='qtl_results_all.txt.gz', chunk_pattern='chr\\d+\\-\\d+\\-\\d+') {
  # list all of the chunks
  chunks <- list.dirs(chunk_directory, recursive = F, full.names = F)
  # keep only the ones with the pattern
  if (!is.null(chunk_pattern)) {
    chunks <- chunks[grep(chunk_pattern, chunks)]
  }
  # we'll save each chunk in a list
  chunks_list <- list()
  # warn if we have no chunks
  if (length(chunks) == 0) {
    warning(paste('no chunks at', chunk_directory))
  }
  # and read each chunk
  for (chunk in chunks) {
    # paste together the full path
    chunk_path <- paste(chunk_directory, chunk, filename_output, sep = '/')
    # check if the file exists
    if (file.exists(chunk_path)) {
      # check the size of the file
      if (file.size(chunk_path) > 0) {
        # read the file
        chunk_output <- fread(chunk_path, header = T, sep = '\t')
        # add chunk info
        chunk_output[['chunk']] <- chunk
        # and put in the list
        chunks_list[[chunk]] <- chunk_output
      } else {
        warning(paste('file exists at', paste(chunk_directory, chunk, filename_output, sep = '/'), 'but no file is of size 0', filename_output))
      }
    } else {
      warning(paste('folder exists at', paste(chunk_directory, chunk, filename_output, sep = '/'), 'but no file is there with name', filename_output))
    }
  }
  # merge all the chunks
  chunks_all <- do.call('rbind', chunks_list)
  return(chunks_all)
}


#' Read and process CRE QTL output for a single cell type
#'
#' This function reads QTL results from a specified folder, applies optional filtering and multiple testing correction (MTC),
#' and computes nominal p-value thresholds for significance. It is designed to work with output from `merge_chunks_in_directory`.
#'
#' @param cre_output_folder Character string. Path to the folder containing CRE QTL output chunks.
#' @param filename_output Character string. Name of the file to read from each chunk subdirectory. Default is `'qtl_results_all.txt.gz'`.
#' @param significance_column Character string. Column name used to assess significance. Default is `'empirical_feature_p_value'`.
#' @param significance_cutoff Numeric. Significance threshold for filtering and threshold calculation. Default is `0.05`.
#' @param add_mtc Logical. Whether to compute and add multiple testing correction (q-values). Default is `TRUE`.
#' @param mtc_column Character string. Column used for multiple testing correction. Default is `'empirical_feature_p_value'`.
#' @param feature_mtc_column Character string. Column identifying features for MTC. Default is `'feature_id'`.
#' @param mtc_column_to_add Character string. Name of the new column to store q-values. Default is `'feature_q_value'`.
#' @param add_global_nominal_threshold Logical. Whether to compute and add a global nominal p-value threshold. Default is `TRUE`.
#' @param add_local_nominal_threshold Logical. Whether to compute and add local nominal p-value thresholds. Default is `TRUE`.
#' @param global_nominal_threshold_column_to_add Character string. Name of the column to store the global nominal threshold. Default is `'pval_nominal_threshold_global'`.
#' @param local_nominal_threshold_column_to_add Character string. Name of the column to store local nominal thresholds. Default is `'pval_nominal_threshold_local'`.
#' @param alpha_column Character string. Column name for the alpha parameter. Default is `'alpha_param'`.
#' @param beta_column Character string. Column name for the beta parameter. Default is `'beta_param'`.
#' @param nominal_p_column Character string. Column name for nominal p-values. Default is `'p_value'`.
#' @param filter_alpha Logical. Whether to filter based on alpha parameter range. Default is `TRUE`.
#' @param alpha_min Numeric. Minimum alpha value for filtering. Default is `0.2`.
#' @param alpha_max Numeric. Maximum alpha value for filtering. Default is `5`.
#' @param filter_significance Logical. Whether to filter the final output based on significance. Default is `FALSE`.
#'
#' @return A `data.table` containing the processed QTL results for the cell type, with optional q-values and nominal thresholds.
#'
#' @import data.table
#' @importFrom qvalue qvalue
#' @export
#'
#' @examples
#' \dontrun{
#' results <- read_cre_output_per_celltype("results/CRE_output")
#' head(results)
#' }
read_cre_output_per_celltype <- function(cre_output_folder, filename_output='qtl_results_all.txt.gz', significance_column='empirical_feature_p_value', significance_cutoff=0.05, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', add_global_nominal_threshold=T, add_local_nominal_threshold=T, global_nominal_threshold_column_to_add='pval_nominal_threshold_global', local_nominal_threshold_column_to_add='pval_nominal_threshold_local', alpha_column='alpha_param', beta_column='beta_param', nominal_p_column='p_value', filter_alpha=T, alpha_min=.2, alpha_max=5, filter_significance=F, gene_frac_exp_loc = NULL, gene_frac_exp_cutoff = NULL, gene_frac_feature_column='feature', gene_frac_frac_column='frac_exp', gene_column='feature_id', gene_frac_cutoff=NULL) {
  # initalize the variable
  cell_type_output <- NULL
  # check if the file exists
  if (dir.exists(cre_output_folder)) {
    # read this file
    cell_type_output <- merge_chunks_in_directory(cre_output_folder)
    # make sure there are no duplicates
    cell_type_output <- unique(cell_type_output)
    # filter on alpha if requested
    if (filter_alpha) {
      cell_type_output <- cell_type_output[!(cell_type_output[[alpha_column]] > alpha_max | cell_type_output[[alpha_column]] < alpha_min), ]
    }
    # add the fraction of expression if we have that information
    if (!is.null(gene_frac_exp_loc)) {
      # read the file
      gene_fracs <- fread(gene_frac_exp_loc, header = T, sep = '\t')
      # add this data
      cell_type_output[['gene_frac_exp']] <- gene_fracs[match(cell_type_output[[gene_column]], gene_fracs[[gene_frac_feature_column]]), ][[gene_frac_frac_column]]
      # filter on this info if asked to
      if (!is.null(gene_frac_cutoff)) {
        cell_type_output <- cell_type_output[
          cell_type_output[['gene_frac_exp']] >= gene_frac_cutoff
        ]
      }
    }
    
    # get the features and the emperical p value
    if (add_mtc) {
      # subset to what we need
      cell_type_output_features <- NULL
      # which is a bit if we care about the nominal threshold
      if (add_global_nominal_threshold) {
        cell_type_output_features <- cell_type_output[, c(..feature_mtc_column, ..mtc_column, ..nominal_p_column, ..alpha_column, ..beta_column), with = F]
      }
      # even less if we don't try to get the nominal threshold as well
      else {
        cell_type_output_features <- cell_type_output[, c(..feature_mtc_column, ..mtc_column), with = F]
      }
      # remove the wherever we dont have our significance
      cell_type_output_features <- cell_type_output_features[!is.na(cell_type_output_features[[significance_column]]) & cell_type_output_features[[significance_column]] >= 0, ]
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
      if (add_local_nominal_threshold) {
        cell_type_output_local_threshold <- calculate_nominal_thresholds(cell_type_output_features, fdr=significance_cutoff, pval_col=nominal_p_column, nominal_threshold_column='nomthres', cutoff_column = 'qvalue', alpha_column = alpha_column, beta_column = beta_column)
        # now add the nominal threshold to the full table
        cell_type_output[[local_nominal_threshold_column_to_add]] <- cell_type_output_local_threshold[match(cell_type_output[[feature_mtc_column]], cell_type_output_local_threshold[[feature_mtc_column]]), 'nomthres'][['nomthres']]
      }
      if(add_global_nominal_threshold) {
        # filter the output to significant MTC hits
        cell_type_output_features_significant <- cell_type_output_features[cell_type_output_features[['qvalue']] < significance_cutoff, ]
        # and get the maximum significant nominal value
        global_p_cutoff <- max(cell_type_output_features_significant[[nominal_p_column]])
        # add that to the table
        cell_type_output[[global_nominal_threshold_column_to_add]] <- global_p_cutoff
      }
    }
    # filter the file if requested
    if (filter_significance) {
      cell_type_output <- cell_type_output[
        cell_type_output[[significance_column]] < significance_cutoff, 
      ]
    }
  }
  else {
    warning(paste('directory not found at', cre_output_folder))
  }
  return(cell_type_output)
}


####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- F


####################
# Main code        #
####################

# make command line options
option_list <- list(
  make_option(c("-i", "--in"), type="character", default=NULL, 
              help="input directory", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file with all data", metavar="character"), 
  make_option(c("-g", "--gene_frac_exp_loc"), type="character", default=NULL, 
              help="location of the file which has the fraction of cells expressing a gene", metavar="character"), 
  make_option(c("-f", "--gene_frac_exp_cutoff"), type="numeric", default=NULL, 
              help="minimum fraction of cells expressing a gene required to keep a gene in CRE results", metavar="numeric")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
qtl_in_loc <- NULL
qtl_out_loc <- NULL
gene_frac_exp_loc <- NULL
gene_frac_exp_cutoff <- NULL

# load debug settings if set to debug mode
if (debug) {
  qtl_in_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/NK/'
  qtl_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/NK/qtl_results_all_frac01.txt.gz'
  gene_frac_exp_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/frac_exp/NK.tsv.gz'
  gene_frac_exp_cutoff <- 0.1
} else {
  # there are some things we cannot allow
  if (is.null(opt[['in']])) {
    stop('input directory must be supplied')
  } else {
    qtl_in_loc <- opt[['in']]
  }
  if (is.null(opt[['out']])) {
    stop('output QTL file must be supplied')
  } else {
    qtl_out_loc <- opt[['out']]
  }
  # stop if we are overwriting our source file
  if (qtl_in_loc == qtl_out_loc) {
    stop('input and output are the same, do not overwrite your source file!')
  }
  gene_frac_exp_loc <- opt[['gene_frac_exp_loc']]
  gene_frac_exp_cutoff <- opt[['gene_frac_exp_cutoff']]
}
# do the actual things
qtl_merged_all <- read_cre_output_per_celltype(qtl_in_loc, gene_frac_exp_loc = gene_frac_exp_loc, gene_frac_exp_cutoff = gene_frac_exp_cutoff)
# make the output location filehandle
output_loc_fh <- qtl_out_loc
# gz filehandle, if the output location ends with .gz
if (grepl('.gz$', output_loc_fh)) {
  # gzip if ends with .gz
  output_loc_fh <- gzfile(output_loc_fh)
}
# final step
message(paste('writing output', qtl_out_loc))
# write the resulting file
write.table(qtl_merged_all, output_loc_fh, row.names = F, col.names = T, sep = '\t', quote = F)
# make a checksum as well
mdfiver::create_md5_for_file(qtl_out_loc)

# and let them know we are done
message('finished')
