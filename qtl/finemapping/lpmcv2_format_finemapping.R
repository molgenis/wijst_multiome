#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lpmcv2_format_finemapping.R
# Function: format the binary finemapping output into tsv format
# Example: 
# Rscript ~/lpmcv2_format_finemapping.R \
#   --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_eqtl/sc-eqtlgen/output/nominal_condition/L1//B_finemapped.rds \
#   --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_eqtl/sc-eqtlgen/output/nominal_condition/L1//B_finemapped.tsv.gz
############################################################################################################################


####################
# libraries        #
####################

library(coloc)
library(mdfiver)
library(optparse)


####################
# Functions        #
####################

#' Convert SuSiE RSS Output to Table
#'
#' This function takes unformatted SuSiE RSS output and converts it into a formatted table.
#'
#' @param unformatted_susie_output A list of unformatted SuSiE RSS output, where each element contains a feature and its corresponding result.
#' @param feature_key A list character string denoting the key used for the feature
#' @return A data frame containing the formatted results for each feature, including feature name, variant, lambda, PIP, residual variance estimation, credible sets, and log Bayes factors.
#' @examples
#' # Example usage:
#' unformatted_output <- list(
#'   list(feature = "feature1", result = list(susie_rss = list(converged = TRUE, alpha = matrix(0.1, ncol = 3), pip = c(0.2, 0.3, 0.5), sets = list(cs = list(1:2)), lbf_variable = matrix(0.4, ncol = 3)), lambda = 0.5)),
#'   list(feature = "feature2", result = list(susie_rss = list(converged = FALSE), lambda = 0.7))
#' )
#' formatted_table <- susie_rss_to_table(unformatted_output)
#' print(formatted_table)
susie_rss_to_table <- function(unformatted_susie_output, feature_key='feature') {
  # create a list of each feature
  formatted_fm_per_feature <- list()
  # and keep track of the max number of sets
  max_sets <- 0
  # check each of the outputs
  for (i in 1 : length(unformatted_susie_output)) {
    # extract the output
    unformatted_susie_output_feature <- unformatted_susie_output[[i]]
    # get the feature
    feature_output <- unformatted_susie_output_feature[[feature_key]]
    # and the result
    feature_result <- unformatted_susie_output_feature[['result']]
    # extract the susie_rss result
    feature_result_susierss <- feature_result[['susie_rss']]
    # check if converged
    if (feature_result_susierss$converged) {
      # create a table with the ensemble ID and the variants
      feature_result_res_table <- data.frame('feature' = rep(feature_output, times = ncol(feature_result_susierss$alpha)), 'variant' = colnames(feature_result_susierss$alpha))
      # add the lambda
      feature_result_res_table[['lambda']] <- feature_result[['lambda']]
      # add the pip
      feature_result_res_table[['pip']] <- feature_result_susierss[['pip']][match(feature_result_res_table[['variant']], names(feature_result_susierss[['pip']]))]
      # add whether we estimated the residual variance
      if ('resvar' %in% names(feature_result)) {
        feature_result_res_table[['resvar']] <- feature_result_susierss[['resvar']]
      }
      else {
        feature_result_res_table[['resvar']] <- NA
      }
      # and if it converged
      feature_result_res_table[['converged']] <- T
      # get the number of iterations we used
      feature_result_res_table[['n_iter']] <- feature_result_susierss[['niter']]
      # initialize the credible sets
      feature_result_res_table[['CS']] <- NA
      # check if we have credible sets
      if(length(feature_result_susierss[['sets']][[1]])!=0){
        # check each set
        for(l in names(feature_result_susierss[['sets']][['cs']])){
          # set for those positions, which set it was
          indices_variant_set <- feature_result_susierss[['sets']][['cs']][[l]]
          # and set for those
          feature_result_res_table[indices_variant_set, 'CS'] <- l
        }
      }
      # get the lbf, and transpose
      lbf_out <- data.frame(t(feature_result_susierss[['lbf_variable']]))
      # replace the X in the name, with CS
      colnames(lbf_out) <- paste0('CS', 1:ncol(lbf_out))
      # add in the correct order
      if (nrow(lbf_out) > 1) {
        feature_result_res_table <- cbind(feature_result_res_table, lbf_out[match(feature_result_res_table[['variant']], rownames(lbf_out)), , drop = F])
      }
      # except when there is only one variant, then this would mess up the column names
      else {
        feature_result_res_table <- cbind(feature_result_res_table, lbf_out)
      }
      # and put into the list
      formatted_fm_per_feature[[feature_output]] <- feature_result_res_table
      # update the max number of sets if this was higher
      if (ncol(lbf_out) > max_sets) {
        max_sets <- ncol(lbf_out)
      }
    }
    else {
      feature_result_res_table <- data.frame(
        'feature' = c(feature_output), 
        'variant' = c(NA), 
        'lambda' = c(NA), 
        'pip' = c(NA), 
        'resvar' = c(NA), 
        'converged' = c(F), 
        'n_iter' = c(feature_result_susierss[['niter']]), 
        'CS' = c(NA)
      )
      # and put into the list
      formatted_fm_per_feature[[feature_output]] <- feature_result_res_table
    }
  }
  # these are the columns that are always present
  permanent_columns <- c('feature', 'variant', 'lambda', 'pip', 'resvar', 'converged', 'n_iter', 'CS')
  # check each per-feature output
  for(feature_output in names(formatted_fm_per_feature)) {
    # and check if the number of columns is less than we would expect based on the max number of sets
    if (ncol(formatted_fm_per_feature[[feature_output]]) < (length(permanent_columns) + max_sets)) {
      # then get get which columns are missing
      missing_columns <- setdiff(paste0('CS', 1 : max_sets), colnames(formatted_fm_per_feature[[feature_output]]))
      # we need to make the number of column the same across the tables, or the rbind will fail. S we create an empty data frame
      missing_data <- data.frame(matrix(NA, nrow = nrow(formatted_fm_per_feature[[feature_output]]), ncol = length(missing_columns)))
      # with the correct row names
      rownames(missing_data) <- rownames(formatted_fm_per_feature[[feature_output]])
      # and the columns that we are missing from the total
      colnames(missing_data) <- missing_columns
      # and merge these missing columns onto the existing data
      formatted_fm_per_feature[[feature_output]] <- cbind(formatted_fm_per_feature[[feature_output]], missing_data)
    }
  }
  # merge all of them together
  formatted_fm_all <- do.call('rbind', formatted_fm_per_feature)
  return(formatted_fm_all)
}


####################
# Settings        #
####################

set.seed(7777)


####################
# Debug            #
####################



####################
# Main Code        #
####################


# make command line options
option_list <- list(
  make_option(c("-i", "--input_file"), type="character", default=NULL, 
              help="path to the binary finemapping output", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="path to the tab formatted file to write to", metavar="character"), 
  make_option(c("-k", "--feature_key"), type="character", default='feature', 
              help="key which has the feature name for each finemapped locus", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize the values we need
input_file <- NULL
output_file <- NULL
# check variables, and stop if these were not supplied
if (is.null(opt[['input_file']])) {
  stop(paste('-i/--input_file is an obligatory parameter\n'))
} else {
  input_file <- opt[['input_file']]
}
if (is.null(opt[['output_file']])) {
  stop(paste('-o/--output_file is an obligatory parameter\n'))
} else {
  output_file <- opt[['output_file']]
}
# check if the directory for the output exists, otherwise we would fail at the very last step
output_dir <- dirname(output_file)
if (!(dir.exists(output_dir))) {
  stop(paste('directory for output table', output_dir, 'does not exist. This would mean the final step writing to file would fail, please create the directory first\n'))
}
# check if we are not accidentally overwriting the source file
if (input_file == output_file) {
  stop(paste('output file', output_file, 'is the same as the input file. You should not overwrite your input! Please select a different output file.'))
}
# check the optional parameter
feature_key <- opt[['feature_key']]

# read the binary output
input_rds <- readRDS(input_file)
# convert to the table format
output_table <- susie_rss_to_table(input_rds, feature_key = feature_key)
# make the output table
output_path_full <- output_file
# gzip if that ends with
if (grepl('.gz$', output_file)) {
  output_path_full <- gzfile(output_file)
}
# write the result
write.table(output_table, output_path_full, row.names = F, col.names = T, sep = '\t', quote = F)
# and make a checksum
mdfiver::create_md5_for_file(output_file)