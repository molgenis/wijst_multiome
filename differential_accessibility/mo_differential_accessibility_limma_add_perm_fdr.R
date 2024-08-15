#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_differential_accessibility_limma_add_perm_fdr.R
# Function: add permutation-based FDR to limma results
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(optparse)


####################
# Functions        #
####################


create_permutation_p_values_distribution <- function(permutation_directory, permutation_prepend, permutation_append, pval_column='P.Value', feature_column='feature') {
  # list all the permutation files
  permutation_file_regex <- paste(permutation_prepend, '\\d+', permutation_append, '$', sep = '')
  # list all the files
  permutation_files <- list.files(permutation_directory, permutation_file_regex)
  # keep the sorted p values in a list
  p_values_per_permutation <- list()
  # let's read the files
  for (perm_file in permutation_files) {
    # read only the p value
    p_values_per_permutation[[perm_file]] <- fread(paste(permutation_directory, '/', perm_file, sep = ''), select = c(feature_column, pval_column))
    # set the P value column to be the name of the file
    colnames(p_values_per_permutation[[perm_file]]) <- c(feature_column, perm_file)
  }
  # now combine them
  p_values_permutations <- Reduce(function(...) merge(..., all = F, by = feature_column), p_values_per_permutation)
  # remove the feature column now
  p_values_permutations[[feature_column]] <- NULL
  # get the mean p per feature
  p_values_mins <- as.vector(apply(p_values_permutations, 1, FUN = min, na.rm = TRUE))
  return(p_values_mins)
}


calculate_permutation_fdr <- function(result_table, permutation_p_vector, result_table_p_column='P.Value', verbose=F) {
  # check the number of rows of the result table
  n_features <- nrow(result_table)
  # add a new column denoting the permutation-based FDR]
  result_table[['perm.FDR']] <- NA
  # check how many permuted P values we have
  n_permuted_ps <- length(permutation_p_vector)
  # sort both of these
  result_table <- result_table[order(result_table[[result_table_p_column]], decreasing = T), ]
  permutation_p_vector <- permutation_p_vector[order(permutation_p_vector, decreasing = T)]
  # set an index of where we are 
  i_permuted_p_vector <- 1
  # and extract that p
  permuted_p_at_index <- permutation_p_vector[i_permuted_p_vector]
  # check each p value in the result table
  for (i_p_true in 1:nrow(result_table)) {
    # extract by index
    p_true <- result_table[i_p_true, result_table_p_column]
    # otherwise keep searching
      # walk through the permuted p values until one is found that is equal or smaller than the true one
      while(p_true < permuted_p_at_index & i_permuted_p_vector <= length(permutation_p_vector)) {
        # update the index  
        i_permuted_p_vector <- i_permuted_p_vector + 1
        # and update the value
        permuted_p_at_index <- permutation_p_vector[i_permuted_p_vector]
        if (verbose & (i_permuted_p_vector %% 1000) == 0) {
          message(paste('processed', as.character(i_permuted_p_vector), 'permuted P-values'))
        }
      }
      # if we stopped because of no more p values left
      if (i_permuted_p_vector == length(permutation_p_vector) & p_true < permuted_p_at_index) {
        result_table[i_p_true, 'perm.FDR'] <- 0
      }
      # if this last one is better
      if (i_permuted_p_vector == length(permutation_p_vector) & p_true == permuted_p_at_index) {
        result_table[i_p_true, 'perm.FDR'] <- 1 / length(permutation_p_vector)
      }
      # otherwise calculate the fraction of permuted p values that was better than the true one
      else {
        result_table[i_p_true, 'perm.FDR'] <- 1 - ((i_permuted_p_vector - 1) / length(permutation_p_vector))
      }
    if (verbose & (i_p_true %% 1000) == 0) {
      message(paste('processed', as.character(i_p_true), 'features'))
    }
  }
  return(result_table)
}

####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="unpermuted results file location", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="location to write permutation-based FDR file to [default= %default]", metavar="character"),
  make_option(c("-d", "--perm_directory"), type="character", default=NULL,
              help="directory with permuted results [default= %default]", metavar="character"),
  make_option(c("-s", "--perm_prepend"), type="character", default=NULL,
              help="what each permuted result starts with [default= %default]", metavar="character"),
  make_option(c("-e", "--perm_append"), type="character", default=NULL,
              help="what each permuted result ends with [default= %default]", metavar="character"),
  make_option(c("-t", "--threads"), type="numeric", default=8,
              help="number of threads to use [default= %default]", metavar="numeric"),
  make_option(c("-v", "--verbose"), type="character", default='TRUE',
              help="print status messages [default= %default]", metavar="character")            
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# check the parameters
if (is.null(opt[['input']])) {
  stop('please supply the input parameter, the unpermuted file location. For example my_results.tsv.gz')
} else if (is.null(opt[['output']])) {
  stop('please supply the output parameter, where to store the result of the permutation-based FDR file. For example my_results_wfdr.tsv.gz')
} else if (is.null(opt[['perm_directory']])) {
  stop('please supply the perm_directory parameter, the directory with permuted results. For example ./')
} else if (is.null(opt[['perm_prepend']])) {
  stop('please supply the perm_prepend parameter, what each permuted result starts with. For example \'perm_res\'')
} else if (is.null(opt[['perm_append']])) {
  stop('please supply the perm_append parameter, what each permuted result end with. For example \'.tsv.gz\'')
} else if (!file.exists(opt[['input']])) {
  stop(paste('input file', opt[['input']], 'does not exist'))
} else if (!dir.exists(opt[['perm_directory']])) {
  stop(paste('permutation directory', opt[['perm_directory']], 'does not exist'))
} else {
    # set verbosity
    verbose = F
    # get verbosity
    if (opt[['verbose']] %in% c('True', 'true', 'TRUE', 't', 'T', '1')) {
        verbose <- T
    }else if (opt[['verbose']] %in% c('False', 'false', 'FALSE', 'f', 'F', '0')) {
        verbose <- F
    }else {
        stop(paste('invalid option for verbosity, valid options are \'TRUE\' or \'FALSE\''))
    }
    if (verbose) {
        message(paste('reading permutations with pattern ', opt[['perm_directory']], opt[['perm_prepend']], '\\d+', opt[['perm_append']], sep = ''))
    }
    # get the minimal P-values
    minimal_permuted_ps <- create_permutation_p_values_distribution(opt[['perm_directory']], opt[['perm_prepend']], opt[['perm_append']])
    if (verbose) {
        message(paste('reading true result at', opt[['input']]))
    }
    # read the results table
    result_table <- read.table(opt[['input']], header = T, sep = '\t')
    # add the permuted p value
    result_table <- calculate_permutation_fdr(result_table = result_table, permutation_p_vector = minimal_permuted_ps, verbose = verbose)
    if(verbose) {
        message(paste('writing result at', opt[['output']]))
    }
    # save result
    write.table(result_table, gzfile(opt[['output']]), row.names = F, col.names = T, sep = '\t')
}