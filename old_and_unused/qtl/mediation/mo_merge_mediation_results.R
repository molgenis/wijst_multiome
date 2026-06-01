#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_merge_mediation_results.R
# Function: merge the results of each mediation job per chromosome into one file with MTC applied
############################################################################################################################


####################
# libraries        #
####################

library(data.table)


####################
# Functions        #
####################

#' Get Merged Table
#'
#' This function reads multiple files from a specified directory, filters them based on a given cell type and file pattern, merges the data, and performs significance filtering and multiple testing correction.
#'
#' @param input_dir Character. The directory containing the input files.
#' @param cell_type Character. The cell type to filter the files.
#' @param output_file Character. The path to the output file.
#' @param file_prepend Character. A string to prepend to the file names (default is '').
#' @param file_append Character. A string to append to the file names (default is '_chr\\d+.tsv.gz').
#' @param drop_insignificant_formula Logical. Whether to drop entries with insignificant formulas (default is TRUE).
#' @param formula_significance_cutoff Numeric. The significance cutoff for filtering formulas (default is 0.05).
#'
#' @return Integer. Returns 0 upon successful completion.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Lists the files in the specified directory.
#'   \item Filters the files based on the provided cell type and file pattern.
#'   \item Reads each filtered file and stores the data in a list.
#'   \item Merges all the data tables into one.
#'   \item Optionally filters out entries with insignificant direct and indirect effects.
#'   \item Applies multiple testing correction to the mediated p-values.
#'   \item Sorts the data by proportion mediated and p-value.
#'   \item Writes the merged data to the specified output file and creates an MD5 checksum.
#' }
#'
#' @examples
#' \dontrun{
#' get_merged_table(input_dir = "/path/to/data/", cell_type = "monocyte", output_file = "/path/to/data/monocyte_merged_data.tsv.gz")
#' }
#'
#' @import data.table
#' @import mdfiver
get_merged_table <- function(input_dir, cell_type, output_file, file_prepend='', file_append='_chr\\d+.tsv.gz', drop_insignificant_formula=T, formula_significance_cutoff=0.05) {
  # list the files in the directory
  files_in_dir <- list.files(input_dir, recursive = F, full.names = F)
  # paste together the regular expression pattern
  files_regex_pattern <- paste('^', file_prepend, cell_type, file_append, '$', sep = '')
  # filter the files in the directory
  files_in_dir <- files_in_dir[grep(files_regex_pattern, files_in_dir)]
  # create a list to store the table of each chromosome
  table_per_file <- list()
  # read each file
  for (file_in_dir in files_in_dir) {
    # read the file
    table_file <- fread(paste(input_dir, '/', file_in_dir, sep = ''), header = T, sep = '\t')
    # put into the list
    table_per_file[[file_in_dir]] <- table_file
  }
  # merge all of the separate tables
  all_files <- do.call('rbind', table_per_file)
  # metrics
  message(paste('tested', as.character(length(unique(all_files[['variant']]))), 'variants and', as.character(length(unique(all_files[['accessibility']]))), 'regions on', as.character(length(unique(all_files[['gene']]))), 'genes'))
  # remove the entries for which one of the formulas was already not significant
  if (drop_insignificant_formula) {
    # get the number of entries
    nr_possible_mediations <- nrow(all_files)
    # filter by direct
    all_files <- all_files[all_files[['p_direct']] < formula_significance_cutoff, ]
    nr_direct_filtered <- nrow(all_files)
    # let the user know
    message(paste('dropped', as.character(nr_possible_mediations - nr_direct_filtered), 'entries due to direct effect being insignificant'))
    # filter by indirected
    # filter by direct
    all_files <- all_files[all_files[['p_indirect']] < formula_significance_cutoff, ]
    nr_indirect_filtered <- nrow(all_files)
    # let the user know
    message(paste('dropped', as.character(nr_direct_filtered - nr_indirect_filtered), 'entries due to indirect effect being insignificant'))
    # # filter by total p
    # all_files <- all_files[all_files[['p_total']] < formula_significance_cutoff, ]
    # nr_ptotal_filtered <- nrow(all_files)
    # # let the user know
    # message(paste('dropped', as.character(nr_possible_mediations - nr_ptotal_filtered), 'entries due to total effect being insignificant'))
  }
  # do mtc on the mediated p-value
  all_files[['p_mediated_bonferroni']] <- p.adjust(all_files[['p_mediated']], method = 'bonferroni')
  all_files[['p_mediated_BH']] <- p.adjust(all_files[['p_mediated']], method = 'BH')
  # sort by proportion mediated first
  all_files <- all_files[order(all_files[['prop_med']], decreasing = T), ]
  # sort by p-value
  all_files <- all_files[order(all_files[['p_mediated']]), ]
  # write output
  output_loc_gz <- output_file
  # gz file ends with .gz
  if (grepl('.gz$', output_file)) {
    # gzip if ends with .gz
    output_loc_gz <- gzfile(output_file)
  }
  write.table(all_files, output_loc_gz, sep = '\t', row.names = F, col.names = T, quote = F)
  # also make a checksum
  mdfiver::create_md5_for_file(output_file)
  return(0)
}


do_debug <- function() {
  # debug list of options
  options <- list()
  options[['input_dir']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/output/regressed_gausnormed_atac_to_gausnormed_expression/UT/'
  options[['cell_type']] <- 'monocyte'
  options[['file_prepend']] <- ''
  options[['file_append']] <- '_chr\\d+.tsv.gz'
  options[['output_file']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/output/regressed_gausnormed_atac_to_gausnormed_expression/UT/monocyte_merged_minus1.tsv.gz'
  # run
  get_merged_table(
    options[['input_dir']], 
    options[['cell_type']], 
    options[['output_file']], 
    options[['file_prepend']], 
    options[['file_append']]
  )
}


####################
# Settings        #
####################


####################
# Debugging        #
####################

# debug function call, commented out when script is finished
do_debug()

####################
# Main Code        #
####################

