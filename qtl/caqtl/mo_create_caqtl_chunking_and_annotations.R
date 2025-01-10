#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_caqtl_chunking_and_annotations.R
# Function: create feature chunking and annotation files for caQTL mapping, based on percentage of cells having a peak as described in beds
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

#' Convert Peak Files to annotation and chunking files
#'
#' This function reads peak files for specified cell types and conditions, merges the conditions for a cell type, filters the data based on a given column and value, and writes the appropriate chunking and annotation files.
#'
#' @param peak_input_dir Character. Directory containing the peak files in bed format.
#' @param output_file_chunking Character. Path to file to save chunking to.
#' @param output_file_annotatations Character. Path to file to save annotations to.
#' @param cell_types Character vector. Cell types to process. Default is c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK').
#' @param conditions Character vector. Conditions to process. Default is c('UT', '24hCA').
#' @param file_prepend Character. Prefix for the input peak file names. Default is 'mo_peaks_lane1to80_'.
#' @param file_append Character. Suffix for the input peak file names. Default is '.bed'.
#' @param filter_column Character. Column name to apply the filter on. Default is 'pct_exp'.
#' @param filter_value Numeric. Minimum value for filtering the data. Default is 0.001.
#' @param feature_column Character. Column name to extract unique features. Default is 'name'.
#' @param strip_chr Character. What to strip from the beginning of the chromosome name column. Default is 'chr'.
#' @return Integer. Returns 0 upon successful completion.
#' @examples
#' \dontrun{
#' peak_files_to_annotation_and_chunks("path/to/peak_input", "path/to/chunkfile.txt.gz", "path/to/annotations.tsv.gz")
#' }
peak_files_to_annotation_and_chunks <- function(peak_input_dir, output_file_chunking, output_file_annotatations, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', '24hCA'), file_prepend='mo_peaks_lane1to80_', file_append='.bed', filter_column='pct_exp', filter_value=0.001, feature_column='name', strip_chr='chr') {
  # we'll store all results
  results_list <- list()
  # check cell type
  for (cell_type in cell_types) {
    # now check each condition
    for (condition in conditions) {
      # paste together the full path
      full_path_condition_celltype <- paste(peak_input_dir, '/', file_prepend, condition, '_', cell_type, file_append, sep = '')
      # check if this file exists
      if (file.exists(full_path_condition_celltype)) {
        # read the file
        condition_celltype <- fread(full_path_condition_celltype, header = T, sep = '\t', check.names = F)
        # put into the list
        results_list[[paste(condition, cell_type)]] <- condition_celltype
      }
      else {
        warning(paste('cell type', cell_type, 'and condition', condition, 'does not have a file at', full_path_condition_celltype, 'and is skipped'))
      }
    }
  }
  # merge conditions
  results_all <- do.call('rbind', results_list)
  # filter on minimum
  results_all <- results_all[results_all[[filter_column]] >= filter_value, ]
  # get the columns we care about
  results_all <- unique(results_all[, c('#chrom','start', 'end', 'name')])
  # just to be save, order them
  results_all <- results_all[order(results_all[['#chrom']], results_all[['start']], results_all[['end']]), ]
  # remove the chr if requested
  if (!is.na(strip_chr)) {
    results_all[['#chrom']] <- gsub(paste('^', strip_chr, sep = ''), '', results_all[['#chrom']])
  }
  
  # create chunking file
  chunking_file <- data.frame(x = paste(results_all[['#chrom']], ':', results_all[['start']], '-', results_all[['end']], sep = ''))
  # gz file ends with .gz
  output_file_chunking_gz <- output_file_chunking
  if (grepl('.gz$', output_file_chunking)) {
    output_file_chunking_gz <- gzfile(output_file_chunking)
  }
  # write
  write.table(chunking_file, output_file_chunking_gz, row.names = F, col.names = F, quote = F)
  # and make an md5
  mdfiver::create_md5_for_file(output_file_chunking)
  
  # create annotation file
  annotation_file <- results_all[, c('name', '#chrom', 'start', 'end')]
  # set colnames as expected
  colnames(annotation_file) <- c('feature_id', 'chromosome', 'start', 'end')
  # gz file ends with .gz
  output_file_annotatations_gz <- output_file_annotatations
  if (grepl('.gz$', output_file_annotatations)) {
    output_file_annotatations_gz <- gzfile(output_file_annotatations)
  }
  # write
  write.table(annotation_file, output_file_annotatations_gz, row.names = F, col.names = T, sep = '\t', quote = F)
  # and make an md5
  mdfiver::create_md5_for_file(output_file_annotatations)
  
  return(0)
}

####################
# Settings        #
####################

# set so that positions are not converted to scientific notations
options(scipen=999)


####################
# Main Code        #
####################

# where the peaks are
peak_pct_dir <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/signac_peaks/output/'
# where the output should be placed
annotation_file_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/annotations/pct01/LimixAnnotationFile.tsv.gz'
chunking_file_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/chunking_files/pct01/ChunkingFile.txt.gz'

peak_files_to_annotation_and_chunks(
  peak_pct_dir, 
  chunking_file_loc, 
  annotation_file_loc
)
