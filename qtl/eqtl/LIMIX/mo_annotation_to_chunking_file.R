#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_annotation_to_chunking_file.R
# Function: convert an annotation file to a chunking file for LIMIX-QTL
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


#' Convert Annotation File to Chunking File
#'
#' This function reads an annotation file and converts it into a chunking file based on specified parameters.
#'
#' @param annotation_file_loc Character. The location of the annotation file.
#' @param chunking_file_loc Character. The location where the chunking file will be saved.
#' @param feature_chromosome_column Character. The column name for chromosome in the annotation file. Default is 'chromosome'.
#' @param feature_start_column Character. The column name for feature start positions in the annotation file. Default is 'start'.
#' @param feature_end_column Character. The column name for feature end positions in the annotation file. Default is 'end'.
#' @param strip_chr Character. A string to strip from the chromosome names. Default is 'chr'.
#' @param features_per_chunk Integer. The number of features per chunk. Default is 1.
#' @param include_chroms character vector. include only specific chromosomes. Default is NULL, meaning no filter.
#'
#' @return Integer. Returns 0 upon successful completion.
#' @export
#'
#' @examples
#' annotation_to_chunking("path/to/annotation_file.txt", "path/to/chunking_file.txt")
annotation_to_chunking <- function(annotation_file_loc, chunking_file_loc, feature_chromosome_column='chromosome', feature_start_column='start', feature_end_column='end', strip_chr='chr', features_per_chunk=1, include_chroms=NULL) {
  # read the annotation file
  annotation_file <- fread(annotation_file_loc, header = T, sep = '\t')
  # if we need to strip something from the chromosome name
  if (!is.null(strip_chr)) {
    annotation_file[[feature_chromosome_column]] <- gsub(strip_chr, '', annotation_file[[feature_chromosome_column]])
  }
  # create chunking file
  chunking_file <- NULL
  if (features_per_chunk == 1) {
    chunking_file <- data.frame(x = paste(annotation_file[[feature_chromosome_column]], ':', annotation_file[[feature_start_column]], '-', annotation_file[[feature_end_column]], sep = ''))
  }
  else {
    # create a data.table per chromosome
    chunking_table_per_chrom <- list()
    # check each chromosome
    for (chrom in unique(annotation_file[[feature_chromosome_column]])) {
      # subset to that chromosome
      results_chrom <- annotation_file[annotation_file[[feature_chromosome_column]] == chrom, ]
      # get the number of features for that chromosome
      n_features_chrom <- nrow(results_chrom)
      # get the starting chunks
      starts_i <- seq(from = 1, to = n_features_chrom, by = features_per_chunk)
      # and the ending chunks
      end_to <- n_features_chrom + features_per_chunk
      # get the stop positions
      ends_i <- seq(from = features_per_chunk, to = end_to, by = features_per_chunk)
      # now get the positions
      starts <- results_chrom[[feature_start_column]][starts_i]
      ends <- results_chrom[[feature_end_column]][ends_i]
      # if the start and end are not of the same length, that means the last chunk will be a bit smaller, it will be the end of the last feature
      if (ends_i[length(ends_i)] > n_features_chrom) {
        # so extract that one
        ends[length(ends_i)] <- results_chrom[[feature_end_column]][n_features_chrom]
      }
      # now put those back into the chunking file
      chunking_table_per_chrom[[chrom]] <- data.table('chrom' = rep(chrom, length(starts)), 'start' = starts, 'end' = ends)
    }
    # merge all
    chunkin_table_all <- do.call('rbind', chunking_table_per_chrom)
    # have only autosomal chromosomes for example if so requested
    if (!is.null(include_chroms)) {
      chunkin_table_all <- chunkin_table_all[chunkin_table_all[['chrom']] %in% include_chroms, ]
    }
    # and make the one we care about
    chunking_file <- data.frame(x = paste(chunkin_table_all[['chrom']], ':', chunkin_table_all[['start']], '-', chunkin_table_all[['end']], sep = ''))
  }
  # gz file ends with .gz
  output_file_chunking_gz <- chunking_file_loc
  if (grepl('.gz$', chunking_file_loc)) {
    output_file_chunking_gz <- gzfile(chunking_file_loc)
  }
  # write
  write.table(chunking_file, output_file_chunking_gz, row.names = F, col.names = F, quote = F)
  # and make an md5
  mdfiver::create_md5_for_file(chunking_file_loc)
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

# where the output should be placed
annotation_file_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/input/LimixAnnotationFile.txt'
chunking_file_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/chunking_files/ChunkingFile.txt.gz'

# convert the annotation file to a chunking file
annotation_to_chunking(
  annotation_file_loc, 
  chunking_file_loc, 
  features_per_chunk=100, 
  include_chroms = as.character(1:22)
)
