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


####################
# Functions        #
####################

#' get the number eGenes per cell type from QTL output
#' 
#' @param unfiltered_loc base location of the folder containing files to filter
#' @param unfiltered_file which output file to read for the results
#' @param filtered_loc base location of the folder to put filtered output
#' @param filtered_file name of filtered file to create
#' @param significance_column column denoting significance
#' @param significance_cutoff cutoff for which to set significance
#' @param verbose print progress
#' @returns 0 if success
#' 
filter_output_by_significance <- function(unfiltered_loc, unfiltered_file='qtl_results_all.txt.gz', filtered_loc=NULL, filtered_file=NULL, significance_column='p_value', significance_cutoff=0.05, verbose=T) {
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
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    
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
