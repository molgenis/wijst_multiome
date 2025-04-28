#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_outer_join_filtered_chromatin_regions.R
# Function: check the overlap of effects for colocalizing eQTLs with caQTLs
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(mdfiver)


####################
# Functions         #
####################

#' Get Outer Join Region Includes
#'
#' This function processes BED files for different cell types and conditions, filters the data based on specified criteria, and extracts unique genomic regions.
#'
#' @param beds_path Character. The directory path where the BED files are located.
#' @param prepend Character. A string to prepend to the file names. Default is 'mo_peaks_lane1to80_'.
#' @param midpend Character. A string to insert between the condition and cell type in the file names. Default is '_'.
#' @param append Character. A string to append to the file names. Default is '.bed'.
#' @param cell_types Character vector. A vector of cell types to process. Default is c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK').
#' @param conditions Character vector. A vector of conditions to process. Default is c('UT', '24hCA').
#' @param filters List. A named list of filters to apply to the data. Default is list('pct_exp' = 0.001).
#' @param chrom_sep Character. A string to separate chromosome names in the output. Default is ':'.
#' @param region_sep Character. A string to separate start and end positions in the output. Default is '-'.
#'
#' @return Character vector. A sorted vector of unique genomic regions in the format 'chromosome;start-end'.
#'
#' @examples
#' \dontrun{
#' regions <- get_outer_join_region_includes(
#'   beds_path = "/path/to/bed/files",
#'   prepend = "mo_peaks_lane1to80_",
#'   midpend = "_",
#'   append = ".bed",
#'   cell_types = c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"),
#'   conditions = c("UT", "24hCA"),
#'   filters = list("pct_exp" = 0.001),
#'   chrom_sep = ":",
#'   region_sep = "-"
#' )
#' }
get_outer_join_region_includes <- function(beds_path, prepend='mo_peaks_lane1to80_', midpend='_', append='.bed', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', '24hCA'), filters=list('pct_exp' = 0.001), chrom_sep=':', region_sep='-') {
  # we will store the regions in a list first
  regions_per_celltype_and_condition <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # next check each condition
    for (condition in conditions) {
      # paste together the full path
      bed_file_path <- paste0(beds_path, '/', prepend, condition, midpend, cell_type, append)
      # check if this file exists
      if (file.exists(bed_file_path)) {
        # read the file
        bed_file_contents <- data.table::fread(bed_file_path, header = T, sep = '\t')
        # now filter on the filters that we have
        for(filter_column in names(filters)) {
          # filter that column on that value
          bed_file_contents <- bed_file_contents[bed_file_contents[[filter_column]] >= filters[[filter_column]], ]
        }
        # extract the regions we care about
        regions_bed_file <- paste0(bed_file_contents[['#chrom']], chrom_sep, bed_file_contents[['start']], region_sep, bed_file_contents[['end']])
        # put into the list
        regions_per_celltype_and_condition[[bed_file_path]] <- regions_bed_file
      }
      # warn if not there
      else {
        warning(paste0('checking ', bed_file_path, ', but does not exist. Will be skipped.\n'))
      }
    }
  }
  # merge all of the regions from the different files
  regions_all <- do.call('c', regions_per_celltype_and_condition)
  # get only the unique ones
  regions_all <- unique(regions_all)
  # and sort them for good measure
  regions_all <- regions_all[order(regions_all)]
  return(regions_all)
}


####################
# Settings         #
####################


####################
# Main Code        #
####################

# location of the bed files
bed_files_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/signac_peaks/output/'
# get all of the regions
outer_join_regions <- get_outer_join_region_includes(bed_files_loc)
# where to place the resulting file
all_regions_list_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/merged_major_celltypes/mo_pct0001_regions.txt.gz'
# write to a file
write.table(data.frame(x = outer_join_regions), gzfile(all_regions_list_loc), row.names = F, col.names = F, quote = F)
# and make md5
mdfiver::create_md5_for_file(all_regions_list_loc)

# also when including the minor cell types
outer_join_regions_also_minor <- get_outer_join_region_includes(bed_files_loc, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK', 'plasmablast', 'T_other'))
all_regions_also_minor_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/merged_major_and_minor_celltypes/mo_pct0001_regions.txt.gz'
write.table(data.frame(x = outer_join_regions_also_minor), gzfile(all_regions_also_minor_loc), row.names = F, col.names = F, quote = F)
mdfiver::create_md5_for_file(all_regions_also_minor_loc)
