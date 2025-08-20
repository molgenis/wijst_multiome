#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_qtl_variant_to_region.R
# Function: get list of variants ever tested for QTLs
############################################################################################################################

####################
# libraries        #
####################

library(mdfiver)
library(data.table)


####################
# Functions        #
####################

read_qtl_files_all <- function(unfiltered_loc, folders=NULL, unfiltered_file='qtl_results_all.txt.gz', variant_id_column='snp_id', variant_chromosome_column='snp_chromosome', variant_position_column='snp_position', verbose = T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(unfiltered_loc, full.names = F, recursive = F)
  # check if overlaps with the folders that we want to take a look at
  if (!is.null(folders)) {
    cell_types <- intersect(cell_types, folders)
  }
  # initialize our final table
  all_variants <- NULL
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(unfiltered_loc, '/', cell_type, '/', unfiltered_file, sep = '')
    # check if the file exists
    if (file.exists(full_cell_type_path)) {
      # log if requested
      if (verbose) {
        print(paste('reading', full_cell_type_path))
      }
      # read the file
      cell_type_output <- fread(full_cell_type_path, sep = '\t', header = T)
      # subset to the variants
      cell_type_output[, c(..variant_id_column, ..variant_chromosome_column, ..variant_position_column)]
      # get only unique entries
      cell_type_output <- unique(cell_type_output)
      # and add to existing table
      if (!is.null(all_variants)) {
        all_variants <- rbind(all_variants, cell_type_output)
        # don't forget to make unique again
        all_variants <- unique(all_variants)
      }
      # if the table doesn't exist, we need to initialize it
      else {
        all_variants <- cell_type_output
      }
    }
    # warn
    else {
      warning(paste('directory for file exists, but file is not there, so skipped:', full_cell_type_path))
    }
  }
  # order the table
  all_variants <- all_variants[order(all_variants[[variant_chromosome_column]], all_variants[[variant_position_column]], all_variants[[variant_id_column]]), ]
  # and return the result
  return(all_variants)
}

####################
# Main Code        #
####################


###################
# mo eQTLs        #
###################

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# location of the caQTL
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'

# get all eQTL variants
eqtl_variants_all <- read_qtl_files_all(eqtl_output_loc)
caqtl_output_all <- read_qtl_files_all(caqtl_output_loc)
