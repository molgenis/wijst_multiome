#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_limix_interaction_qtl_confinements_nc2022.R
# Function: create the limix-QTL compatible variant-feature files for the NC2022 dataset
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(mdfiver)


####################
# Functions        #
####################


get_effects_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all.txt.gz', keep_columns=NULL, significance_column_and_cutoffs=list('feature_q_value' = 0.05), verbose=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  effects_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(qtl_output_loc, '/', cell_type, '/', output_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    # filter the results on significance
    if (!is.null(significance_column_and_cutoffs)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
      }
      # filter based on each column
      for (significance_column in names(significance_column_and_cutoffs)) {
        # get the cutoff
        significance_cutoff <- significance_column_and_cutoffs[[significance_column]]
        # and filter
        cell_type_output <- cell_type_output[
          !is.na(cell_type_output[[significance_column]]) &
            cell_type_output[[significance_column]] < significance_cutoff, 
        ]
      }
      
      if (verbose) {
        print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
      }
    }
    # get only the columns we care about
    if(!is.null(keep_columns)) {
      cell_type_output <- cell_type_output[, keep_columns, drop = F]
    }
    # add to the list
    effects_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(effects_per_celltype)
}


write_confinements <- function(significant_variant_gene_list, confinements_loc, confinement_file_prepend='', confinement_file_append='_confinement.tsv.gz') {
  # check all the cell types
  for (cell_type in names(significant_variant_gene_list)) {
    # paste together the file location
    file_loc_out <- paste(confinements_loc, '/', confinement_file_prepend, cell_type, confinement_file_append, sep = '')
    # check if we need to gzip
    file_loc_string <- file_loc_out
    # gz file ends with .gz
    if (grepl('.gz$', file_loc_string)) {
      file_loc_out <- gzfile(file_loc_string)
    }
    # set colnames
    colnames(significant_variant_gene_list[[cell_type]]) <- c('snp_id', 'feature_id')
    write.table(significant_variant_gene_list[[cell_type]], file_loc_out, row.names = F, col.names = T, sep = '\t', quote = F)
  }
}


####################
# Settings        #
####################

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# location of the current interaction-eQTL output
interaction_eqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/L1/ut_and_24hca_significant/'
# get all significant effects
interaction_eqtl_effects <- get_effects_per_celltype_limix(interaction_eqtl_loc, output_file = 'iqtl_results_all_eigenmt_qval.tsv.gz', significance_column_and_cutoffs = list('feature_bf_eigen' = 0.05, 'feature_q_value' = 0.05), keep_columns = c('snp_id', 'feature'))

# location to put the confinements
confinement_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/confinements/replication/'
# write confinements
write_confinements(interaction_eqtl_effects, confinement_ieqtl_loc, confinement_file_append = '_mo_ut_and24hca_ieqtl_significant.tsv.gz')
