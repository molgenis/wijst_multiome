#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_significant_variant_feature_pairs.R
# Function: check each cell type and condition, and create confinements of significant variant-feature pairs
############################################################################################################################

####################
# libraries        #
####################

library(data.table)


####################
# Functions        #
####################

get_significant_pairs_per_cell_type <- function(qtl_output_loc, qtl_prepend='qtl_results_all_qval_', qtl_append='_nominally_significant.txt.gz$', variant_column='snp_id', feature_column='feature_id', significance_list=list('feature_q_value' = 0.05, 'empirical_feature_p_value' = 0.05)) {
  # store per cell type
  results_per_cell_type <- list()
  # list all the folders
  cell_type_folders <- list.dirs(qtl_output_loc, recursive = F, full.names = F)
  # go through each folder
  for (cell_type in cell_type_folders) {
    # get the cell type folder
    cell_type_folder <- paste(qtl_output_loc, '/', cell_type, sep = '')
    # paste together the regex
    file_pattern <- paste(qtl_prepend, '.*', qtl_append, sep = '')
    # now list the files
    qtl_files <- list.files(cell_type_folder, pattern = file_pattern, recursive = F)
    # only continue if we have files
    if (length(qtl_files) > 0) {
      # get each confinement in a list first
      qlt_files_per_file <- list()
      # now check each file
      for (qtl_file in qtl_files) {
        # read the fle
        qtl_output_file <- fread(paste(cell_type_folder, '/', qtl_file, sep = ''), header = T, sep = '\t')
        # now subset based no our parameters
        for (significance_column in names(significance_list)) {
          # subset the data
          qtl_output_file <- qtl_output_file[qtl_output_file[[significance_column]] < significance_list[[significance_column]], ]
        }
        # subset to the variant-feature combinations
        var_feature_qtl_output <- qtl_output_file[, c(..variant_column, ..feature_column)]
        # put that into the list
        qlt_files_per_file[[qtl_file]] <- var_feature_qtl_output
      }
      # merge all of the files
      qtls_celltype <- do.call('rbind', qlt_files_per_file)
      # put in list for this cell type
      results_per_cell_type[[cell_type]] <- qtls_celltype
    }
  }
  return(results_per_cell_type)
}

get_significant_pairs_per_condition <- function(qtl_output_loc, qtl_prepend='qtl_results_all_qval_', qtl_append='_nominally_significant.txt.gz$', variant_column='snp_id', feature_column='feature_id', significance_list=list('feature_q_value' = 0.05, 'empirical_feature_p_value' = 0.05)) {
  # store per condition
  results_per_condition <- list()
  # list all the folders
  condition_folders <- list.dirs(qtl_output_loc, recursive = F, full.names = F)
  # go through each folder
  for (condition in condition_folders) {
    # get the cell type folder
    condition_folder <- paste(qtl_output_loc, '/', condition, sep = '')
    # get the results for this condition
    results_condition <- get_significant_pairs_per_cell_type(condition_folder, qtl_prepend = qtl_prepend, qtl_append = qtl_append, variant_column = variant_column, feature_column = feature_column, significance_list = significance_list)
    # put into list
    results_per_condition[[condition]] <- results_condition
  }
  return(results_per_condition)
}

get_significant_pairs_per_celltype_merged_conditions <- function(qtl_output_loc, qtl_prepend='qtl_results_all_qval_', qtl_append='_nominally_significant.txt.gz$', variant_column='snp_id', feature_column='feature_id', significance_list=list('feature_q_value' = 0.05, 'empirical_feature_p_value' = 0.05)) {
  # store per condition first
  results_per_condition <- get_significant_pairs_per_condition(qtl_output_loc, qtl_prepend = qtl_prepend, qtl_append = qtl_append, variant_column = variant_column, feature_column = feature_column, significance_list = significance_list)
  # we will save the outer join of the conditions
  results_per_celltype <- list()
  # we'll go through the conditions
  for (condition in names(results_per_condition)) {
    # get the results of that condition
    results_condition <- results_per_condition[[condition]]
    # check each cell type
    for (cell_type in names(results_condition)) {
      # check if this celltype was already in the new list
      if (cell_type %in% names(results_per_celltype)) {
        # get the current results
        cell_type_results_cur <- results_per_celltype[[cell_type]]
        # get the new ones
        cell_type_results_new <- results_condition[[cell_type]]
        # merge these
        cell_types_results_both <- rbind(cell_type_results_new, cell_type_results_cur)
        # but make them unique
        cell_types_results_both <- unique(cell_types_results_both)
        # then put that back into the list
        results_per_celltype[[cell_type]] <- cell_types_results_both
      }
      # otherwise just add it to the list
      else {
        results_per_celltype[[cell_type]] <- results_condition[[cell_type]]
      }
    }
  }
  return(results_per_celltype)
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



####################
# Main Code        #
####################

# location of eQTL result files
eqtl_results_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/'

# get the var-feature links for the eQTL cell types
eqtl_var_feature_celltypes <- get_significant_pairs_per_celltype_merged_conditions(eqtl_results_loc)

# get the location where to put the confinements
confinement_eqt_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/confinements/'
# write them
write_confinements(eqtl_var_feature_celltypes, confinement_eqt_loc)


# get the var-feature links for caQTL cell types
caqtl_results_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output_to_sigs/L1/'
caqtl_var_feature_celltypes <- get_significant_pairs_per_celltype_merged_conditions(caqtl_results_loc, qtl_prepend = '', qtl_append = '.tsv', significance_list=list())
# write the confinements
confinement_caqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/confinements/'
write_confinements(caqtl_var_feature_celltypes, confinement_caqtl_loc)
