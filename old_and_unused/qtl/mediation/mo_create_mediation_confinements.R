#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_mediation_confinements.R
# Function: create confinements of variant-region-feature for mediation analysis, based on finemapped caQTLs and eQTLs
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

#' Create Mediation Confinements
#'
#' This function processes fine-mapped caQTLs and eQTLs to create mediation confinements for specified cell types.
#'
#' @param caqtls_finemapped Data frame or data.table containing fine-mapped caQTLs.
#' @param eqtls_dir Directory containing eQTL results.
#' @param output_dir Directory to save the output files.
#' @param cell_types_to_do Vector of cell types to process. If NULL, all cell types in caQTLs will be processed. Default is NULL.
#' @param output_file_prepend String to prepend to the output file names. Default is an empty string.
#' @param output_file_append String to append to the output file names. Default is '.confinement.tsv.gz'.
#' @param finemapped_celltype_column Column name for cell types in the fine-mapped caQTLs. Default is 'ct'.
#' @param finemapped_variant_column Column name for variants in the fine-mapped caQTLs. Default is 'snp_id'.
#' @param finemapped_feature_column Column name for features in the fine-mapped caQTLs. Default is 'feature_id'.
#' @param qtl_results_prepend String to prepend to the eQTL result file names. Default is 'qtl_results_all_qval_'.
#' @param qtl_results_append String to append to the eQTL result file names. Default is '_fdr01_significant.txt.gz'.
#' @param qtl_results_variant_column Column name for variants in the eQTL results. Default is 'snp_id'.
#' @param qtl_results_feature_column Column name for features in the eQTL results. Default is 'feature_id'.
#' @param qtl_resuls_filter_significances List of significance thresholds for filtering eQTL results, keys are columns, values are the thresholds. Default is list('empirical_feature_p_value' = 0.05, 'p_value' = 0.05).
#' @return Returns 0 upon successful completion.
#'
#' @examples
#' create_mediation_confinements(caqtls_finemapped, eqtls_dir, output_dir)
create_mediation_confinements <- function(caqtls_finemapped, eqtls_dir, output_dir, cell_types_to_do=NULL, output_file_prepend='', output_file_append='.confinement.tsv.gz', finemapped_celltype_column='ct', finemapped_variant_column='snp_id', finemapped_feature_column='feature_id', qtl_results_prepend='qtl_results_all_qval_', qtl_results_append='_fdr01_significant.txt.gz', qtl_results_variant_column='snp_id', qtl_results_feature_column='feature_id', qtl_resuls_filter_significances=list('empirical_feature_p_value' = 0.05, 'p_value' = 0.05)) {
  # get the cell types in the finemapped data
  cell_types_finemapped <- unique(caqtls_finemapped[[finemapped_celltype_column]])
  # only some cell types if asked
  if (!is.null(cell_types_to_do)) {
    cell_types_finemapped <- intersect(cell_types_finemapped, cell_types_to_do)
  }
  # check each of these
  for (cell_type in cell_types_finemapped) {
    # store the resuls of each file
    results_each_file_list <- list()
    # get the folder for that cell type in the eqtls
    eqtl_celltype_dir <- paste(eqtls_dir, '/', cell_type, '/', sep = '')
    # only do this cell type if the directory exists
    if (dir.exists(eqtl_celltype_dir)) {
      # subset the caqtls for this cell type
      caqtls_finemapped_celltype <- caqtls_finemapped[caqtls_finemapped[[finemapped_celltype_column]] == cell_type, ]
      # get the eqtl file pattern
      eqtl_file_pattern <- paste('^', qtl_results_prepend, '.*', qtl_results_append, '$', sep = '')
      # list the files in the directory
      files_directory <- list.files(eqtl_celltype_dir, full.names = F, recursive = F)
      # and filter those by our pattern
      eqtl_files_directory <- files_directory[grepl(eqtl_file_pattern, files_directory)]
      # go through each file
      for (eqtl_file in eqtl_files_directory) {
        # get the full path to the eqtl file
        eqtl_file_full_loc <- paste(eqtl_celltype_dir, '/', eqtl_file, sep = '')
        # read the file
        eqtl_file_full <- fread(eqtl_file_full_loc, header = T, sep = '\t')
        # subset to the variants we have in our caQTLs to reduce search space
        eqtl_file_full <- eqtl_file_full[eqtl_file_full[[qtl_results_variant_column]] %in% caqtls_finemapped[[finemapped_variant_column]], ]
        # now filter on significance
        for (significance_column in names(qtl_resuls_filter_significances)) {
          eqtl_file_full <- eqtl_file_full[eqtl_file_full[[significance_column]] < qtl_resuls_filter_significances[[significance_column]], ]
        }
        # put in the list
        results_each_file_list[[eqtl_file]] <- eqtl_file_full
      }
      # merge outputs of files
      results_each_file <- do.call('rbind', results_each_file_list)
      # extract the unique feature and variant combinations
      qtl_variant_feature <- unique(results_each_file[, c(..qtl_results_variant_column, ..qtl_results_feature_column)])
      # same for the finemapped data
      finemapped_variant_feature <- unique(caqtls_finemapped_celltype[, c(..finemapped_variant_column, ..finemapped_feature_column)])
      # now merge them
      confinement <- merge(finemapped_variant_feature, qtl_variant_feature, by.x = finemapped_variant_column, by.y = qtl_results_variant_column)
      # paste together the output path
      output_path_full <- paste(output_dir, '/', output_file_prepend, cell_type, output_file_append, sep = '')
      # and a zipped one if necessary
      output_path_full_gz <- output_path_full
      # gz file ends with .gz
      if (grepl('.gz$', output_path_full)) {
        # gzip if ends with .gz
        output_path_full_gz <- gzfile(output_path_full)
      }
      write.table(confinement, output_path_full_gz, sep = '\t', row.names = F, col.names = F, quote = F)
      # also make a checksum
      mdfiver::create_md5_for_file(output_path_full)
    }
    # otherwise give a warning
    else {
      warning(paste('cell type', cell_type, 'in finemapping, but no folder found for it in', eqtl_celltype_dir))
    }
  }
  return(0)
}


####################
# Settings        #
####################



####################
# Main Code        #
####################

# location of the finemapped caQTLs
caqtls_finemapped_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output_cis_50kb_updated_features/finemapping/susie_output/aggregated/ALL_Finemapped_UT_snp_vs_peak.txt'
# location of the eQTLs
eqtls_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/UT/'
# where to store the confinements
confinements_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/confinements/'

# read the finemapped results file
caqtls_finemapped_ut <- fread(caqtls_finemapped_ut_loc, header = T, sep = '\t')

# write the confinements
create_mediation_confinements(
  caqtls_finemapped_ut, 
  eqtls_ut_loc, 
  confinements_ut_loc
)