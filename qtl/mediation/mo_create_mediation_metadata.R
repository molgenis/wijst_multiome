#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_mediation_metadata.R
# Function: create metadata for the mediation analysis
############################################################################################################################

####################
# libraries        #
####################



####################
# Functions        #
####################

#' Add Covariates to Metadata for Each Cell Type
#'
#' This function adds specified covariates to the metadata of each cell type. It reads covariate data from files, matches them to the metadata based on sample identifiers, and appends the covariates to the metadata.
#'
#' @param metadata_per_celltype A list where each element is a data frame containing metadata for a specific cell type. The names of the list elements should correspond to cell type names.
#' @param covariates_loc A string specifying the directory location where the covariate files are stored.
#' @param covariates_to_take A vector of strings specifying the names of the covariates to be added. Default is `paste('PC', 1:10, sep='')`, which corresponds to `PC1` to `PC10`.
#' @param sample_column_metata A string specifying the column name in the metadata that contains the sample identifiers. Default is `'sample'`.
#' @param sample_column_covariates A string specifying the column name in the covariates file that contains the sample identifiers. If `NULL`, the row names of the covariates file are used. Default is `NULL`.
#' @param covariates_prepend_to_add A string to prepend to the names of the covariates when adding them to the metadata. Default is an empty string.
#' @param covariates_file_prepend A string to prepend to the covariate file names. Default is an empty string.
#' @param covariates_file_append A string to append to the covariate file names. Default is `'.qtlInput.Pcs.txt.gz'`.
#' @return A list with the updated metadata for each cell type, including the added covariates.
#' @examples
#' result <- add_covariates_each_celltype(metadata_per_celltype, covariates_loc)
add_covariates_each_celltype <- function(metadata_per_celltype, covariates_loc, covariates_to_take=paste('PC', 1:10, sep=''), sample_column_metata='sample', sample_column_covariates=NULL, covariates_prepend_to_add='', covariates_file_prepend='', covariates_file_append='.qtlInput.Pcs.txt.gz') {
  # check each cell type
  for (cell_type in names(metadata_per_celltype)) {
    # get that metadata table
    metadata_celltype <- metadata_per_celltype[[cell_type]]
    # now get the path to the covariates file
    covariates_path <- paste(covariates_loc, covariates_file_prepend, cell_type, covariates_file_append, sep = '')
    # check if the file is there
    if (file.exists(covariates_path)) {
      # read the file
      covariates <- NULL
      # if the sample column is NULL, it must be the rownames
      if (is.null(sample_column_covariates)) {
        covariates_input <- read.table(covariates_path, header = T, sep = '\t', row.names = 1)
        # get the matching covariates
        covariates <- covariates_input[match(metadata_celltype[[sample_column_metata]], rownames(covariates_input)), covariates_to_take]
      }
      else {
        covariates_input <- read.table(covariates_path, header = T, sep = '\t')
        # get the matching covariates
        covariates <- covariates_input[match(metadata_celltype[[sample_column_metata]], covariates_input[[sample_column_covariates]]), covariates_to_take]
      }
      # rename the covariates with the prepend
      colnames(covariates) <- paste(covariates_prepend_to_add, colnames(covariates), sep = '')
      # add to the metadata
      metadata_celltype <- cbind(metadata_celltype, covariates)
      # put back into list
      metadata_per_celltype[[cell_type]] <-metadata_celltype
    }
    else {
      warning(paste('no covariates for', cell_type, 'at', covariates_path, ', no covariates added!'))
    }
  }
  return(metadata_per_celltype)
}



####################
# Settings        #
####################



####################
# Main Code        #
####################

# the seurat object
seurat_object_metadata_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
# read that metadata
seurat_object_metadata <- read.table(seurat_object_metadata_loc, header = T, sep = '\t')
# add the lane and sample combination
seurat_object_metadata[['sample']] <- paste(seurat_object_metadata[['sample_final']], seurat_object_metadata[['lane']], sep = ';;')
# subset to all metadata that is not cell-specific
seurat_object_metadata_aggregated <- unique(seurat_object_metadata[, c('sample', 'lane', 'sample_final', 'condition_final', 'age', 'sex', 'LONG_COVID_final', 'celltype_imputed_lowerres')])
# rename the columns now
colnames(seurat_object_metadata_aggregated) <- c('sample', 'lane', 'donor', 'condition', 'age', 'sex', 'LONG_COVID', 'celltype')
# make into a list with the metadata per cell type
metadata_per_celltype <- list()
for (cell_type in unique(seurat_object_metadata_aggregated[['celltype']])) {
  metadata_per_celltype[[cell_type]] <- seurat_object_metadata_aggregated[seurat_object_metadata_aggregated[['celltype']] == cell_type, ]
}

# the location of the joint expression PCs
pcs_exp_joint_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/input/L1/combined/'
# the location of the joint atac PCs
pcs_atac_joint_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/input/L1/combined/'
# the location of the UT and 24hCA expression PCs separately
pcs_exp_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/L1/UT/'
pcs_exp_24hca_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/L1/24hCA/'
# the location of the UT and 24hCA atac PCs separately
pcs_atac_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/'
pcs_atac_24hca_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/24hCA/'

# add each set of covariates
metadata_per_celltype <- add_covariates_each_celltype(metadata_per_celltype, pcs_exp_joint_loc, covariates_prepend_to_add = 'RNA_J_')
metadata_per_celltype <- add_covariates_each_celltype(metadata_per_celltype, pcs_atac_joint_loc, covariates_prepend_to_add = 'ATAC_J_')
metadata_per_celltype <- add_covariates_each_celltype(metadata_per_celltype, pcs_exp_ut_loc, covariates_prepend_to_add = 'RNA_UT_')
metadata_per_celltype <- add_covariates_each_celltype(metadata_per_celltype, pcs_exp_24hca_loc, covariates_prepend_to_add = 'RNA_24hCA_')
metadata_per_celltype <- add_covariates_each_celltype(metadata_per_celltype, pcs_atac_ut_loc, covariates_prepend_to_add = 'ATAC_UT_', covariates_file_append = '.qtlInput.Pcs.txt')
metadata_per_celltype <- add_covariates_each_celltype(metadata_per_celltype, pcs_atac_24hca_loc, covariates_prepend_to_add = 'ATAC_24hCA_', covariates_file_append = '.qtlInput.Pcs.txt')

# write each set of metadata files
metadata_tables_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/metadata/'
for (cell_type in names(metadata_per_celltype)) {
  write.table(
    metadata_per_celltype[[cell_type]],
    gzfile(paste(metadata_tables_loc, cell_type, '.metadata.tsv.gz', sep = '')),
    row.names = F, col.names = T, sep = '\t'
  )
}
