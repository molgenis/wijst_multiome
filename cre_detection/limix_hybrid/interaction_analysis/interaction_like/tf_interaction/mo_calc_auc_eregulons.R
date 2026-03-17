#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_calc_auc_eregulons.R
# Function: recalculate the AUC excluding the eGenes
#
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Matrix)
library(data.table)
library(AUCell)


####################
# Functions        #
####################


####################
# Main code        #
####################

# location of the AUC matrix
auc_mtx_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc.mtx.gz'
# location of the barcodes
auc_barcodes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_barcodes.txt.gz'
# and the eregulon names
ereg_names_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_eregnames.txt.gz'

# location of the SCENIC outout
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both_filtered.tsv.gz'
# read that
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# order by the extended
scenic_output <- scenic_output[order(scenic_output[['is_extended']]), ]
# get which are non-extended if it was both extended and non-extended
scenic_eregs <- unique(scenic_output[, c('TF', 'Gene_signature_direction', 'Gene_signature_name', 'source')])
scenic_eregs <-scenic_eregs[order(scenic_eregs[['source']]), ]
scenic_eregs_to_keep <- scenic_eregs[!duplicated(paste(scenic_eregs[['TF']], scenic_eregs[['Gene_signature_direction']])), ]
# then use that to filer
scenic_output <- scenic_output[scenic_output[['Gene_signature_name']] %in% scenic_eregs_to_keep[['Gene_signature_name']], ]
# read the matrix
auc_mtx <- Matrix::readMM(auc_mtx_loc)
# read the barcodes and eregnames
barcodes <- fread(auc_barcodes_loc, header = F)[[1]]
eregnames <- fread(ereg_names_loc, header = F)[[1]]
# keep only the eregulons that we kept in our output
auc_mtx <- auc_mtx[eregnames %in% scenic_output[['Gene_signature_name']], ]
eregnames <- eregnames[eregnames %in% scenic_output[['Gene_signature_name']]]
# construct matrix in Seurat format
colnames(auc_mtx) <- barcodes
rownames(auc_mtx) <- eregnames

# location of the Seurat object
seurat_object_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240223_seuratv5_normalized.rds'
# laod the data
seurat_object <- readRDS(seurat_object_loc)

# make a list of all the genes per ereg
ereg_to_genes <- list()
for (ereg in unique(scenic_output[['Gene_signature_name']])) {
  # get the genes
  genes_ereg <- scenic_output[
    scenic_output[['Gene_signature_name']] == ereg, 
  ][['Gene']]
  # put in the list
  ereg_to_genes[[ereg]] <- genes_ereg
}
# extract expression matrix
expr_mat <- seurat_object@assays$RNA@layers$counts
# set the barcodes as column names
colnames(expr_mat) <- colnames(expr_mat)
# extract genes
gene_names <- data.frame(seurat_object@assays$RNA@features)
gene_names_counts <- rownames(gene_names[gene_names[['counts']], , drop = F])
# set those as row names
rownames(expr_mat) <- gene_names_counts
# Calculate enrichment scores
cells_AUC <- AUCell_run(expr_mat, ereg_to_genes)


    