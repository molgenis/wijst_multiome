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
  unique(genes_ereg <- scenic_output[
    scenic_output[['Gene_signature_name']] == ereg, 
  ][['Gene']])
  # put in the list
  ereg_to_genes[[ereg]] <- genes_ereg
}
# extract expression matrix
expr_mat <- seurat_object@assays$RNA@layers$counts
# set the barcodes as column names
colnames(expr_mat) <- colnames(seurat_object)
# extract genes
gene_names <- data.frame(seurat_object@assays$RNA@features)
gene_names_counts <- rownames(gene_names[gene_names[['counts']], , drop = F])
# set those as row names
rownames(expr_mat) <- gene_names_counts
# Calculate enrichment scores
cells_AUC <- AUCell_run(expr_mat, ereg_to_genes)
# set the path to save this one
cells_AUC_loc <- '~/multiome/rds/mo_auc_mtx_full.rds'
# save the file
saveRDS(cells_AUC, cells_AUC_loc)
# make a checksum
mdfiver::create_sha256_for_file(cells_AUC_loc)

# subset the expression matrix
expr_mat_included <- expr_mat[, colnames(expr_mat) %in% barcodes]
# Calculate enrichment scores
cells_AUC_included <- AUCell_run(expr_mat_included, ereg_to_genes)
# set the path to save this one
cells_AUC_included_loc <- '~/multiome/rds/mo_auc_mtx_included.rds'
# save the file
saveRDS(cells_AUC_included, cells_AUC_included_loc)
# make a checksum
mdfiver::create_sha256_for_file(cells_AUC_included_loc)

# reload data
cells_AUC_loc <- '~/multiome/rds/mo_auc_mtx_full.rds'
cells_AUC_included_loc <- '~/multiome/rds/mo_auc_mtx_included.rds'
cells_AUC <- readRDS(cells_AUC_loc)
cells_AUC_included <- readRDS(cells_AUC_included_loc)

# first do AUC we had done before vs the new full one
barcodes_before_vs_full <- intersect(barcodes, colnames(cells_AUC))
# then subset
cells_AUC_vs_before <- getAUC(cells_AUC[, barcodes_before_vs_full])
cells_before_vs_AUC <- auc_mtx[, barcodes_before_vs_full]
# then do AUC we had before vs the new limited one
barcodes_before_vs_included <- intersect(barcodes, colnames(cells_AUC_included))
# then subset
cells_AUC_included_vs_before <- getAUC(cells_AUC_included[, barcodes_before_vs_included])
cells_before_vs_AUC_included <- auc_mtx[, barcodes_before_vs_included]

# create a list to store each df
auc_cors <- list()
# get correlation for each gene set
for (gene_set in names(ereg_to_genes)) {
  print(gene_set)
  # do the full matrix vs what we had before
  auc_vs_before_cor <- cor(
    x = as.vector(unlist(cells_AUC_vs_before[gene_set, ])), 
    y = as.vector(unlist(cells_before_vs_AUC[gene_set, ]))
  )
  # do the limited matrix vs what we had before
  auc_included_vs_before_cor <- cor(
    x = as.vector(unlist(cells_AUC_included_vs_before[gene_set, ])), 
    y = as.vector(unlist(cells_before_vs_AUC_included[gene_set, ]))
  )
  # make into a df
  cor_df_gs <- data.frame('regulon' = c(gene_set), 'cor_full' = c(auc_vs_before_cor), 'cor_included' = c(auc_included_vs_before_cor))
  # put in the list
  auc_cors[[gene_set]] <- cor_df_gs
}
# merge all
auc_cors_all <- do.call('rbind', auc_cors)

# read the TF-interaction-QTL output
tf_i_eqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls/merged/results_fdr.tsv.gz'
tf_i_eqtl <- fread(tf_i_eqtl_loc,  header = T, sep = '\t')
# get what is significant
sc_tf_ieqtl_sig <- tf_i_eqtl[
  tf_i_eqtl[['region:genotype_bh']] < 0.05 & 
    tf_i_eqtl[['anova_bh']] < 0.05 & 
    tf_i_eqtl[['region_bh']] < 0.05 & 
    tf_i_eqtl[['genotype_bh']] < 0.05, 
]
# get the top TF for each gene
sc_tf_ieqtl_sig_top_tf_per_gene <- sc_tf_ieqtl_sig[order(sc_tf_ieqtl_sig[['region:genotype_p']]), ]
sc_tf_ieqtl_sig_top_tf_per_gene <- sc_tf_ieqtl_sig_top_tf_per_gene[!duplicated(sc_tf_ieqtl_sig_top_tf_per_gene[['gene']]), ]
# remove the top gene from the eregulon lists
ereg_to_genes_notop <- list()
for (ereg in names(ereg_to_genes)) {
  # check if there is an interaction
  if (ereg %in% sc_tf_ieqtl_sig_top_tf_per_gene[['region']]) {
    # get the genes
    genes_ereg <- ereg_to_genes[[ereg]]
    # get the top gene
    top_gene <- unlist(as.vector(sc_tf_ieqtl_sig_top_tf_per_gene[sc_tf_ieqtl_sig_top_tf_per_gene[['region']] == ereg, ][['gene']][[1]]))
    # and remove one
    genes_ereg_notop <- setdiff(genes_ereg, top_gene)
    # put in the list
    ereg_to_genes_notop[[ereg]] <- genes_ereg_notop
  }
}

# Calculate enrichment scores
cells_AUC_notop <- AUCell_run(expr_mat, ereg_to_genes_notop)

# first do AUC we had done before vs the new full one
barcodes_full_vs_notop <- intersect(colnames(cells_AUC), colnames(cells_AUC_notop))
# then subset
cells_AUC_notop_vs_full <- getAUC(cells_AUC_notop[, barcodes_full_vs_notop])
cells_AUC_full_vs_notop <- getAUC(cells_AUC[, barcodes_full_vs_notop])
