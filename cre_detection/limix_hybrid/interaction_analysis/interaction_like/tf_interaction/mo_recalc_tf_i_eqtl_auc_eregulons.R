#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_recalc_tf_i_eqtl_auc_eregulons.R
# Function: recalculate the AUC excluding the eGenes without the TF-i-eGene in them to see if that is a confounding effect
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

# location of the eregulon from scenics
scenic_eregs_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_signatures.tsv.gz'
# read that
scenic_eregs <- fread(scenic_eregs_loc, header = T, sep = '\t')
# check each of the gene based ones
for (ereg in unique(scenic_eregs[scenic_eregs[['modality']] == 'Gene_based', ][['signature_name']])) {
  # get the genes
  genes_ereg <- unique(scenic_eregs[
    scenic_eregs[['signature_name']] == ereg,
  ][['gene_or_region']])
  # put in the list
  ereg_to_genes[[ereg]] <- unique(genes_ereg)
}

# location of the Seurat object
seurat_object_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240223_seuratv5_normalized.rds'
# laod the data
seurat_object <- readRDS(seurat_object_loc)
# extract expression matrix
expr_mat <- seurat_object@assays$RNA@layers$counts
# set the barcodes as column names
colnames(expr_mat) <- colnames(seurat_object)
# extract genes
gene_names <- data.frame(seurat_object@assays$RNA@features)
gene_names_counts <- rownames(gene_names[gene_names[['counts']], , drop = F])
# set those as row names
rownames(expr_mat) <- gene_names_counts
# subset the expression matrix
expr_mat_included <- expr_mat[, colnames(expr_mat) %in% barcodes]

# read the TF-interaction-QTL output
tf_i_eqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls/merged/results_fdr.tsv.gz'
tf_i_eqtl <- fread(tf_i_eqtl_loc,  header = T, sep = '\t')
# get what is significant
sc_tf_ieqtl_sig <- tf_i_eqtl[
  tf_i_eqtl[['region:genotype_bh']] < 0.05 & 
    tf_i_eqtl[['anova_bh']] < 0.05 & 
    tf_i_eqtl[['region_bh']] < 0.05, # & 
    # tf_i_eqtl[['genotype_bh']] < 0.05, 
]
# remove the top gene from the eregulon lists
ereg_to_genes_iegenes <- list()
for (ereg in names(ereg_to_genes)) {
  # check if there is an interaction
  if (ereg %in% sc_tf_ieqtl_sig[['region']]) {
    # get the genes
    genes_ereg <- ereg_to_genes[[ereg]]
    # get the genes
    tf_i_egenes <- unlist(as.vector(sc_tf_ieqtl_sig_top_tf_per_gene[sc_tf_ieqtl_sig_top_tf_per_gene[['region']] == ereg, ][['gene']]))
    # check each gene
    for (tf_i_egene in tf_i_egenes) {
      # and thate one
      genes_ereg_no_tf_i_egene <- setdiff(genes_ereg, tf_i_egene)
      # put in the list
      ereg_to_genes_iegenes[[paste(ereg, tf_i_egene, sep = '_')]] <- genes_ereg_no_tf_i_egene
    }
  }
}
# recalculate enrichment scores
cells_AUC_notop <- AUCell_run(expr_mat_included, ereg_to_genes_iegenes)
# extract the auc matrix
cells_AUC_notop_mtx <- getAUC(cells_AUC_notop)
# as matrix
cells_AUC_notop_mtx <- as.matrix(cells_AUC_notop_mtx)
# extract dimension names
cells_AUC_notop_eregs <- rownames(cells_AUC_notop_mtx)
cells_AUC_notop_barcodes <- colnames(cells_AUC_notop_mtx)
# make into data.table
cells_AUC_notop_dt <- data.table(cells_AUC_notop_mtx)
# but add the eregulon as the first column
cells_AUC_notop_dt <- cbind(data.table('eregulon' = cells_AUC_notop_eregs), cells_AUC_notop_dt)
# set export location
cells_AUC_notop_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/mo_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed.tsv.gz'
# write the output
write.table(cells_AUC_notop_dt, gzfile(cells_AUC_notop_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# with a checksum
mdfiver::create_sha256_for_file(cells_AUC_notop_loc)

# make a new confinement file
sc_tf_ieqtl_sig_confinement <- sc_tf_ieqtl_sig[, c('variant', 'region', 'gene')]
# and replace the region with the eregulon name and the gene, as we generated for the new AUC matrix
sc_tf_ieqtl_sig_confinement[['region']] <- paste(sc_tf_ieqtl_sig[['region']], sc_tf_ieqtl_sig[['gene']], sep = '_')
# update the column names to actually reflect what we do
colnames(sc_tf_ieqtl_sig_confinement) <- c('variant', 'eregulon', 'gene')
# set output loc
sc_tf_ieqtl_sig_confinement_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement_inclcaqtls_significant.tsv.gz'
# write output
write.table(sc_tf_ieqtl_sig_confinement, gzfile(sc_tf_ieqtl_sig_confinement_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# with a checksum
mdfiver::create_sha256_for_file(sc_tf_ieqtl_sig_confinement_loc)
