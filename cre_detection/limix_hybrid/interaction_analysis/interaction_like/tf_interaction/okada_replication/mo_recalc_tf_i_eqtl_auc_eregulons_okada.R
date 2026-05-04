#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_recalc_tf_i_eqtl_auc_eregulons_okada.R
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

# location of okada cell-type objects
okada_rds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_okada/input/L1/'
# which cell types to use
okada_cts <- c('B', 'CD4_T', 'CD8_T', 'DC', 'Mono', 'NK', 'other_T')
# append for each file
okada_rds_append <- '.Qced.Normalized.SCs.Rds'
# make vector of objects
okada_object_list <- list()
# read each cell type
for (i in 1 : length(okada_cts)) {
  # get cell type
  cell_type<- okada_cts[i]
  # paste object location together
  okada_ct_object_loc <- paste0(okada_rds_loc, '/', cell_type, okada_rds_append)
  # read that object
  okada_ct_object <- readRDS(okada_ct_object_loc)
  # place in the list
  okada_object_list[[i]] <- okada_ct_object
}
# merge all
okada_object <- merge(okada_object_list[[1]], as.vector(okada_object_list[2:length(okada_object_list)]), add.cell.ids = okada_cts)
# upgrade to latest Seurat
okada_object <- UpdateSeuratObject(okada_object)
# set location to save object
okada_merged_major_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/okada/seurat_objects/okada_major_cts.rds'
# save the object
saveRDS(okada_object, okada_merged_major_loc)
# create checksum
mdfiver::create_sha256_for_file(okada_merged_major_loc)
# do SCT normalization
okada_object <- SCTransform(okada_object)
# extract the expression matrix
# expr_mat <- seurat_object@assays$RNA@layers$counts
expr_mat <- okada_object@assays$SCT@counts
# set the barcodes as column names
colnames(expr_mat) <- colnames(okada_object)
# extract genes
gene_names <- data.frame(okada_object@assays$SCT@features)
gene_names_counts <- rownames(gene_names[gene_names[['counts']], , drop = F])


# location of the SCENIC outout
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both_filtered.tsv.gz'
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

# read the TF-interaction-QTL output
tf_i_eqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz'
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
    tf_i_egenes <- unique(unlist(as.vector(sc_tf_ieqtl_sig[sc_tf_ieqtl_sig[['region']] == ereg, ][['gene']])))
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
# save result somewhere
cells_AUC_notop_loc <- '~/multiome/rds/okada_auc_mtx_noiegene_fromscenic.rds'
saveRDS(cells_AUC_notop, cells_AUC_notop_loc)
mdfiver::create_sha256_for_file(cells_AUC_notop_loc)
# extract the auc matrix
cells_AUC_notop_mtx <- getAUC(cells_AUC_notop)
# as matrix
cells_AUC_notop_mtx <- as.matrix(cells_AUC_notop_mtx)
# extract dimension names
cells_AUC_notop_eregs <- rownames(cells_AUC_notop_mtx)
cells_AUC_notop_barcodes <- colnames(cells_AUC_notop_mtx)
# make into data.table
cells_AUC_notop_dt <- data.table(cells_AUC_notop_mtx)
# extract the repressors
cells_AUC_notop_eregs_repressors_i <- grep('-/+', cells_AUC_notop_eregs)
# subset both the eregs and the dt on those
cells_AUC_notop_eregs_repressors <- cells_AUC_notop_eregs[cells_AUC_notop_eregs_repressors_i]
cells_AUC_notop_dt_repressors <- cells_AUC_notop_dt[cells_AUC_notop_eregs_repressors_i, ]
# and modify
# cells_AUC_notop_dt_repressors <- 1 - cells_AUC_notop_dt_repressors
cells_AUC_notop_dt_repressors <- -1 * cells_AUC_notop_dt_repressors
# take the originals as well
cells_AUC_notop_eregs_activators <- cells_AUC_notop_eregs[-cells_AUC_notop_eregs_repressors_i]
cells_AUC_notop_dt_activators <- cells_AUC_notop_dt[-cells_AUC_notop_eregs_repressors_i, ]
# merge them again
cells_AUC_notop_dt <- rbind(cells_AUC_notop_dt_repressors, cells_AUC_notop_dt_activators)
cells_AUC_notop_eregs <- c(cells_AUC_notop_eregs_repressors, cells_AUC_notop_eregs_activators)
# but add the eregulon as the first column
cells_AUC_notop_dt <- cbind(data.table('eregulon' = cells_AUC_notop_eregs), cells_AUC_notop_dt)
# set export location
cells_AUC_notop_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/okada//tf_interaction/okada_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed.tsv.gz'
# write the output
write.table(cells_AUC_notop_dt, gzfile(cells_AUC_notop_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# with a checksum
mdfiver::create_sha256_for_file(cells_AUC_notop_loc)


