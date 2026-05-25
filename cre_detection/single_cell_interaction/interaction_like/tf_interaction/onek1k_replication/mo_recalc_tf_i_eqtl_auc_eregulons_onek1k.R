#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_recalc_tf_i_eqtl_auc_eregulons_onek1k.R
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
library(R.utils)

####################
# Functions        #
####################


####################
# Main code        #
####################

# location of onek1k cell-type objects
onek1k_rds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/input/L1/'
# which cell types to use
onek1k_cts <- c('B', 'CD4_T', 'CD8_T', 'DC', 'Mono', 'NK', 'other_T')
# append for each file
onek1k_rds_append <- '.Qced.Normalized.SCs.Rds'
# make vector of objects
onek1k_object_list <- list()
# read each cell type
for (i in 1 : length(onek1k_cts)) {
  # get cell type
  cell_type<- onek1k_cts[i]
  # paste object location together
  onek1k_ct_object_loc <- paste0(onek1k_rds_loc, '/', cell_type, onek1k_rds_append)
  # read that object
  onek1k_ct_object <- readRDS(onek1k_ct_object_loc)
  # place in the list
  onek1k_object_list[[i]] <- onek1k_ct_object
}
# merge all
onek1k_object <- merge(onek1k_object_list[[1]], as.vector(onek1k_object_list[2:length(onek1k_object_list)]), add.cell.ids = onek1k_cts)
# upgrade to latest Seurat
onek1k_object <- UpdateSeuratObject(onek1k_object)
# set location to save object
onek1k_merged_major_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/onek1k/seurat_objects/onek1k_major_cts.rds'
# save the object
saveRDS(onek1k_object, onek1k_merged_major_loc)
# create checksum
mdfiver::create_sha256_for_file(onek1k_merged_major_loc)
# do SCT normalization
onek1k_object <- readRDS(onek1k_merged_major_loc)
#onek1k_object <- SCTransform(onek1k_object)
# extract the expression matrix
expr_mat <- onek1k_object@assays$RNA@counts
#expr_mat <- onek1k_object@assays$SCT@counts
# set the barcodes as column names
#colnames(expr_mat) <- colnames(onek1k_object)
# extract genes
#gene_names <- data.frame(onek1k_object@assays$SCT@features)
#gene_names_counts <- rownames(gene_names[gene_names[['counts']], , drop = F])

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

# location of the eregulon from scenics
scenic_eregs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_signatures.tsv.gz'
# read that
scenic_eregs <- fread(scenic_eregs_loc, header = T, sep = '\t')
# store in a list
ereg_to_genes <- list()
# check each of the gene based ones
for (ereg in unique(scenic_eregs[scenic_eregs[['modality']] == 'Gene_based', ][['signature_name']])) {
  # get the genes
  genes_ereg <- unique(scenic_eregs[
    scenic_eregs[['signature_name']] == ereg,
  ][['gene_or_region']])
  # put in the list
  ereg_to_genes[[ereg]] <- unique(genes_ereg)
}


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
cells_AUC_notop <- AUCell_run(expr_mat, ereg_to_genes_iegenes)
# save result somewhere
cells_AUC_notop_loc <- '~/multiome/rds/onek1k_auc_mtx_noiegene_fromscenic.rds'
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
cells_AUC_notop_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/onek1k//tf_interaction/onek1k_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed.tsv.gz'
# write the output
write.table(cells_AUC_notop_dt, gzfile(cells_AUC_notop_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# with a checksum
mdfiver::create_sha256_for_file(cells_AUC_notop_loc)

# write cells
cells_AUC_notop_barcodes_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/onek1k//tf_interaction/onek1k_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed_barcodes.txt.gz'
# write the barcode names
write.table(data.frame(x = cells_AUC_notop_barcodes), gzfile(cells_AUC_notop_barcodes_loc), row.names = F, col.names = F, quote = F)
# checksum
mdfiver::create_sha256_for_file(cells_AUC_notop_barcodes_loc)
# chunk  max of 1000 eregs at a time
chunk_i_start <- 1
while(chunk_i_start <= nrow(cells_AUC_notop_mtx)) {
  # get the end
  chunk_i_end <- chunk_i_start + 999
  # except if bigger than the end
  if (chunk_i_end > nrow(cells_AUC_notop_mtx)) {
    chunk_i_end <- nrow(cells_AUC_notop_mtx)
  }
  # write matrix as well
  cells_AUC_notop_mtx_loc <- paste0('/groups/umcg-franke-scrna/tmp04/external_datasets/onek1k//tf_interaction/onek1k_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed_', chunk_i_start, '_' , chunk_i_end, '.mtx')
  cells_AUC_notop_eregs_loc <- paste0('/groups/umcg-franke-scrna/tmp04/external_datasets/onek1k//tf_interaction/onek1k_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed_', chunk_i_start, '_' , chunk_i_end, '_ergenames.txt.gz')
  # write the matrix
  writeMM(Matrix(cells_AUC_notop_mtx[chunk_i_start:chunk_i_end, ], sparse = TRUE), cells_AUC_notop_mtx_loc)
  # zip the file
  gzip(filename = cells_AUC_notop_mtx_loc)
  # make checksum
  mdfiver::create_sha256_for_file(paste0(cells_AUC_notop_mtx_loc, '.gz'))
  # write the ereg names
  write.table(data.frame(x = cells_AUC_notop_eregs[chunk_i_start:chunk_i_end]), gzfile(cells_AUC_notop_eregs_loc), row.names = F, col.names = F, quote = F)
  # checksum
  mdfiver::create_sha256_for_file(cells_AUC_notop_eregs_loc)
  
  # increase count
  chunk_i_start <- chunk_i_start + 1000
}
