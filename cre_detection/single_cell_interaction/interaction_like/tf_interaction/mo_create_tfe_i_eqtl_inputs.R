#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_tfe_i_eqtl_inputs.R
# Function: create input files for TF-i-eQTL mapping, but with expression of TF instead of activity
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(Seurat)
library(Matrix)


####################
# Functions        #
####################


normalize_mj <- function(seurat_object) {
  # get the count matrix where we have the correct cell type
  count_matrix <- GetAssayData(seurat_object, slot = "counts")
  # ignore genes that are never expressed
  count_matrix <-  count_matrix[which(rowSums(count_matrix) != 0), ]
  # create new object to store the counts in
  norm_count_matrix <- count_matrix
  # do mean sample-sum normalization
  sample_sum_info = colSums(norm_count_matrix)
  mean_sample_sum = mean(sample_sum_info)
  sample_scale = sample_sum_info / mean_sample_sum
  # divide each column by sample_scale
  norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
  if ('layers' %in% slotNames(seurat_object[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    seurat_object[['MJ']] <- CreateAssay5Object(data = norm_count_matrix)
    
  } else {
    print('using Seurat v3/4 style \'slot\'')
    seurat_object[['MJ']] <- CreateAssayObject(data = norm_count_matrix)
  }
  return(seurat_object)
}


#############
# Main code #
#############

# location of the Seurat object
seurat_object_loc <- paste0('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds')
# read the object
seurat_object <- readRDS(seurat_object_loc)
# add MJ normalization
seurat_object <- normalize_mj(seurat_object)

# location of the AUC matrix
auc_mtx_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc.mtx.gz'
# location of the barcodes
auc_barcodes_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_barcodes.txt.gz'
# and the eregulon names
ereg_names_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_eregnames.txt.gz'
# read the matrix
auc_mtx <- Matrix::readMM(auc_mtx_loc)
# read the barcodes and eregnames
barcodes <- fread(auc_barcodes_loc, header = F)[[1]]
eregnames <- fread(ereg_names_loc, header = F)[[1]]
# construct matrix in Seurat format
colnames(auc_mtx) <- barcodes
rownames(auc_mtx) <- eregnames
auc_assay <- CreateAssay5Object(data = auc_mtx)

# keep in the seurat object only what we have in the AUC matrix
seurat_object <- seurat_object[, barcodes]

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
# rename the regions
scenic_output[['region_cpeaks']] <- gsub(':', '-', scenic_output[['Region']])

# create a mapping
mapping <- data.frame(
  'TF' = gsub("_(extended|direct).*", "", eregnames), 
  'eRegulon' = eregnames
)

# check which of the TFs are in the expression data for the slot I care about
features_mj_data <- data.frame(seurat_object@assays$MJ@features)
# get which TFs are in the features of the MJ assay and in the mapping
tfs_in_mj_data <- intersect(mapping[['TF']], rownames(features_mj_data[features_mj_data[['data']], , drop = F]))
# make those unique to be safe
tfs_in_mj_data <- unique(tfs_in_mj_data)
# get the indices of these TFs in the features
tf_indices <- match(tfs_in_mj_data, rownames(features_mj_data))

# extract the MJ counts
exp_mtx <- seurat_object@assays$MJ@layers$data[
  tf_indices, 
]

# make into data.table
exp_dt <- cbind(
  data.table('tf' = tfs_in_mj_data), 
  data.table(as.matrix(exp_mtx))
)
# set column names
colnames(exp_dt) <- c('tf', barcodes)
# save the result
tf_exp_dt_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/mo_tf_interaction_expression_all_nonsparse_transposed.tsv.gz'
write.table(exp_dt, gzfile(tf_exp_dt_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(tf_exp_dt_loc)

# read the confinement file
confinement_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement_inclcaqtls_varinregion.tsv.gz'
# read this file
confinement <- fread(confinement_loc, header = T, sep = '\t')
# get the tf for the eregulons
confinement[['eregulon']] <- mapping$TF[match(confinement[['eregulon']], mapping$eRegulon)]
# keep only cases where we were able to map
confinement <- confinement[!is.na(confinement[['eregulon']]), ]
# and unique entries
confinement <- unique(confinement)
# save the result
confinement_loc_out <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/mo_var_tfexpression_gene_confinement_inclcaqtls_varinregion.tsv.gz'
write.table(confinement, gzfile(confinement_loc_out), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(confinement_loc_out)
