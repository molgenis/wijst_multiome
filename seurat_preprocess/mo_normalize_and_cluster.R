#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_normalize_and_cluster.R
# Function: perform normalization
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)


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
  seurat_object[['MJ']] <- CreateAssay5Object(data = norm_count_matrix)
  return(seurat_object)
}




####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
set.seed(7777)


####################
# Main Code        #
####################

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'

# load object
mo_loc <- paste(seurat_objects_loc, 'mo_all_20240223_seuratv5_filtered.rds', sep = '')
mo <- readRDS(mo_loc)

# do MJ normalization
mo <- normalize_mj(mo)

# do normalization
mo <- SCTransform(mo, vst.flavor = 'v2')
# explicitly set this assay
DefaultAssay(mo) <- 'SCT'
# do clustering pipeline
mo <- FindVariableFeatures(mo)
mo <- RunPCA(mo)
mo <- RunUMAP(mo, dims = 1:30, return.model = T)
mo <- FindNeighbors(mo, dims = 1:30)
mo <- FindClusters(mo, resolution = 1.5)

# save the result
mo_normalized_loc <- paste(seurat_objects_loc, 'mo_all_20240223_seuratv5_normalized.rds', sep = '')
saveRDS(mo, mo_normalized_loc)
