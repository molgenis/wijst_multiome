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


# add metadata that is based on existing incomplete metadata in the seurat object
add_imputed_meta_data <- function(seurat_object, column_to_transform, column_to_reference, column_to_create){
  # add the column
  seurat_object@meta.data[column_to_create] <- NA
  # go through the grouping we have for the entire object
  for(group in unique(seurat_object@meta.data[[column_to_transform]])){
    # subset to get only this group
    seurat_group <- seurat_object[,seurat_object@meta.data[[column_to_transform]] == group]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(seurat_group@meta.data[[column_to_reference]])){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(seurat_group@meta.data[seurat_group@meta.data[[column_to_reference]] == reference & is.na(seurat_group@meta.data[[column_to_reference]]) == F,])
        correctpercent <- number_of_reference_in_group/ncol(seurat_group)
        print(paste(group,"matches",reference,correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    seurat_object@meta.data[seurat_object@meta.data[[column_to_transform]] == group,][column_to_create] <- best_group
    # force cleanup
    rm(seurat_group)
  }
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

# try to impute the stimulation condition
mo <- add_imputed_meta_data(seurat_object = mo, column_to_transform = 'seurat_clusters', column_to_reference = 'condition_final', column_to_create = 'condition_imputed')
# and celltypes where missing
mo <- add_imputed_meta_data(seurat_object = mo, column_to_transform = 'seurat_clusters', column_to_reference = 'predicted.mo_10x_cell_type', column_to_create = 'celltype_imputed')
mo <- add_imputed_meta_data(seurat_object = mo, column_to_transform = 'seurat_clusters', column_to_reference = 'predicted.mo_10x_cell_type.lowerres', column_to_create = 'celltype_imputed_lowerres')

# save the result
mo_normalized_loc <- paste(seurat_objects_loc, 'mo_all_20240223_seuratv5_normalized.rds', sep = '')
saveRDS(mo, mo_normalized_loc)
