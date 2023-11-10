#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cluster_seurat_objects.R
# Function: merge the seurat objects
############################################################################################################################

####################
# libraries        #
####################

# required to create object
library(Seurat)


####################
# Functions        #
####################


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')


####################
# Main Code        #
####################

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/'

# get the object
object_loc_souped <- paste(seurat_objects_loc, 'mo_all_souped_20231109.rds', sep = '')
object_all <- readRDS(object_loc_souped)

# merge layers
object_all <- JoinLayers(object_all)

# do normalization
object_all <- NormalizeData(object_all)
object_all <- ScaleData(object_all)

# do PCA
object_all <- FindVariableFeatures(object_all, layer = 'data')
object_all <- RunPCA(object_all)
# set seed
set.seed(1337)
# umap
object_all <- RunUMAP(object_all, dims = 1:30, return.model = T)
# knn
object_all <- FindNeighbors(object_all, dims = 1:30)
# find clusters
object_all <- FindClusters(object_all, resolution = 1.2)

# save result
object_all_cluster_unfiltered_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_unfiltered_20231109.rds', sep = '')
saveRDS(object_all, object_all_cluster_unfiltered_loc)

# subset to only singlets
nrow(object_all@meta.data)
# 1143976
object_all <- object_all[, !is.na(object_all@meta.data[['soup_status']]) & object_all@meta.data[['soup_status']] == 'singlet']
nrow(object_all@meta.data)
# 817661

# subset by the nCount MAD
object_all <- object_all[, !is.na(object_all@meta.data[['nCount_RNA_mad']]) & object_all@meta.data[['nCount_RNA_mad']] == "NotOutlier"]
nrow(object_all@meta.data)
# 790355

# subset by nFeature
object_all <- object_all[, !is.na(object_all@meta.data[['nFeature_RNA_mad']]) & object_all@meta.data[['nFeature_RNA_mad']] == "NotOutlier"]
nrow(object_all@meta.data)

# backup the old clusters and reductions
object_all@meta.data[['RNA_snn_res.1.2_unfiltered']] <- object_all@meta.data[['RNA_snn_res.1.2']]
object_all[['umap_unfiltered']] <- object_all[['umap']]

# redo clustering pipeline
object_all <- ScaleData(object_all)
object_all <- FindVariableFeatures(object_all, layer = 'data')
object_all <- RunPCA(object_all)
object_all <- RunUMAP(object_all, dims = 1:30, return.model = T)
object_all <- FindNeighbors(object_all, dims = 1:30)
object_all <- FindClusters(object_all, resolution = 1.2)

# save result
object_all_cluster_filtered_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_filtered_20231109.rds', sep = '')
saveRDS(object_all, object_all_cluster_filtered_loc)
