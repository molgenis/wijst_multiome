#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_filter_qc.R
# Function: filter the cells
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


ref10xmo_predictions_to_lower_res_mapping <- function() {
  high_to_low <- list()
  high_to_low[['CD4 Naive']] <- 'CD4T'
  high_to_low[['CD4 TCM']] <- 'CD4T'
  high_to_low[['CD8 Naive']] <- 'CD8T'
  high_to_low[['CD16 Mono']] <- 'monocyte'
  high_to_low[['NK']] <- 'NK'
  high_to_low[['Treg']] <- 'T_other'
  high_to_low[['CD14 Mono']] <- 'monocyte'
  high_to_low[['cDC']] <- 'DC'
  high_to_low[['CD8 TEM_1']] <- 'CD8T'
  high_to_low[['Intermediate B']] <- 'B'
  high_to_low[['Naive B']] <- 'B'
  high_to_low[['Plasma']] <- 'plasmablast'
  high_to_low[['CD4 TEM']] <- 'CD4T'
  high_to_low[['MAIT']] <- 'T_other'
  high_to_low[['Memory B']] <- 'B'
  high_to_low[['gdT']] <- 'T_other'
  high_to_low[['pDC']] <- 'DC'
  high_to_low[['CD8 TEM_2']] <- 'CD8T'
  high_to_low[['HSPC']] <- 'hemapoietic_stem'
  high_to_low[['unannotated']] <- 'unannotated'
  return(high_to_low)
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
mo_loc <- paste(seurat_objects_loc, 'mo_all_20240214_seuratv5_assigned.rds', sep = '')
mo <- readRDS(mo_loc)
ncol(mo)
# [1] 1423526

# remove doublets
mo <- mo[, !is.na(mo@meta.data[['soup_status']]) & mo@meta.data[['soup_status']] == 'singlet']
ncol(mo)
# [1] 1037166

# get location of the cell type assignments
cell_type_assignments_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cell_type_assignment/azimuth/10x_multiome_PBMCs/mo_azimuth_ct_10xmultiome.tsv'
cell_type_assignments <- read.table(cell_type_assignments_loc, sep = '\t', header = T, row.names = 1)

# add this data as well
mo <- AddMetaData(mo, cell_type_assignments[, colnames(cell_type_assignments)])

# add MJ normalization
mo <- normalize_mj(mo)

# remove lower than 200 UMIs
mo <- mo[, mo@meta.data[['nCount_RNA']] >= 200]
ncol(mo)
# [1] 945298

# ncount outliers at max
mo <- mo[, !is.na(mo@meta.data[['nCount_RNA_mad']]) & mo@meta.data[['nCount_RNA_mad']] == 'NotOutlier']
ncol(mo)
# [1] 924637

# nfeature outliers at max
mo <- mo[, !is.na(mo@meta.data[['nFeature_RNA_mad']]) & mo@meta.data[['nFeature_RNA_mad']] == 'NotOutlier']
ncol(mo)
# [1] 924636

# set NA to unannotated for cell types
mo@meta.data[is.na(mo@meta.data[['predicted.mo_10x_cell_type']]), 'predicted.mo_10x_cell_type'] <- 'unannotated'

# add lower res cell type classification
mo@meta.data[['predicted.mo_10x_cell_type.lowerres']] <- as.vector(unlist(ref10xmo_predictions_to_lower_res_mapping()[as.character(mo@meta.data[['predicted.mo_10x_cell_type']])]))

# write the annotated and filtered file
mo_filtered_loc <- paste(seurat_objects_loc, 'mo_all_20240214_seuratv5_filtered.rds', sep = '')
saveRDS(mo, mo_filtered_loc)
