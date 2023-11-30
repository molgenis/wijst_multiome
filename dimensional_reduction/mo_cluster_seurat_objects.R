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
# for plots
library(ggplot2)
library(cowplot)


####################
# Functions        #
####################

nc2022_predictions_to_lower_res_mapping <- function() {
  high_to_low <- list()
  high_to_low[['reg CD4T']] <- 'CD4T'
  high_to_low[['mono 1']] <- 'monocyte'
  high_to_low[['B']] <- 'B'
  high_to_low[['memory CD8T']] <- 'CD8T'
  high_to_low[['unknown']] <- 'unknown'
  high_to_low[['mono 4']] <- 'monocyte'
  high_to_low[['mono 2']] <- 'monocyte'
  high_to_low[['plasma B']] <- 'plasmablast'
  high_to_low[['mono 3']] <- 'monocyte'
  high_to_low[['naive CD4T']] <- 'CD4T'
  high_to_low[['pDC']] <- 'DC'
  high_to_low[['mDC']] <- 'DC'
  high_to_low[['NKbright']] <- 'NK'
  high_to_low[['hemapoietic stem']] <- 'hemapoietic_stem'
  high_to_low[['NKdim']] <- 'NK'
  high_to_low[['th2 CD4T']] <- 'CD4T'
  high_to_low[['memory CD8T left and naive CD8T right']] <- 'CD8T'
  high_to_low[['th1 CD4T']] <- 'CD4T'
  high_to_low[['megakaryocyte']] <- 'megakaryocyte'
  high_to_low[['naive CD4T transitioning to stim']] <- 'CD4T'
  high_to_low[['naive CD8T']] <- 'CD8T'
  high_to_low[['NK']] <- 'NK'
  high_to_low[['double negative T']] <- 'dnT'
  high_to_low[['T helper']] <- 'T_helper'
  high_to_low[['megakaryocytes']] <- 'megakaryocyte'
  high_to_low[['cyto CD4T']] <- 'CD4T'
  return(high_to_low)
}


# add metadata that is based on existing incomplete metadata in the seurat object
add_imputed_meta_data <- function(seurat_object, column_to_transform, column_to_reference, column_to_create){
  # add the column
  seurat_object@meta.data[[column_to_create]] <- NA
  # go through the grouping we have for the entire object
  for(group in unique(seurat_object@meta.data[!is.na(seurat_object@meta.data[[column_to_transform]]), column_to_transform])){
    # subset to get only this group
    seurat_group <- seurat_object[, !is.na(seurat_object@meta.data[[column_to_transform]]) & seurat_object@meta.data[[column_to_transform]] == group]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(seurat_group@meta.data[[column_to_reference]])){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(seurat_group@meta.data[!(is.na(seurat_group@meta.data[[column_to_reference]])) & seurat_group@meta.data[[column_to_reference]] == reference,])
        correctpercent <- number_of_reference_in_group/ncol(seurat_group)
        print(paste(group,"matches", reference, correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    seurat_object@meta.data[!is.na(seurat_object@meta.data[[column_to_transform]]) & seurat_object@meta.data[[column_to_transform]] == group, column_to_create] <- best_group
    # force cleanup
    rm(seurat_group)
  }
  return(seurat_object)
}


add_conditions <- function(seurat_object, condition_mapping, lane_column_mapping='lane', sample_column_mapping='sample', condition_column_mapping='condition', lane_column_metadata='lane', sample_column_metadata='soup_best_match_sample', condition_column_metadata='condition') {
  # add the metadata column
  seurat_object@meta.data[[condition_column_metadata]] <- NA
  # check each row
  for (row_i in 1 : nrow(condition_mapping)) {
    # get the lane
    lane <- condition_mapping[row_i, lane_column_mapping]
    sample <- condition_mapping[row_i, sample_column_mapping]
    condition <- condition_mapping[row_i, condition_column_mapping]
    if (!is.na(lane) & !is.na(sample)) {
      # set the condition
      seurat_object@meta.data[!is.na(seurat_object@meta.data[[lane_column_metadata]]) & seurat_object@meta.data[[lane_column_metadata]] == lane &
                                !is.na(seurat_object@meta.data[[sample_column_metadata]]) & seurat_object@meta.data[[sample_column_metadata]] == sample
                              , condition_column_metadata] <- condition
    }
  }
  return(seurat_object)
}


get_scrublet_output <- function(scrublet_output_loc) {
  # get all files
  files_scrublet <- list.files(scrublet_output_loc)
  # we'll add them in a list first
  scrublet_per_lane <- list()
  # check each file
  for (scrublet_file in files_scrublet) {
    # read the file
    scrublet_output_lane <- read.table(paste(scrublet_output_loc, '/', scrublet_file, sep = ''), header = T, sep = '\t')
    scrublet_per_lane[[scrublet_file]] <- scrublet_output_lane
  }
  # merge all the lanes
  scrublet_all <- do.call('rbind', scrublet_per_lane)
  return(scrublet_all)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  # stimulations
  color_coding[['UT']] <- 'lightgrey'
  color_coding[['24hCa']] <- 'forestgreen'
  return(color_coding)
}


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')


####################
# Main Code        #
####################

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/seurat_preprocess_samples//objects/'

# get the object
object_loc_souped <- paste(seurat_objects_loc, 'mo_all_souped_20231129_seuratv5.rds', sep = '')
object_all <- readRDS(object_loc_souped)

# merge layers
#object_all <- JoinLayers(object_all)

# do normalization
object_all <- NormalizeData(object_all)
#object_all <- ScaleData(object_all)

# do PCA
# object_all <- FindVariableFeatures(object_all, layer = 'data')
# object_all <- RunPCA(object_all)
# set seed
set.seed(1337)
# umap
# object_all <- RunUMAP(object_all, dims = 1:30, return.model = T)
# knn
# object_all <- FindNeighbors(object_all, dims = 1:30)
# find clusters
# object_all <- FindClusters(object_all, resolution = 1.2)

# save result
# object_all_cluster_unfiltered_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_unfiltered_20231129.rds', sep = '')
# saveRDS(object_all, object_all_cluster_unfiltered_loc)

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
# 790233

# and with minimal expression
object_all <- object_all[, object_all@meta.data[['nCount_RNA']] >= 200 & object_all@meta.data[['nFeature_RNA']] >= 3]
nrow(object_all@meta.data)
# 781748


# backup the old clusters and reductions
object_all@meta.data[['RNA_snn_res.1.2_unfiltered']] <- object_all@meta.data[['RNA_snn_res.1.2']]
object_all[['umap_unfiltered']] <- object_all[['umap']]

# redo clustering pipeline
object_all <- ScaleData(object_all)
object_all <- FindVariableFeatures(object_all, layer = 'data')
object_all <- RunPCA(object_all)
object_all <- RunUMAP(object_all, dims = 1:30, return.model = T)
object_all <- FindNeighbors(object_all, dims = 1:30)
object_all <- FindClusters(object_all, resolution = 1.5)

# save the result
object_all_cluster_filtered_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_filtered_20231129.rds', sep = '')
saveRDS(object_all, object_all_cluster_filtered_loc)

# get location of cell type annotation
cell_type_predictions_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/cell_type_assignment/azimuth/NC2022_v3/mo_azimuth_ct_nc2022_v2.tsv'
# read annotation
cell_type_predictions <- read.table(cell_type_predictions_loc, header = T, sep = '\t', row.names = 1)
# add lower resolution predictions
cell_type_predictions[['cell_type_lowerres']] <- as.vector(unlist(nc2022_predictions_to_lower_res_mapping()[cell_type_predictions[['cell_type']]]))
# change column name for score
colnames(cell_type_predictions)[2] <- 'cell_type_score'
# add to the object
object_all <- AddMetaData(object_all, cell_type_predictions['cell_type'])
object_all <- AddMetaData(object_all, cell_type_predictions['cell_type_score'])
object_all <- AddMetaData(object_all, cell_type_predictions['cell_type_lowerres'])

# add imputed annotation for cell types
object_all <- add_imputed_meta_data(object_all, column_to_transform = 'seurat_clusters', column_to_reference = 'cell_type_lowerres', column_to_create = 'cell_type_lowerres_imputed')

# set the NA values
object_all@meta.data[is.na(object_all@meta.data[['cell_type']]), 'cell_type'] <- 'unmapped'

# the output of scrublet
scrublet_output_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/demultiplexing/scrublet/scrublet_output/corrected/'
scrublet_output <- get_scrublet_output(scrublet_output_loc)
# set the lane barcode as the column name
rownames(scrublet_output) <- scrublet_output[['lane_barcode']]
# rename columns
colnames(scrublet_output) <- paste('scrublet_', colnames(scrublet_output), sep = '')
# add the ones we care about
object_all <- AddMetaData(object_all, metadata = scrublet_output[, c('scrublet_doublet', 'scrublet_doublet_score')])

# get the condition mapping
condition_mapping_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'
condition_mapping <- read.table(condition_mapping_loc, header = T, sep = '\t')
# add this information
object_all <- add_conditions(object_all, condition_mapping)
# but only where we are sure of the sample assignment and thus the condition
object_all@meta.data[!is.na(object_all@meta.data$soup_best_match_correlation) & object_all@meta.data$soup_best_match_correlation > 0.7, 'filtered_condition'] <- object_all@meta.data[!is.na(object_all@meta.data$soup_best_match_correlation) & object_all@meta.data$soup_best_match_correlation > 0.7, 'condition']
# we will add an imputed version as well
object_all <- add_imputed_meta_data(object_all, column_to_transform = 'seurat_clusters', column_to_reference = 'condition', column_to_create = 'condition_imputed')
object_all <- add_imputed_meta_data(object_all, column_to_transform = 'seurat_clusters', column_to_reference = 'filtered_condition', column_to_create = 'filtered_condition_imputed')


# save the result
object_all_cluster_filtered_ctd_cond_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_filtered_ctd_cond_20231129.rds', sep = '')
saveRDS(object_all, object_all_cluster_filtered_ctd_cond_loc)

# now repeat for just UT
object_ut <- object_all[, object_all@meta.data$soup_best_match_correlation > 0.7 & !is.na(object_all@meta.data[['condition']]) & object_all@meta.data[['condition']] == 'UT']
object_ut <- ScaleData(object_ut)
object_ut <- FindVariableFeatures(object_ut, layer = 'data')
object_ut <- RunPCA(object_ut)
object_ut <- RunUMAP(object_ut, dims = 1:30, return.model = T)
object_ut <- FindNeighbors(object_ut, dims = 1:30)
object_ut <- FindClusters(object_ut, resolution = 1.5)
# save result
object_ut_cluster_filtered_ctd_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_filtered_20231121.rds', sep = '')
saveRDS(object_ut, object_ut_cluster_filtered_ctd_loc)
# make plots
p_dim_ut <- DimPlot(object_ut)
p_markers_ut <- plot_grid(FeaturePlot(object_ut, features=c('CD14', 'CD19', 'CD3G', 'CD3D')), FeaturePlot(object_ut, features=c('CD4', 'CD74', 'CD8A', 'CST7')), FeaturePlot(object_ut, features=c('CTSS', 'NCAM1', 'FCGR3A', 'NKG7')), nrow = 1, ncol=3)
ggsave('~/mo_ut_souped_clus_filtered_ctd_cond_20231117_clusters.pdf', plot = p_dim_ut, width = 10, height = 10)
ggsave('~/mo_ut_marker_genes.pdf', plot = p_markers_ut, width = 20, height = 10)


# now repeat for just ca
object_ca <- object_all[, object_all@meta.data$soup_best_match_correlation > 0.7 & object_all@meta.data[['condition']] == '24hCa']
object_ca <- ScaleData(object_ca)
object_ca <- FindVariableFeatures(object_ca, layer = 'data')
object_ca <- RunPCA(object_ca)
object_ca <- RunUMAP(object_ca, dims = 1:30, return.model = T)
object_ca <- FindNeighbors(object_ca, dims = 1:30)
object_ca <- FindClusters(object_ca, resolcaion = 1.5)
plot_grid(FeaturePlot(object_ca, features=c('CD14', 'CD19', 'CD3G', 'CD3D')), FeaturePlot(object_ca, features=c('CD4', 'CD74', 'CD8A', 'CST7')), FeaturePlot(object_ca, features=c('CTSS', 'NCAM1', 'FCGR3A', 'NKG7')), nrow = 1, ncol=3)
# save result
object_ca_cluster_filtered_ctd_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_filtered_20231129.rds', sep = '')
saveRDS(object_ca, object_ca_cluster_filtered_ctd_loc)
# make plots
p_dim_ca <- DimPlot(object_ca)
p_markers_ca <- plot_grid(FeaturePlot(object_ca, features=c('CD14', 'CD19', 'CD3G', 'CD3D')), FeaturePlot(object_ca, features=c('CD4', 'CD74', 'CD8A', 'CST7')), FeaturePlot(object_ca, features=c('CTSS', 'NCAM1', 'FCGR3A', 'NKG7')), nrow = 1, ncol=3)
ggsave('~/mo_24hca_souped_clus_filtered_ctd_cond_20231117_clusters.pdf', plot = p_dim_ca, width = 10, height = 10)
ggsave('~/mo_24hca_marker_genes.pdf', plot = p_markers_ca, width = 20, height = 10)

