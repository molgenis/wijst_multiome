#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_azimuth_reference_map.R
# Function: do cell type annotation 
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)


####################
# Functions        #
####################

#' do Azimuth reference mapping of query onto a reference
#' 
#' @param reference the reference Seurat object
#' @param query the query Seurat object
#' @returns a Seurat object
#' 
do_reference_mapping <- function(reference, query){
  # find transfer anchors
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "SCT",
    reference.reduction = "pca",
    dims = 1:50
  )
  # do the reference mapping
  query <- MapQuery(
    anchorset = anchors,
    query = query,
    reference = reference,
    refdata = list(
      cell_type_1M_2022 = "cell_type"
    ),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
  return(query)
}


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')


####################
# Main Code        #
####################

# these are the lanes
lanes <- c('230105_lane1', '230105_lane2', '230105_lane3', '230105_lane4',
           '230105_lane5', '230105_lane6', '230105_lane7', '230105_lane8',
           '230112_lane1', '230112_lane2', '230112_lane3', '230112_lane4',
           '230112_lane5', '230112_lane6', '230112_lane7', '230112_lane8',
           '230120_lane1', '230120_lane2', '230120_lane3', '230120_lane4',
           '230120_lane5', '230120_lane6', '230120_lane7', '230120_lane8',
           '230127_lane1', '230127_lane2', '230127_lane3', '230127_lane4',
           '230127_lane5', '230127_lane6', '230127_lane7', '230127_lane8',
           '230202_lane1', '230202_lane2', '230202_lane3', '230202_lane4',
           '230202_lane5', '230202_lane6', '230202_lane7', '230202_lane8',
           '230209_lane1', '230209_lane2', '230209_lane3', '230209_lane4',
           '230209_lane5', '230209_lane6', '230209_lane7', '230209_lane8',
           '230216_lane1', '230216_lane2', '230216_lane3', '230216_lane4',
           '230216_lane5', '230216_lane6', '230216_lane7', '230216_lane8',
           '230223_lane1', '230223_lane2', '230223_lane3', '230223_lane4',
           '230223_lane5', '230223_lane6', '230223_lane7', '230223_lane8',
           '230302_lane1', '230302_lane2', '230302_lane3', '230302_lane4',
           '230302_lane5', '230302_lane6', '230302_lane7', '230302_lane8',
           '230316_lane1', '230316_lane2', '230316_lane3', '230316_lane4',
           '230316_lane5', '230316_lane6', '230316_lane7'
)

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects'

# location of the reference
reference_loc <- '/groups/umcg-franke-scrna/tmp02/releases/wijst-2020-hg38/v1/seurat/1M_v3_mediumQC_sct_b38.rds'
# load reference
reference <- readRDS(reference_loc)
# pca
reference <- RunPCA(reference)
# set seed
set.seed(1337)
# umap
reference <- RunUMAP(reference, dims = 1:30, return.model = T)
# knn
reference <- FindNeighbors(reference, dims = 1:30)
# find clusters
reference <- FindClusters(reference, resolution = 1.2)
# save result
saveRDS(reference, '/groups/umcg-franke-scrna/tmp02/releases/wijst-2020-hg38/v1/seurat/1M_v3_mediumQC_sct_clus_b38.rds')

# save the cell type annotation
cell_type_predictions_list <- list()

# check each lane
for (lane in lanes) {
  # print progress
  print(paste('processing lane: ', lane, sep = ''))
  # read the object
  try({
    object_loc <- paste(seurat_objects_loc, '/', 'mo_', lane, '.rds', sep = '')
    # read object
    object_lane <- readRDS(object_loc)
    # normalize
    object_lane <- SCTransform(object_lane, vst.flavor = 'v2')
    # pca
    object_lane <- RunPCA(object_lane)
    # umap
    object_lane <- RunUMAP(object_lane, dims = 1:30, return.model = T)
    # knn
    object_lane <- FindNeighbors(object_lane, dims = 1:30)
    # find clusters
    object_lane <- FindClusters(object_lane, resolution = 1.2)
    # do reference mapping
    object_lane <- do_reference_mapping(reference, object_lane)
    # we'll save the cell type
    cell_type_prediction_lane <- data.frame(barcode = rownames(object_lane@meta.data), cell_type = object_lane@meta.data[['predicted.cell_type_1M_2022']], score = object_lane@meta.data[['predicted.cell_type_1M_2022.score']])
    cell_type_predictions_list[[lane]] <- cell_type_prediction_lane
    # save result
    result_loc <- paste(paste(seurat_objects_loc, '/', 'mo_', lane, '_azi_1m_v3.rds', sep = ''))
    saveRDS(object_lane, result_loc)
  })
}

# merge all predictions
cell_type_predictions <- do.call('rbind', cell_type_predictions_list)
# write the result
write.table(cell_type_predictions, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cell_type_assignment/azimuth/NC2022_v3/mo_azimuth_ct_nc2022_v2.tsv', sep = '\t', row.names = F, col.names = T, quote = T)
