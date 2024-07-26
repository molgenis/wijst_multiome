#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_azimuth_reference_map_pbmc.R
# Function: reference mapping of multiome data using Azimuth reference PBMC dataset
############################################################################################################################

####################
# libraries        #
####################

library(SeuratData)
library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(cowplot)

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
    reference.reduction = "spca",
    dims = 1:50,
    recompute.residuals = FALSE
  )

  # do the reference mapping
  query <- MapQuery(
    anchorset = anchors,
    query = query,
    reference = reference,
    refdata = list(
      predicted.azipbmc.l1 = "celltype.l1",
      predicted.azipbmc.l2 = "celltype.l2",
      predicted.azipbmc.ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  return(query)
}


process_lane <- function(seurat_objects_full, lane, lane_column='lane') {
  # read the object
  tryCatch({
    #object_loc <- paste(seurat_objects_dir, '/', seurat_prepend, lane, seurat_append, sep = '')
    # read RNA object
    #cbmc <- readRDS(object_loc)
    cbmc <- seurat_objects_full[, !is.na(seurat_objects_full@meta.data[[lane_column]]) & seurat_objects_full@meta.data[[lane_column]] == lane]
    # recreate the object
    cbmc <- CreateSeuratObject(counts = CreateAssayObject(cbmc@assays$RNA$counts), meta.data = cbmc@meta.data)
    # filter by UMI
    cbmc <- cbmc[, cbmc@meta.data[['nCount_RNA']] >= 200]
    cbmc <- cbmc[, cbmc@meta.data[['nFeature_RNA']] >= 200]
    # first do RNA
    DefaultAssay(cbmc) <- "RNA"
    # perform visualization and clustering steps
    cbmc <- NormalizeData(cbmc)
    cbmc <- FindVariableFeatures(cbmc)
    cbmc <- ScaleData(cbmc)
    cbmc <- RunPCA(cbmc, verbose = FALSE)
    cbmc <- FindNeighbors(cbmc, dims = 1:30, reduction = 'pca')
    cbmc <- FindClusters(cbmc, resolution = 1.2, verbose = FALSE)
    cbmc <- RunUMAP(cbmc, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", return.model = T)
    
    # now do SCT
    cbmc <- SCTransform(cbmc)
    cbmc <- RunPCA(cbmc, verbose = FALSE)
    cbmc <- FindNeighbors(cbmc, dims = 1:30, reduction = 'pca')
    cbmc <- FindClusters(cbmc, resolution = 1.2, verbose = FALSE)
    cbmc <- RunUMAP(cbmc, dims = 1:30, reduction.name = "umap.sct", reduction.key = "sctUMAP_", return.model = T)

    return(cbmc)
  }, error = function(err) {
    print(paste('error in processing', err))
    return(NULL)
  })
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


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
set.seed(7777)

####################
# Main Code        #
####################

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
           '230316_lane5', '230316_lane6', '230316_lane7', '230316_lane8'
)

# location of per-lane objects
seurat_objects_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
# the full object loc
full_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240204_seuratv5.rds'
# read the object
full_object <- readRDS(full_object_loc)

# location of the reference
reference_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/azimuth_pbmc_reference/pbmc_multimodal.h5seurat'
reference_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/azimuth_pbmc_reference/multi.h5seurat'
# read the reference
reference <- LoadH5Seurat(reference_loc)
reference <- UpdateSeuratObject(reference)

# do each lane
for (lane in lanes) {
  print(paste('starting on reading lane', lane))
  processed_object <- process_lane(seurat_objects_full = full_object, lane = lane)
  # try to do mapping
  if (!is.null(processed_object)) {
    print(paste('starting reference map of lane', lane))
    processed_object <- do_reference_mapping(reference = reference, query = processed_object)
    # save result
    result_loc <- paste(seurat_objects_dir, 'mo_', lane, '_pbmcref_azi_mapped.rds', sep = '')
    print(paste('saving object for lane', lane))
    saveRDS(processed_object, result_loc)
  }
}

# now summarize all
all_ct_predictions_per_lane <- list()
for (lane in lanes) {
  # load the object
  result_loc <- paste(seurat_objects_dir, 'mo_', lane, '_pbmcref_azi_mapped.rds', sep = '')
  try({
    processed_object <- readRDS(result_loc)
    # extract the data we want
    annotation_lane <- data.frame(barcode = processed_object@meta.data[['barcode_lane']],
                                  predicted.azipbmc.l1 = processed_object@meta.data[['predicted.predicted.azipbmc.l1']],
                                  predicted.azipbmc.l1.score = processed_object@meta.data[['predicted.predicted.azipbmc.l1.score']],
                                  predicted.azipbmc.l2 = processed_object@meta.data[['predicted.predicted.azipbmc.l2']],
                                  predicted.azipbmc.l2.score = processed_object@meta.data[['predicted.predicted.azipbmc.l2.score']]
                                  #predicted.azipbmc.ADT = processed_object@meta.data[['predicted.azipbmc.ADT']],
                                  #predicted.azipbmc.ADT.score = processed_object@meta.data[['predicted.azipbmc.ADT.score']]
                                  )
    # add to the list
    all_ct_predictions_per_lane[[lane]] <- annotation_lane
  })
}
# merge all
all_ct_predictions <- do.call('rbind', all_ct_predictions_per_lane)
# write the results
write.table(all_ct_predictions, '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cell_type_assignment/azimuth/azimuth_PBMCs/mo_azimuth_ct_azipbmc.tsv', sep = '\t', row.names = F, col.names = T)
