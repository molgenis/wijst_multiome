#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_rna_qced_barcodes.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# required to create object
library(Seurat)

####################
# Functions        #
####################

#' get a vector of as distinct possible colours
#' 
#' @param seurat_metadata the seurat metadata
#' @param barcode_output_loc where to place the barcode files
#' @param lane_column which column in the metadata denotes the lane
#' @param barcodes_column which column in the metadata denotes the barcode to write
#' @returns 0 is succesful
#' 
write_barcodes_per_lane <- function(seurat_metadata, barcode_output_loc, lane_column='lane', barcodes_column='barcode_1') {
  # check each lane
  for (lane in unique(seurat_metadata[[lane_column]])) {
    # get the barcodes
    barcodes <- seurat_metadata[!is.na(seurat_metadata[[lane_column]]) & seurat_metadata[[lane_column]] == lane, barcodes_column]
    # set an output directory
    barcode_output_file <- paste(barcode_output_loc, '/', lane, '.tsv', sep = '')
    # and write the file
    write.table(data.frame(x = barcodes), barcode_output_file, quote = F, row.names = F, col.names = F)
  }
  return(0)
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
object_all_cluster_unfiltered_loc <- paste(seurat_objects_loc, 'mo_all_souped_clus_unfiltered_20231109.rds', sep = '')

# read the object
object_all <- readRDS(object_all_cluster_unfiltered_loc)

# subset by the nCount MAD
object_all <- object_all[, !is.na(object_all@meta.data[['nCount_RNA_mad']]) & object_all@meta.data[['nCount_RNA_mad']] == "NotOutlier"]
# subset by nFeature
object_all <- object_all[, !is.na(object_all@meta.data[['nFeature_RNA_mad']]) & object_all@meta.data[['nFeature_RNA_mad']] == "NotOutlier"]
# and with minimal expression
object_all <- object_all[, object_all@meta.data[['nCount_RNA']] >= 200 & object_all@meta.data[['nFeature_RNA']] >= 3]

# get the barcodes per lane
write_barcodes_per_lane(object_all@meta.data, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/demultiplexing/valid_barcodes/')
