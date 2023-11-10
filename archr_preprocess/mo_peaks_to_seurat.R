#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_peaks_to_seurat.R
# Function: peaks to Seurat object
############################################################################################################################


####################
# libraries        #
####################

# required to create object
library(Seurat)


####################
# Functions        #
####################

#' get a Seurat object from a combination of matrix, barcodes and features files
#' 
#' @param lane the lane to load
#' @param base_counts_dir base location of the input per lane output
#' @param matrix_dir possible subdirectory or prepend for each input file
#' @param min.cells minimal cells for a gene to be included
#' @param min.features minimal unique transcript for a cell to be included
#' @returns a Seurat object
#' 
add_data <- function(lane, base_counts_dir, matrix_dir='', min.cells = 0, min.features = 0) {
  # get the count directory
  counts_dir <- paste0(base_counts_dir, '/', lane, '/', matrix_dir, sep = '')
  
  # read the actual counts
  counts <- Read10X(counts_dir)[['Peaks']]
  
  # create some metadata, for now, we'll first just store the lane here
  barcodes_short <- gsub('(-\\d+)', '', colnames(counts))
  barcodes_lane <- paste(barcodes_short, rep(lane, times = length(barcodes_short)), sep = '_')
  
  metadata <- data.frame(barcode=barcodes_short,
                         barcode_1=colnames(counts),
                         barcode_lane=barcodes_lane)
  # save the lane the samples are from
  metadata$lane <- lane
  metadata$batch <- lane
  rownames(metadata) <- metadata$barcode_lane
  # change to a barcode unique across lanes
  colnames(counts) <- metadata$barcode_lane
  
  # create the actual object
  seurat_new <- Seurat::CreateSeuratObject(counts = counts,
                                           min.cells = min.cells,
                                           min.features = min.features,
                                           project = "wijst_multiome",
                                           meta.data = metadata)
  return(seurat_new)
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

# location of the deconstructed matrices
matrices_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/'

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/archr_preprocess_samples/objects'

# check each lane
for (lane in lanes) {
  # print progress
  print(paste('processing lane: ', lane, sep = ''))
  
  try({
    # create the object
    object_lane <- add_data(lane, matrices_loc, matrix_dir = 'outs/filtered_feature_bc_matrix/')
    # write the result
    object_loc <- paste(seurat_objects_loc, '/', 'mo_', lane, '_peaks.rds', sep = '')
    saveRDS(object_lane, object_loc)
  })
}