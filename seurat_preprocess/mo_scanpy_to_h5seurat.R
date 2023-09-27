#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_scanpy_to_h5seurat.R
# Function: preprocess the count data
############################################################################################################################


####################
# libraries        #
####################

# required to create object
library(Seurat)
library(SeuratDisk)

####################
# Functions        #
####################

convert_scanpy_to_h5seurat <- function(scanpy_objects_loc, h5seurat_objects_loc, lanes, scanpy_object_prepend='mo_', scanpy_object_append='.h5ad', h5seurat_object_prepend='mo_', h5seurat_object_append='.h5seurat') {
  # check each lane
  for (lane in lanes) {
    # get the path to the scanpy file
    scanpy_object_loc <- paste(scanpy_objects_loc, scanpy_object_prepend, lane, scanpy_object_append, sep = '')
    # and the seurat file
    h5seurat_object_loc <- paste(h5seurat_objects_loc, h5seurat_object_prepend, lane, h5seurat_object_append, sep = '')
    # check if the original exists
    if (file.exists(scanpy_object_loc)) {
      # actually convert
      try({
        # an error will be thrown when closing, but the file created is still valid: https://github.com/mojaveazure/seurat-disk/issues/82
        Convert(scanpy_object_loc, dest = h5seurat_object_loc, overwrite = TRUE)
      })
    }
    else {
      warning(paste('missing for lane ', lane, ', at ', scanpy_object_loc, '. Skipping', sep = ''))
    }
  }
  return(0)
}
 

####################
# Main Code        #
####################

# the location of the scanpy objects
scanpy_objects_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scanpy_preprocess_samples/objects/'
# the location where we want the h5Seurat objects
h5seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/'

# these are the lanes we care about
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

# now let's try to convert
convert_scanpy_to_h5seurat(scanpy_objects_loc, h5seurat_objects_loc, lanes)
