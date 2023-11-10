#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Drew Neavin
# Name: mo_merge_seurat_objects.R
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

#' add the MADs
#' 
#' @param seurat the Seurat object
#' @param column the column to calculate MAD for
#' @param number_mad how many mads to be considered an outlier
#' @returns a Seurat object
#' 
mad_function <- function(seurat, column, number_mad=3){
  mad <- mad(seurat@meta.data[,column])
  low <- median(seurat@meta.data[,column]) - number_mad*mad
  high <- median(seurat@meta.data[,column]) + number_mad*mad
  print("The lower bound is:")
  print(low)
  print("The upper bound is:")
  print(high)
  seurat@meta.data[,paste0(column,"_mad")] <- ifelse((seurat@meta.data[,column] > low & seurat@meta.data[,column] < high),"NotOutlier", "Outlier")
  return(seurat)
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
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/seurat_preprocess_samples/objects/'

# initialize the object
object_all <- NULL

# check each lane
for (lane in lanes) {
  # print progress
  print(paste('processing lane: ', lane, sep = ''))
  # read the object
  try({
    object_loc <- paste(seurat_objects_loc, '/', 'mo_', lane, '.rds', sep = '')
    # read object
    object_lane <- readRDS(object_loc)
    # merge
    if (!is.null(object_all)) {
      object_all <- merge(object_all, object_lane)
    }
    else{
      object_all <- object_lane
    }
  })
}

# add the percentage of MT
object_all[["percent.mt"]] <- PercentageFeatureSet(object_all, pattern = "^MT-")
# add MADs
object_all <- mad_function(seurat = object_all, column = "percent.mt", number_mad = 3)
object_all <- mad_function(seurat = object_all, column = "nCount_RNA", number_mad = 3)
object_all <- mad_function(seurat = object_all, column = "nFeature_RNA", number_mad = 3)
# get the souporcell output
soup_out_all_rna_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_sample_matched.tsv'
soup_out_all_rna <- read.table(soup_out_all_rna_loc, header = T, sep = '\t', row.names = 3)
# remove columns we don't need
soup_out_all_rna <- soup_out_all_rna[, setdiff(colnames(soup_out_all_rna), c('lane', 'barcode', 'barcode_original'))]
# modify column names
colnames(soup_out_all_rna) <- paste('soup', colnames(soup_out_all_rna), sep = '_')
# add to object
object_all <- AddMetaData(object_all, soup_out_all_rna)
# save the object somewhere
object_loc_souped <- paste(seurat_objects_loc, 'mo_all_souped_20231109.rds', sep = '')
saveRDS(object_all, object_loc_souped)
