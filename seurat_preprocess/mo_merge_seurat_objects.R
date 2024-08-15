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


get_souporcell_output <- function(souporcell_output_loc, lanes, souporcell_file_name='clusters.tsv'){
  # init table
  soupor_table <- NULL
  # check each lane
  for(lane in lanes){
    # paste the full path
    full_output_loc <- paste(souporcell_output_loc, lane, '/', souporcell_file_name, sep = '')
    # read the file
    try({
      soupor_output <- read.table(full_output_loc, sep = '\t', header = T, stringsAsFactors = F)
      # we'll add the lane
      soupor_output[['lane']] <- lane
      # order the way we like it, with the barcode and lane first
      soupor_output <- soupor_output[, c('barcode', 'lane', setdiff(colnames(soupor_output), c('barcode', 'lane')))]
      # add to the rest
      if(is.null(soupor_table)){
        soupor_table <- soupor_output
      }
      else{
        # check if there are columns missing in the new table we just read
        only_existing_columns <- setdiff(colnames(soupor_table), colnames(soupor_output))
        # add those to the table we just read
        if(length(only_existing_columns) > 0){
          soupor_output[, only_existing_columns] <- NA
        }
        # check if there are columns missing in the existing table
        only_new_columns <- setdiff(colnames(soupor_output), colnames(soupor_table))
        # add those to the table we already have
        if(length(only_new_columns) > 0){
          soupor_table[, only_new_columns] <- NA
        }
        # now we can safely combine these
        soupor_table <- rbind(soupor_table, soupor_output)
      }
    })
  }
  return(soupor_table)
}

add_soup_assignments <- function(seurat_object, soupor_output){
  # create a regex to get the last index of -
  last_dash_pos <- "\\-[^\\-]*$"
  # remove the '-1' from the barcode
  soupor_output[['barcode']] <- substr(soupor_output[['barcode']], 1, regexpr(last_dash_pos, soupor_output[['barcode']])-1)
  # add a combination of the barcode and the lane
  rownames(soupor_output) <- paste(soupor_output[['barcode']], soupor_output[['lane']], sep = '_')
  # now remove the barcode and lane columns
  soupor_output[['barcode']] <- NULL
  soupor_output[['lane']] <- NULL
  # add the 'soup' prepend everywhere
  colnames(soupor_output) <- paste('soup', colnames(soupor_output), sep = '_')
  # now add each column to the object
  for(column in colnames(soupor_output)){
    seurat_object <- AddMetaData(seurat_object, soupor_output[column])
  }
  return(seurat_object)
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
           '230316_lane5', '230316_lane6', '230316_lane7', '230316_lane8'
)

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'

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
object_all_ <- mad_function(seurat = object_all, column = "percent.mt", number_mad = 3)
object_all <- mad_function(seurat = object_all, column = "nCount_RNA", number_mad = 3)
object_all <- mad_function(seurat = object_all, column = "nFeature_RNA", number_mad = 3)
# update to Seurat v5
object_loc_v5 <- paste(seurat_objects_loc, 'mo_all_20240204_seuratv5.rds', sep = '')
# update
object_all <- JoinLayers(object_all)
DefaultAssay(object_all) <- 'RNA'
saveRDS(object_all, object_loc_v5)


# get the soup per lane
soup_per_lane_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/gex/barcode_filtered/'
# initialize the object
object_all_azi <- NULL
# check each lane
for (lane in lanes) {
  # print progress
  print(paste('processing lane: ', lane, sep = ''))
  # read the object
  try({
    object_loc <- paste(seurat_objects_loc, '/', 'mo_', lane, '_multimodal_azi_mapped.rds', sep = '')
    # read object
    object_lane <- readRDS(object_loc)
    # update the path
    object_lane$peaks@fragments[[1]]@path <- gsub('tmp02', 'tmp03', object_lane$peaks@fragments[[1]]@path)
    # get the full souporcell output
    full_soupor_output <- get_souporcell_output(soup_per_lane_loc, c(lane))
    # add the souporcell assignments
    object_lane <- add_soup_assignments(object_lane, full_soupor_output)
    # remove doublets
    object_lane <- object_lane[, !is.na(object_lane@meta.data[['soup_status']]) & object_lane@meta.data[['soup_status']] == 'singlet']
    # merge
    if (!is.null(object_all_azi)) {
      object_all_azi <- merge(object_all_azi, object_lane)
    }
    else{
      object_all_azi <- object_lane
    }
  })
}
DefaultAssay(object_all_azi) <- 'RNA'
# add the percentage of MT
object_all_azi[["percent.mt"]] <- PercentageFeatureSet(object_all_azi, pattern = "^MT-")
# add MADs
object_all_azi <- mad_function(seurat = object_all_azi, column = "percent.mt", number_mad = 3)
object_all_azi <- mad_function(seurat = object_all_azi, column = "nCount_RNA", number_mad = 3)
object_all_azi <- mad_function(seurat = object_all_azi, column = "nFeature_RNA", number_mad = 3)
# update to Seurat v5
object_loc_v5 <- paste(seurat_objects_loc, 'mo_all_azi_20240204_seuratv5.rds', sep = '')
# update
object_all_azi <- JoinLayers(object_all_azi)
DefaultAssay(object_all_azi) <- 'RNA'
saveRDS(object_all_azi, object_loc_v5)
