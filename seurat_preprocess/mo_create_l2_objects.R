#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_l2_objects.R
# Function: create L2 azimuth specific Seurat objects
############################################################################################################################

####################
# libraries        #
####################

# read the object
library(Seurat)

####################
# Functions        #
####################

make_celltypes_safe <- function(cell_types){
  # make the NA unannotated
  cell_types[is.na(cell_types)] <- 'unannotated'
  # get a safe file name
  cell_type_safes <- gsub(' |/', '_', cell_types)
  cell_type_safes <- gsub('-', '_negative', cell_type_safes)
  cell_type_safes <- gsub('\\+', '_positive', cell_type_safes)
  cell_type_safes <- gsub('\\)', '', cell_type_safes)
  cell_type_safes <- gsub('\\(', '', cell_type_safes)
  return(cell_type_safes)
}


####################
# Main Code        #
####################

# location of the Seurat objects
seurat_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds'
# load object
mo <- readRDS(seurat_object_loc)

# location of annotation for cell types
azi_pbmc_annotations_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cell_type_assignment/azimuth/azimuth_PBMCs/mo_azimuth_ct_azipbmc.tsv'
azi_pbmc_annotations <- read.table(azi_pbmc_annotations_loc, header =T, sep = '\t', row.names = 1)

# add these annotations
mo <- AddMetaData(mo, azi_pbmc_annotations[, colnames(azi_pbmc_annotations)])

# make the l2 cell types without illegal characters
mo@meta.data[['predicted.azipbmc.l2.safe']] <- make_celltypes_safe(mo@meta.data[['predicted.azipbmc.l2']])

# now go through each cell type
for (cell_type in unique(mo@meta.data$predicted.azipbmc.l2.safe)) {
  # check for NA
  if (!is.na(cell_type)) {
    # subset to this celltype and condition
    mo_covid_ct <- mo[, !is.na(mo@meta.data$predicted.azipbmc.l2.safe) &
                        !is.na(mo@meta.data$condition_final) &
                        mo@meta.data$predicted.azipbmc.l2.safe == cell_type &
                        mo@meta.data$condition_final == 'UT']
    # write this file
    saveRDS(mo_covid_ct,
            paste('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240619_seuratv5_annotated_agesexcovid_', cell_type, '_UT.rds', sep = ''))
    
  }
}

# now go through each cell type without taking UT only
for (cell_type in unique(mo@meta.data$predicted.azipbmc.l2.safe)) {
  # check for NA
  if (!is.na(cell_type)) {
    # subset to this celltype and condition
    mo_covid_ct <- mo[, !is.na(mo@meta.data$predicted.azipbmc.l2.safe) &
                        !is.na(mo@meta.data$condition_final) &
                        mo@meta.data$predicted.azipbmc.l2.safe == cell_type]
    # write this file
    saveRDS(mo_covid_ct,
            paste('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240619_seuratv5_annotated_agesexcovid_', cell_type, '.rds', sep = ''))
    
  }
}
