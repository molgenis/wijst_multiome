#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_add_rna_metadata.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(Seurat)
library(matrixStats)
library(textTinyR) # NOT IN CONTAINER
library(pbapply)

####################
# Functions        #
####################




####################
# Options          #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory, we'll take 95 percent of the max memory in GB
mem_gb_requested <- 1900
options(future.globals.maxSize = (0.95 * mem_gb_requested) * 1000 * 1024^2)

# set the max parallel threads
#register(MulticoreParam(4))

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# location of the condition assignment
condition_assignment_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_monocyte_based_condition_numbers.tsv'
# read the conditions
condition_assignments <- read.table(condition_assignment_loc, header = T, sep = '\t')
# location of the LONG-CoVID assignments
longcovid_assignments_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_longcovid_assignments.tsv'
# read the assignments
longcovid_assignments <- read.table(longcovid_assignments_loc, header = T, sep = '\t')


# locations of objects
objects_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
mo_object_loc <- paste(objects_loc, 'mo_all_20240223_seuratv5_normalized.rds', sep = '')
# read the Seurat object
seurat_object <- readRDS(mo_object_loc)
# get the lane remapping
lane_remapping <- combine_lanes(unique(seurat_object@meta.data$lane))
# add to the object
seurat_object@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[seurat_object@meta.data[['lane']]]))
# do the condition assignment
seurat_object <- add_inflammation_status(
  seurat_object,
  sample_sheet=condition_assignments, 
  seurat_lane_column='lane',
  sheet_lane_column='lane', 
  seurat_participant_column='sample_final', 
  sheet_participants_column='sample_final', 
  seurat_inflammation_column='condition_sheet', 
  sheet_inflammation_column='condition'
)
seurat_object <- add_inflammation_status(
  seurat_object,
  sample_sheet=condition_assignments, 
  seurat_lane_column='lane',
  sheet_lane_column='lane', 
  seurat_participant_column='sample_final', 
  sheet_participants_column='sample', 
  seurat_inflammation_column='condition_prev', 
  sheet_inflammation_column='cond_prev'
)
# rename 24hCA
seurat_object@meta.data[!is.na(seurat_object@meta.data[['condition_sheet']]) &
                          seurat_object@meta.data[['condition_sheet']] == '24hCa', 'condition_sheet'] <- '24hCA'
# now also set the final inflammation assignment
seurat_object@meta.data[['condition_final']] <- seurat_object@meta.data[['condition_sheet']]
seurat_object@meta.data[is.na(seurat_object@meta.data[['condition_final']]), 'condition_final'] <- seurat_object@meta.data[is.na(seurat_object@meta.data[['condition_final']]), 'condition_prev']

# add the first part to the name of the long covid ID
longcovid_assignments[['mo']] <- paste('MO', longcovid_assignments$Project_ID, sep = '')
# add the longcovid assignment first
seurat_object@meta.data[['LONG_COVID']] <- NA
# now add the ones we have
seurat_object@meta.data[seurat_object@meta.data[['sample_final']] %in% longcovid_assignments[['mo']], 'LONG_COVID'] <- longcovid_assignments[match(seurat_object@meta.data[seurat_object@meta.data[['sample_final']] %in% longcovid_assignments[['mo']], 'sample_final'], longcovid_assignments[['mo']]), 'Case.Control']
# set empty to NA
seurat_object@meta.data[!is.na(seurat_object@meta.data[['LONG_COVID']]) & seurat_object@meta.data[['LONG_COVID']] == '', 'LONG_COVID'] <- NA
# save result
saveRDS(seurat_object, paste(objects_loc, 'mo_all_20240517_seuratv5_annotated.rds', sep = ''))

# save the metadata as well
write.table(seurat_object@meta.data, gzfile('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'), sep = '\t', row.names = F, col.names = T)