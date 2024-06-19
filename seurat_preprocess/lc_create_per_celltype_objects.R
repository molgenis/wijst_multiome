#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lc_create_per_celltype_objects.R
# Function: 
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)


####################
# Functions        #
####################

#' add the inflammation assignments  to the Seurat object
#' 
#' @param seurat_object The Seurat object to add the inflammation status to
#' @param sample_sheet The sample sheet containing lanes, participants and inflammation statuses
#' @param seurat_lane_column The column in the Seurat metadata denoting the 10x lane
#' @param sheet_lane_column The column in the sample sheet denoting the 10x lane
#' @param seurat_participant_column The column in the Seurat metadata denoting the participant assignment
#' @param sheet_participants_column The column in the sample sheet containing the participants per lane
#' @param seurat_inflammation_column The column in the Seurat metadata to add the inflammation status in
#' @param sheet_inflammation_column The column in the sample sheet containing the inflammation statuses per lane
#' @returns the Seurat object with the inflammation status added
#' lpmcv2 <- add_inflammation_status(lpmcv2, sample_sheet)
add_inflammation_status <- function(seurat_object, sample_sheet, seurat_lane_column='lane', sheet_lane_column='lane', seurat_participant_column='soup_final_sample_assignment', sheet_participants_column='genoid', seurat_inflammation_column='inflammation_status', sheet_inflammation_column='inflammation_status') {
  # create a mapping of lane+sample to inflammation status
  mapping_per_lane_list <- list()
  for (i in 1:nrow(sample_sheet)) {
    # extract lane
    lane <- sample_sheet[i, sheet_lane_column]
    # extract the participants
    participant <- sample_sheet[i, sheet_participants_column]
    # and the inflammation condition
    condition <- sample_sheet[i, sheet_inflammation_column]
    # if not set, set to unknown
    if (is.null(condition)) {
      condition <- 'unknown'
    }
    # subset to the barcodes which have are this lane and participant
    barcodes_match <- rownames(seurat_object@meta.data[!is.na(seurat_object@meta.data[[seurat_lane_column]]) &
                                                         seurat_object@meta.data[[seurat_lane_column]] == lane &
                                                         !is.na(seurat_object@meta.data[[seurat_participant_column]]) &
                                                         seurat_object@meta.data[[seurat_participant_column]] == participant, ])
    # only add if there are matching barcodes
    if (length(barcodes_match) > 0) {
      # create dataframe
      df_lane_part <- data.frame(barcode = barcodes_match, condition = rep(condition, times = length(barcodes_match)))
      # set the colname to be the one we chose
      colnames(df_lane_part) <- c('barcode', seurat_inflammation_column)
      # then add to the list
      mapping_per_lane_list[[paste(lane, participant, sep = ':')]] <- df_lane_part
    }
  }
  # now merge all together
  mapping_all <- do.call('rbind', mapping_per_lane_list)
  # set the barcode as rownames
  rownames(mapping_all) <- mapping_all[['barcode']]
  # finally add to the object
  seurat_object <- AddMetaData(seurat_object, mapping_all[seurat_inflammation_column])
  return(seurat_object)
}



####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 120 * 2000 * 2024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# location of the objects
seurat_object_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated.rds'

# load object
mo <- readRDS(seurat_object_loc)

# location of the condition assignment
condition_assignment_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_monocyte_based_condition_numbers.tsv'
# read the conditions
condition_assignments <- read.table(condition_assignment_loc, header = T, sep = '\t')
# location of the age/sex assignments
age_sex_assignments_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_age_sex_batch12.tsv'
# read the assignments
age_sex_assigments <- read.table(age_sex_assignments_loc, header = T, sep = '\t')
# the ugli ones as well
age_sex_assigments_ugli_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_age_sex_ugli.tsv'
age_sex_assigments_ugli <- read.table(age_sex_assigments_ugli_loc, header = T, sep = '\t')
# merge
age_sex_assigments <- rbind(age_sex_assigments, age_sex_assigments_ugli)
# and the 'realids'
realid_assignments_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'
# read the realids
realid_assignments <- read.table(realid_assignments_loc, header = T, sep = '\t')

# do the condition assignment
mo <- add_inflammation_status(
  mo,
  sample_sheet=condition_assignments, 
  seurat_lane_column='lane',
  sheet_lane_column='lane', 
  seurat_participant_column='sample_final', 
  sheet_participants_column='sample_final', 
  seurat_inflammation_column='condition_sheet', 
  sheet_inflammation_column='condition'
)
mo <- add_inflammation_status(
  mo,
  sample_sheet=condition_assignments, 
  seurat_lane_column='lane',
  sheet_lane_column='lane', 
  seurat_participant_column='sample_final', 
  sheet_participants_column='sample', 
  seurat_inflammation_column='condition_prev', 
  sheet_inflammation_column='cond_prev'
)
# rename 24hCA
mo@meta.data[!is.na(mo@meta.data[['condition_sheet']]) &
                          mo@meta.data[['condition_sheet']] == '24hCa', 'condition_sheet'] <- '24hCA'
# now also set the final inflammation assignment
mo@meta.data[['condition_final']] <- mo@meta.data[['condition_sheet']]
mo@meta.data[is.na(mo@meta.data[['condition_final']]), 'condition_final'] <- mo@meta.data[is.na(mo@meta.data[['condition_final']]), 'condition_prev']

# remove any empty entries
mo <- mo[, !is.na(mo@meta.data[['celltype_imputed_lowerres']]) &
                                 !is.na(mo@meta.data[['sample_final']]) &
                                 !is.na(mo@meta.data[['condition_final']])]

# fix the sample IDs in the realid table
realid_assignments[['sample_final']] <- realid_assignments$SampleID_CORRECT
realid_assignments[is.na(realid_assignments$SampleID_CORRECT) | realid_assignments$SampleID_CORRECT == '', 'sample_final']  <- realid_assignments[is.na(realid_assignments$SampleID_CORRECT) | realid_assignments$SampleID_CORRECT == '', 'sample']
realid_assignments[['realid_final']] <- realid_assignments$RealID
realid_assignments[is.na(realid_assignments$RealID) | realid_assignments$RealID == '', 'realid_final']  <- realid_assignments[is.na(realid_assignments$RealID) | realid_assignments$RealID == '', 'sample_final']
# add the realid
mo@meta.data['realid'] <- realid_assignments[match(mo@meta.data$sample_final, realid_assignments$sample_final), 'realid_final']

# check which samples are not in the age/sex file
not_age_sex <- data.frame(mo = setdiff(unique(mo@meta.data$realid), age_sex_assigments$sample))
write.table(not_age_sex, '~/mo_age_sex_missing.tsv', row.names = F, col.names = T, quote = F, sep = '\t')
# add the age and sex
mo@meta.data[['age']] <- age_sex_assigments[match(mo@meta.data$realid, age_sex_assigments$sample), 'age']
mo@meta.data[['sex']] <- age_sex_assigments[match(mo@meta.data$realid, age_sex_assigments$sample), 'sex']
# add the final LONG_COVID assignment
mo@meta.data[['LONG_COVID_final']] <- mo@meta.data$LONG_COVID
mo@meta.data[is.na(mo@meta.data[['LONG_COVID']]), 'LONG_COVID_final'] <- 'control'
mo@meta.data[['LONG_COVID_method']] <- NA
mo@meta.data[is.na(mo@meta.data[['LONG_COVID']]), 'LONG_COVID_method'] <- 'inferred'
mo@meta.data[!is.na(mo@meta.data[['LONG_COVID']]), 'LONG_COVID_method'] <- 'assigned'
# save the result
saveRDS(mo, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds')
# create a new sample mapping file
full_sample_mapping <- unique(mo@meta.data[, c('lane', 'sample_final', 'realid', 'condition_final', 'age', 'sex', 'LONG_COVID_final', 'LONG_COVID_method')])
write.table(full_sample_mapping, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_full_sample_mapping.tsv', row.names = F, col.names = T, quote = F, sep = '\t')
# now go through each cell type
for (cell_type in unique(mo@meta.data$celltype_imputed_lowerres)) {
  # check for NA
  if (!is.na(cell_type)) {
    # subset to this celltype and condition
    mo_covid_ct <- mo[, !is.na(mo@meta.data$celltype_imputed_lowerres) &
                        !is.na(mo@meta.data$condition_final) &
                        mo@meta.data$celltype_imputed_lowerres == cell_type &
                        mo@meta.data$condition_final == 'UT']
    # write this file
    saveRDS(mo_covid_ct,
            paste('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240619_seuratv5_annotated_agesexcovid_', cell_type, '_UT.rds', sep = ''))
    
  }
}
