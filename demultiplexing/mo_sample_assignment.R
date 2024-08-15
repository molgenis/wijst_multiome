#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_sample_assignment.R
# Function: assign samples from demultiplexing
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
    # subset to the barcodes which have are this lane and participant
    barcodes_match <- rownames(seurat_object@meta.data[!is.na(seurat_object@meta.data[[seurat_lane_column]]) &
                                                         seurat_object@meta.data[[seurat_lane_column]] == lane &
                                                         !is.na(seurat_object@meta.data[[seurat_participant_column]]) &
                                                         seurat_object@meta.data[[seurat_participant_column]] == participant, ])
    # create dataframe
    df_lane_part <- data.frame(barcode = barcodes_match, condition = rep(condition, times = length(barcodes_match)))
    # set the colname to be the one we chose
    colnames(df_lane_part) <- c('barcode', seurat_inflammation_column)
    # then add to the list
    mapping_per_lane_list[[paste(lane, participant, sep = ':')]] <- df_lane_part
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

# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'

# load object
mo_loc <- paste(seurat_objects_loc, 'mo_all_20240204_seuratv5.rds', sep = '')
mo <- readRDS(mo_loc)

# get the assignment matrices
correlation_mapping_per_barcode <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched.tsv', header = T, sep = '\t')
correlation_mapping_per_barcode_all <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched_vs_all.tsv', header = T, sep = '\t')
# set barcodes and remove data we already have
rownames(correlation_mapping_per_barcode) <- correlation_mapping_per_barcode[['barcode_lane']]
correlation_mapping_per_barcode[, c('lane', 'barcode_lane', 'barcode', 'barcode_original')] <- NULL
rownames(correlation_mapping_per_barcode_all) <- correlation_mapping_per_barcode_all[['barcode_lane']]
correlation_mapping_per_barcode_all[, c('lane', 'barcode_lane', 'barcode', 'barcode_original')] <- NULL
# now let's get the souporcell data specifically, which would be the same for all and per-lane
soup_only <- correlation_mapping_per_barcode[, setdiff(colnames(correlation_mapping_per_barcode), c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation'))]
# and the correlation data
correlations_confined <- correlation_mapping_per_barcode[, c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')]
correlations_unconfined <- correlation_mapping_per_barcode_all[, c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')]
# rename the columns to denote they are soup specific
colnames(soup_only) <- paste('soup', colnames(soup_only), sep = '_')
# and specify the correlation confinement for the sample assignments as well
colnames(correlations_confined) <- paste('confined', colnames(correlations_confined), sep = '_')
colnames(correlations_unconfined) <- paste('unconfined', colnames(correlations_unconfined), sep = '_')

# get the inflammation assignments
inflammation_assignments_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'
inflammation_assignments <- read.table(inflammation_assignments_loc, header = T, sep = '\t', stringsAsFactors = F)

# fix the sample sheet
inflammation_assignments[['final_sample']] <- inflammation_assignments[['sample']]
inflammation_assignments[!is.na(inflammation_assignments[['SampleID_CORRECT']]) & inflammation_assignments[['SampleID_CORRECT']] != '', 'final_sample'] <- inflammation_assignments[!is.na(inflammation_assignments[['SampleID_CORRECT']]) & inflammation_assignments[['SampleID_CORRECT']] != '', 'SampleID_CORRECT']

# add the confined sample
mo <- AddMetaData(mo, correlations_confined[, colnames(correlations_confined)])
# and inflammation status
mo <- add_inflammation_status(mo, sample_sheet = inflammation_assignments, 
                        seurat_lane_column = 'lane', 
                        sheet_lane_column = 'lane', 
                        seurat_participant_column = 'confined_best_match_sample',
                        sheet_participants_column='final_sample', 
                        seurat_inflammation_column='confined_condition', 
                        sheet_inflammation_column='condition')
# and for unconfined
mo <- AddMetaData(mo, correlations_unconfined[, colnames(correlations_unconfined)])
mo <- add_inflammation_status(mo, sample_sheet = inflammation_assignments, 
                        seurat_lane_column = 'lane', 
                        sheet_lane_column = 'lane', 
                        seurat_participant_column = 'unconfined_best_match_sample',
                        sheet_participants_column='sample', 
                        seurat_inflammation_column='unconfined_condition', 
                        sheet_inflammation_column='condition')

# finally just the soup data
mo <- AddMetaData(mo, soup_only)

# set the final assignments
mo@meta.data[['sample_final']] <- NA
# set those with high enough correlations
mo@meta.data[!is.na(mo@meta.data[['unconfined_best_match_correlation']]) & mo@meta.data[['unconfined_best_match_correlation']] > 0.5, 'sample_final'] <- mo@meta.data[!is.na(mo@meta.data[['unconfined_best_match_correlation']]) & mo@meta.data[['unconfined_best_match_correlation']] > 0.5, 'unconfined_best_match_sample']
# this one has the best correlation for MO408 and MO1001. MO1001 should be present, so that one we will assign it to
mo@meta.data[mo@meta.data[['lane']] == '230112_lane1' & !is.na(mo@meta.data[['unconfined_best_match_correlation']]) & mo@meta.data[['unconfined_best_match_correlation']] == 'MO408', 'sample_final'] <- 'MO1001'
# this one has the best correlation for MO97 and MO203. MO203 should be present, so that one we will assign it to
mo@meta.data[mo@meta.data[['lane']] == '230202_lane8' & !is.na(mo@meta.data[['unconfined_best_match_correlation']]) & mo@meta.data[['unconfined_best_match_correlation']] == 'MO203', 'sample_final'] <- 'MO97'
mo <- add_inflammation_status(mo, sample_sheet = inflammation_assignments, 
                              seurat_lane_column = 'lane', 
                              sheet_lane_column = 'lane', 
                              seurat_participant_column = 'sample_final',
                              sheet_participants_column='sample', 
                              seurat_inflammation_column='final_condition', 
                              sheet_inflammation_column='condition')


# save the result
mo_assigned_loc <- paste(seurat_objects_loc, 'mo_all_20240223_seuratv5_assigned.rds', sep = '')
saveRDS(mo, mo_assigned_loc)
