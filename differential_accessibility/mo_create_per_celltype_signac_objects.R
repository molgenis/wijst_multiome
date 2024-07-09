#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_per_celltype_signac_objects.R
# Function: 
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(tools)
library(doParallel)


####################
# Functions        #
####################


create_anonymized_mapping <- function(metadata, source_column, target_column, target_prepend='') {
  # get the unique values
  source_values <- unique(as.character(metadata[[source_column]]))
  # turn NA into 'unknown'
  source_values[is.na(source_values)] <- 'unknown'
  # randomly sort your values, by sampling all values
  source_values <- source_values[sample(1:length(source_values), length(source_values))]
  # create a mapping
  target_values <- paste(target_prepend, 1:length(source_values), sep = '')
  # now put that into a dataframe
  mapping_table <- data.frame(x = source_values, y = target_values)
  # now make the column names as expected
  colnames(mapping_table) <- c(source_column, target_column)
  return(mapping_table)
}

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



process_and_write_subset <- function(signac_object, correlations_confined, correlations_unconfined, condition_assignments, age_sex_assigments, realid_assignments, longcovid_assignments, output_loc, only_covid=F, recode_id=F) {
  # add the confined sample
  signac_object <- AddMetaData(signac_object, correlations_confined[, colnames(correlations_confined)])
  # and for unconfined
  signac_object <- AddMetaData(signac_object, correlations_unconfined[, colnames(correlations_unconfined)])
  # set the final assignments
  signac_object@meta.data[['sample_final']] <- NA
  # set those with high enough correlations
  signac_object@meta.data[!is.na(signac_object@meta.data[['unconfined_best_match_correlation']]) & signac_object@meta.data[['unconfined_best_match_correlation']] > 0.5, 'sample_final'] <- signac_object@meta.data[!is.na(signac_object@meta.data[['unconfined_best_match_correlation']]) & signac_object@meta.data[['unconfined_best_match_correlation']] > 0.5, 'unconfined_best_match_sample']
  # this one has the best correlation for MO408 and MO1001. MO1001 should be present, so that one we will assign it to
  signac_object@meta.data[signac_object@meta.data[['lane']] == '230112_lane1' & !is.na(signac_object@meta.data[['unconfined_best_match_correlation']]) & signac_object@meta.data[['unconfined_best_match_correlation']] == 'MO408', 'sample_final'] <- 'MO1001'
  # this one has the best correlation for MO97 and MO203. MO203 should be present, so that one we will assign it to
  signac_object@meta.data[signac_object@meta.data[['lane']] == '230202_lane8' & !is.na(signac_object@meta.data[['unconfined_best_match_correlation']]) & signac_object@meta.data[['unconfined_best_match_correlation']] == 'MO203', 'sample_final'] <- 'MO97'
  # do the condition assignment
  signac_object <- add_inflammation_status(
    signac_object,
    sample_sheet=condition_assignments, 
    seurat_lane_column='lane',
    sheet_lane_column='lane', 
    seurat_participant_column='sample_final', 
    sheet_participants_column='sample_final', 
    seurat_inflammation_column='condition_sheet', 
    sheet_inflammation_column='condition'
  )
  signac_object <- add_inflammation_status(
    signac_object,
    sample_sheet=condition_assignments, 
    seurat_lane_column='lane',
    sheet_lane_column='lane', 
    seurat_participant_column='sample_final', 
    sheet_participants_column='sample', 
    seurat_inflammation_column='condition_prev', 
    sheet_inflammation_column='cond_prev'
  )
  # rename 24hCA
  signac_object@meta.data[!is.na(signac_object@meta.data[['condition_sheet']]) &
                            signac_object@meta.data[['condition_sheet']] == '24hCa', 'condition_sheet'] <- '24hCA'
  # now also set the final inflammation assignment
  signac_object@meta.data[['condition_final']] <- signac_object@meta.data[['condition_sheet']]
  signac_object@meta.data[is.na(signac_object@meta.data[['condition_final']]), 'condition_final'] <- signac_object@meta.data[is.na(signac_object@meta.data[['condition_final']]), 'condition_prev']
  
  # remove any empty entries
  signac_object <- signac_object[, !is.na(signac_object@meta.data[['sample_final']]) &
             !is.na(signac_object@meta.data[['condition_final']])]
  
  # add the realid
  signac_object@meta.data['realid'] <- realid_assignments[match(signac_object@meta.data$sample_final, realid_assignments$sample_final), 'realid_final']
  # add the age and sex
  signac_object@meta.data[['age']] <- age_sex_assigments[match(signac_object@meta.data$realid, age_sex_assigments$sample), 'age']
  signac_object@meta.data[['sex']] <- age_sex_assigments[match(signac_object@meta.data$realid, age_sex_assigments$sample), 'sex']
  
  # add annotated LONG_COVID assignment
  # add the first part to the name of the long covid ID
  longcovid_assignments[['mo']] <- paste('MO', longcovid_assignments$Project_ID, sep = '')
  # add the longcovid assignment first
  signac_object@meta.data[['LONG_COVID']] <- NA
  # now add the ones we have
  signac_object@meta.data[signac_object@meta.data[['sample_final']] %in% longcovid_assignments[['mo']], 'LONG_COVID'] <- longcovid_assignments[match(signac_object@meta.data[signac_object@meta.data[['sample_final']] %in% longcovid_assignments[['mo']], 'sample_final'], longcovid_assignments[['mo']]), 'Case.Control']
  # set empty to NA
  signac_object@meta.data[!is.na(signac_object@meta.data[['LONG_COVID']]) & signac_object@meta.data[['LONG_COVID']] == '', 'LONG_COVID'] <- NA
  # add the final LONG_COVID assignment
  signac_object@meta.data[['LONG_COVID_final']] <- signac_object@meta.data$LONG_COVID
  signac_object@meta.data[is.na(signac_object@meta.data[['LONG_COVID']]), 'LONG_COVID_final'] <- 'control'
  signac_object@meta.data[['LONG_COVID_method']] <- NA
  signac_object@meta.data[is.na(signac_object@meta.data[['LONG_COVID']]), 'LONG_COVID_method'] <- 'inferred'
  signac_object@meta.data[!is.na(signac_object@meta.data[['LONG_COVID']]), 'LONG_COVID_method'] <- 'assigned'
  # check if only covid
  if (only_covid) {
    signac_object <- signac_object[, 
                                   !is.na(signac_object@meta.data[['condition_final']]) &
                                     signac_object@meta.data[['condition_final']] == 'UT' &
                                     !is.na(signac_object@meta.data[['LONG_COVID_final']])]
  }
  if (recode_id) {
    # create anonimized mapping
    p_mapping <- create_anonymized_mapping(signac_object@meta.data, 'sample_final', 'sample_number', 's')
    # add these new columns
    signac_object@meta.data[['sample_final']] <- p_mapping[match(signac_object@meta.data[['sample_final']], p_mapping[['sample_final']]), 'sample_number']
    signac_object@meta.data[['best_match_sample']] <- p_mapping[match(signac_object@meta.data[['best_match_sample']], p_mapping[['sample_final']]), 'sample_number']
    signac_object@meta.data[['second_match_sample']] <- p_mapping[match(signac_object@meta.data[['second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
    signac_object@meta.data[['confined_best_match_sample']] <- p_mapping[match(signac_object@meta.data[['confined_second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
    signac_object@meta.data[['confined_second_match_sample']] <- p_mapping[match(signac_object@meta.data[['confined_second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
    signac_object@meta.data[['unconfined_best_match_sample']] <- p_mapping[match(signac_object@meta.data[['unconfined_best_match_sample']], p_mapping[['sample_final']]), 'sample_number']
    signac_object@meta.data[['unconfined_second_match_sample']] <- p_mapping[match(signac_object@meta.data[['unconfined_second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
    # write the sample mapping
    write.table(p_mapping, gzfile(paste(output_loc, '.sample_mapping.tsv.gz', sep = ''), sep = '\t', row.names = F, col.names = T))
    # remove all other IDs
    signac_object@meta.data[, c('realid')] <- NULL
  }
  # save the result
  saveRDS(signac_object, output_loc)
  # make md5
  md5_object <- md5sum(output_loc)[[1]]
  # write the md5
  write.table(data.frame(x = c(md5_object)), paste(output_loc, '.md5', sep = ''), row.names = F, col.names = F, quote = F)
}



####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 185 * 2000 * 1024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# correlation files
correlation_mapping_per_barcode_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched.tsv'
correlation_mapping_per_barcode_all_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched_vs_all.tsv'
# get the assignment matrices
correlation_mapping_per_barcode <- read.table(correlation_mapping_per_barcode_loc, header = T, sep = '\t')
correlation_mapping_per_barcode_all <- read.table(correlation_mapping_per_barcode_all_loc, header = T, sep = '\t')
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

# location of the condition assignment
condition_assignment_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_monocyte_based_condition_numbers.tsv'
# read the conditions
condition_assignments <- read.table(condition_assignment_loc, header = T, sep = '\t')
# location of the age/sex assignments
age_sex_assignments_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_age_sex_batch12.tsv'
# read the assignments
age_sex_assigments <- read.table(age_sex_assignments_loc, header = T, sep = '\t')
# the ugli ones as well
age_sex_assigments_ugli_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_age_sex_ugli.tsv'
age_sex_assigments_ugli <- read.table(age_sex_assigments_ugli_loc, header = T, sep = '\t')
# merge
age_sex_assigments <- rbind(age_sex_assigments, age_sex_assigments_ugli)
# and the 'realids'
realid_assignments_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'
# read the realids
realid_assignments <- read.table(realid_assignments_loc, header = T, sep = '\t')
# fix the sample IDs in the realid table
realid_assignments[['sample_final']] <- realid_assignments$SampleID_CORRECT
realid_assignments[is.na(realid_assignments$SampleID_CORRECT) | realid_assignments$SampleID_CORRECT == '', 'sample_final']  <- realid_assignments[is.na(realid_assignments$SampleID_CORRECT) | realid_assignments$SampleID_CORRECT == '', 'sample']
realid_assignments[['realid_final']] <- realid_assignments$RealID
realid_assignments[is.na(realid_assignments$RealID) | realid_assignments$RealID == '', 'realid_final']  <- realid_assignments[is.na(realid_assignments$RealID) | realid_assignments$RealID == '', 'sample_final']

# location of the LONG-CoVID assignments
longcovid_assignments_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_longcovid_assignments.tsv'
# read the assignments
longcovid_assignments <- read.table(longcovid_assignments_loc, header = T, sep = '\t')

# location of the full object
mo_all_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltypemajor_wstatus_1_80.rds'
# read object
mo_all <- readRDS(mo_all_loc)
# subset to monocytes, because we can't handle all of it at once
mo_mono <- mo_all[['monocyte']]
# remove the full set to make sure we have enough memory
rm(mo_all)
# add celltype as explicit column
mo_mono@meta.data[['cell_type']] <- 'monocyte'
# do this cell type
process_and_write_subset(
  mo_mono,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments, 
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_monocyte_wstatus_1_80_20240709.rds'
)
process_and_write_subset(
  mo_mono,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_monocyte_UT_wstatus_1_80_20240709.rds',
  only_covid = T,
  recode_id = T
)
# clear memory
rm(mo_mono)

# same for CD4T
mo_all <- readRDS(mo_all_loc)
mo_cd4t <- mo_all[['CD4T']]
rm(mo_all)
mo_cd4t@meta.data[['cell_type']] <- 'CD4T'
process_and_write_subset(
  mo_cd4t,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_cd4t_wstatus_1_80_20240709.rds'
)
process_and_write_subset(
  mo_cd4t,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_cd4t_UT_wstatus_1_80_20240709.rds',
  only_covid = T, 
  recode_id = T
)
rm(mo_cd4t)
# and CD8T
mo_all <- readRDS(mo_all_loc)
mo_cd8t <- mo_all[['CD8T']]
rm(mo_all)
mo_cd8t@meta.data[['cell_type']] <- 'CD8T'
process_and_write_subset(
  mo_cd8t,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_cd8t_wstatus_1_80_20240709.rds'
)
process_and_write_subset(
  mo_cd8t,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_cd8t_UT_wstatus_1_80_20240709.rds',
  only_covid = T, 
  recode_id = T
)
rm(mo_cd8t)
# and NK
mo_all <- readRDS(mo_all_loc)
mo_nk <- mo_all[['NK']]
rm(mo_all)
mo_nk@meta.data[['cell_type']] <- 'NK'
process_and_write_subset(
  mo_nk,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_nk_wstatus_1_80_20240709.rds'
)
process_and_write_subset(
  mo_nk,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_nk_UT_wstatus_1_80_20240709.rds',
  only_covid = T, 
  recode_id = T
)
rm(mo_nk)
# and B
mo_all <- readRDS(mo_all_loc)
mo_b <- mo_all[['B']]
rm(mo_all)
mo_b@meta.data[['cell_type']] <- 'B'
process_and_write_subset(
  mo_b,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_b_wstatus_1_80_20240709.rds'
)
process_and_write_subset(
  mo_b,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_b_UT_wstatus_1_80_20240709.rds',
  only_covid = T, 
  recode_id = T
)
rm(mo_b)
# and finally DC
mo_all <- readRDS(mo_all_loc)
mo_dc <- mo_all[['DC']]
rm(mo_all)
mo_dc@meta.data[['cell_type']] <- 'DC'
process_and_write_subset(
  mo_dc,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_dc_wstatus_1_80_20240709.rds'
)
process_and_write_subset(
  mo_dc,
  correlations_confined = correlations_confined,
  correlations_unconfined = correlations_unconfined,
  condition_assignments=condition_assignments, 
  age_sex_assigments=age_sex_assigments, 
  realid_assignments=realid_assignments,  
  longcovid_assignments=longcovid_assignments,
  output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_dc_UT_wstatus_1_80_20240709.rds',
  only_covid = T, 
  recode_id = T
)
rm(mo_dc)