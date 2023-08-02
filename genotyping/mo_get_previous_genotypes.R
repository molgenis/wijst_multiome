#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Maryna Korshevniuk
# Name: mo_get_previous_genotypes.R
# Function: update the sample sheet to the correct format, collect previously genotyped individuals, create age/sex metadata
############################################################################################################################

####################
# libraries        #
####################


####################
# Functions        #
####################

#' split the rows up where the lane column has multiple entires
#' 
#' @param unsplit_table the table where the lane column can have multiple lanes
#' @param lane_column which column contains the lanes
#' @param split the character to split by
#' @returns the same table, but with the lane column only containing one lane, with a row for each lane
#' sample_to_lane_condition <- split_lanes(unsplit_table  = sample_to_lane_condition_unsplit, lane_column = 'Lane', split = ', ')
split_lanes <- function(unsplit_table, lane_column='lane', split=',') {
  # check how many original rows we have
  original_nrow <- nrow(unsplit_table)
  # get the original columns
  column_names <- colnames(unsplit_table)
  # check how many columns there were
  original_ncol <- length(column_names)
  # create a new dataframe, we'll start with twice the size of the original
  split_table <- data.frame(matrix(data = NA, nrow = original_nrow * 2, ncol = original_ncol))
  # set the column names
  colnames(split_table) <- column_names
  # set the new index
  row_i_new <- 1
  # check each of the rows in the original table
  for (row_i_original in 1:original_nrow) {
    # get the lane value
    lane_unsplit_value <- unsplit_table[row_i_original, lane_column]
    # split by the split character
    lane_split <- strsplit(lane_unsplit_value, split)[[1]]
    # check each in the split
    for (individual_lane in lane_split) {
      # check if the new table is still large enough
      if (row_i_new > nrow(split_table)) {
        # we will double it in size
        split_table_added <- data.frame(matrix(data = NA, nrow = nrow(split_table), ncol = ncol(split_table)))
        colnames(split_table_added) <- column_names
        split_table <- rbind(split_table, split_table_added)
      }
      # set the original row first
      split_table[row_i_new, column_names] <- unsplit_table[row_i_original, column_names]
      # but replace the lane entry
      split_table[row_i_new, lane_column] <- individual_lane
      # update the index
      row_i_new <- row_i_new + 1
    }
  }
  # remove where we don't have the lane
  split_table <- split_table[!is.na(split_table[[lane_column]]), ]
  return(split_table)
}

#' for entries in the table missing the GSA_ID, see if the table contains another entry for that sample, where the GSA_ID was given, and use that GSA_ID
#' 
#' @param sample_condition_table the table containing the samples and the GSA IDs
#' @param gsa_column the column containing the GSA ID
#' @param sample_column the column containing the sample name
#' @returns the same table, with the GSA IDs filled where possible
#' sample_to_lane_condition <- get_matched_gsa_id(sample_to_lane_condition)
get_matched_gsa_id <- function(sample_condition_table, gsa_column='GSA_ID', sample_column='Sample') {
  # get the gsa IDs we have
  gsa_ids <- unique(sample_condition_table[[gsa_column]])
  # removing the NA ones
  gsa_ids <- gsa_ids[!is.na(gsa_ids)]
  # check each row
  for (gsa_id in gsa_ids) {
    # get the sample for that gsa IDs
    rows_sample <- sample_condition_table[!is.na(sample_condition_table[[gsa_column]]) & sample_condition_table[[gsa_column]] == gsa_id, ]
    # we have to have a sample
    if (nrow(rows_sample) > 0) {
      # check the samples
      samples <- rows_sample[[sample_column]]
      # filter
      samples <- samples[!is.na(samples)]
      samples <- unique(samples)
      # if we have a sample, then we can add the gsa where it is missing
      if (length(samples) > 0) {
        sample_condition_table[!is.na(sample_condition_table[[sample_column]]) &
                                 is.na(sample_condition_table[[gsa_column]]) &
                                 sample_condition_table[[sample_column]] == samples[1],
                               gsa_column] <- gsa_id
      }
    }
  }
  return(sample_condition_table)
}

#' for samples from the LifeLines DEEP cohort, add the 8-number pseudo ID
#' 
#' @param sample_condition_table the table containing the samples and the unformatted LifeLines DEEP IDs
#' @param deep_to_pseudo_id the table mapping LifeLines DEEP IDs to the the 8-number pseudo ID
#' @param deep_column the column containing the sample name in the unformatted LifeLines DEEP format
#' @param genotype_column the column to place the 8-number pseudo ID
#' @param deep_to_pseudo_id_deep_column the column in deep_to_pseudo_id that contains the formatted LifeLines DEEP ID
#' @param deep_to_pseudo_id_pseudo_column the column in deep_to_pseudo_id that contains the 8-number pseudo ID
#' @returns the same table, with the formatted LineLines DEEP name in 'deep_safe_column', and the 'genotype_column' updated with the 8-number pseudo ID where possible
#' sample_to_lane_condition <- add_sample_id_for_deep(sample_to_lane_condition, lldeep_to_psuedo_mapping)
add_sample_id_for_deep <- function(sample_condition_table, deep_to_pseudo_id, deep_column='Sample', genotype_column='Genotype', deep_to_pseudo_id_deep_column='LLDEEP_ID', deep_to_pseudo_id_pseudo_column='PSEUDOIDEXT') {
  # we need to pad zeroes to the beginning of the NEXT IDs
  sample_condition_table[['deep_safe_column']] <- apply(sample_condition_table, 1, function(x){
    # extract the value from the right column
    unpadded_id <- x[deep_column]
    # first we check if it is a DEEP ID
    if (grepl('DEEP-', unpadded_id)) {
      # split on the dash
      unpadded_id_split <- strsplit(unpadded_id, '-')[[1]]
      # we need the second entry
      deep_id_number <- unpadded_id_split[2]
      # we need four characters, so the current length minus that is how much to pad
      nr_to_pad <- 4 - nchar(deep_id_number)
      # pad some zeroes
      pad <- ''
      if (nr_to_pad > 0) {
        pad <- paste(rep('0', times = nr_to_pad), collapse = '')
      }
      # now add that to the original, plus the actual ID start string
      padded_id <- paste('LLDeep_', pad, deep_id_number, sep = '')
    }
    else {
      return('')
    }
  })
  # now add the 'Genotype' for those the entries that have this deep column
  sample_condition_table[!is.na(sample_condition_table[['deep_safe_column']]) & sample_condition_table[['deep_safe_column']] != '', genotype_column] <- deep_to_pseudo_id[match(sample_condition_table[!is.na(sample_condition_table[['deep_safe_column']]) & sample_condition_table[['deep_safe_column']] != '', 'deep_safe_column'], deep_to_pseudo_id[[deep_to_pseudo_id_deep_column]]), deep_to_pseudo_id_pseudo_column]
  return(sample_condition_table)
}


####################
# Main Code        #
####################

# the per lane samples
sample_to_lane_condition_unsplit_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_sample_to_condition.tsv'
# the mapping of the missing genotypes to UGLI IDs
sample_to_ugli_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_id_to_otherids.tsv'
# lifelines COVID IDs to project pseudo IDs
lifelines_mapping_loc <- '/groups/umcg-lifelines/tmp01/projects/ov22_0441/ongoing/metadata/OV20_00100_3OR_samples_linkage_file_20210609_merged.csv'
# pseudo into to ext
pseudo_int_to_ext_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/phenotype_linkage_file_project_pseudo_id.txt'
# the UGLI mapping
pseudo_int_to_ugli_mapping_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'
# and finally the DEEP mapping
lldeep_to_psuedo_mapping_loc <- '/groups/umcg-lifelines/rsc01/releases/deep_linkage_files/v1/DEEP_linkage_file_cluster.dat'

# read the files
sample_to_lane_condition_unsplit <- read.table(sample_to_lane_condition_unsplit_loc, header = T, sep = '\t')
sample_to_ugli <- read.table(sample_to_ugli_loc, header = T, sep = '\t')
lifelines_mapping <- read.table(lifelines_mapping_loc, sep = ',', header = T, comment.char = '')
lldeep_to_psuedo_mapping <- read.table(lldeep_to_psuedo_mapping_loc, sep = '\t', header = T, comment.char = '')
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T)
pseudo_int_to_ugli_mapping <- read.table(pseudo_int_to_ugli_mapping_loc, sep = '\t', header = T)
# the ugli mapping columns seem to be missing
colnames(pseudo_int_to_ugli_mapping) <- c('PSEUDOIDEXT', 'UGLI_ID', 'DNA')

# drop the sampling date
lifelines_mapping[['sampledate']] <- NULL
# and name
lifelines_mapping[['name']] <- NULL
# and sample type
lifelines_mapping[['sampletype']] <- NULL
# so we can make it so we only have one entry per participant
lifelines_mapping <- unique(lifelines_mapping)


# split the lanes up
sample_to_lane_condition <- split_lanes(unsplit_table  = sample_to_lane_condition_unsplit, lane_column = 'Lane', split = ', ')
# add the gsa ID where it is missing
sample_to_lane_condition <- get_matched_gsa_id(sample_to_lane_condition)
# now add the sample to UGLI mapping
sample_to_lane_condition <- merge(x = sample_to_lane_condition, y = sample_to_ugli, by.x = 'Sample', by.y = 'Other.ID', all = T)
# and the lifelines mapping
sample_to_lane_condition <- merge(x = sample_to_lane_condition, y = lifelines_mapping, by.x = 'Sample', by.y = 'projectparticipant.', all = T)
# and the pseudo ext, sanity check. Here the 'Genotype' column should match the 'll_pseudo_ext' column
sample_to_lane_condition[['ll_pseudo_ext']] <- pseudo_int_to_ext[match(sample_to_lane_condition[['project_pseudo_id']], pseudo_int_to_ext[['PROJECT_PSEUDO_ID']]), 'PSEUDOIDEXT']
# add the Genotype, so the pseudoextid, based on the deep ID
sample_to_lane_condition <- add_sample_id_for_deep(sample_to_lane_condition, lldeep_to_psuedo_mapping)
# seems good, so do the opposite as well. Which also looks good
sample_to_lane_condition[['ll_pseudo_int']] <- pseudo_int_to_ext[match(sample_to_lane_condition[['Genotype']], pseudo_int_to_ext[['PSEUDOIDEXT']]), 'PROJECT_PSEUDO_ID']
# now add the UGLI ID, based on the PSUEDOIEXT
sample_to_lane_condition[['UGLI_ID']] <- pseudo_int_to_ugli_mapping[match(sample_to_lane_condition[['Genotype']], pseudo_int_to_ugli_mapping[['PSEUDOIDEXT']]), 'UGLI_ID']
# remove entries we are not using
sample_to_lane_condition <- sample_to_lane_condition[!is.na(sample_to_lane_condition[['Lane']]), ]
# turn the '<NA>' into a ''
sample_to_lane_condition[is.na(sample_to_lane_condition[['availability']]) | sample_to_lane_condition[['availability']] == '<NA>', 'availability'] <- ''
sample_to_lane_condition[is.na(sample_to_lane_condition[['UGLI_ID']]) | sample_to_lane_condition[['UGLI_ID']] == '<NA>', 'UGLI_ID'] <- ''
# so we can easily check which annotation is more informative
sample_to_lane_condition[sample_to_lane_condition[['availability']] == '' & sample_to_lane_condition[['UGLI_ID']] != '', ]
sample_to_lane_condition[sample_to_lane_condition[['UGLI_ID']] == '' & sample_to_lane_condition[['availability']] != '', ]
# usually the 'availability' is more informative, but the LLDEEP, we need the UGLI_ID one. We will make joint column
sample_to_lane_condition[['UGLI_FINAL']] <- sample_to_lane_condition[['availability']]
sample_to_lane_condition[!is.na(sample_to_lane_condition[['UGLI_ID']]) & sample_to_lane_condition[['UGLI_ID']] != '', 'UGLI_FINAL'] <- sample_to_lane_condition[!is.na(sample_to_lane_condition[['UGLI_ID']]) & sample_to_lane_condition[['UGLI_ID']] != '', 'UGLI_ID']
# see if we are missing anything still 
sample_to_lane_condition[sample_to_lane_condition[['UGLI_FINAL']] == '' & is.na(sample_to_lane_condition[['GSA_ID']]), ]

# check which one we are missing
for (id in unique(sample_to_lane_condition[
  (is.na(sample_to_lane_condition[['UGLI_FINAL']]) | sample_to_lane_condition[['UGLI_FINAL']] == '') & 
  (is.na(sample_to_lane_condition[['GSA_ID']]) | is.na(sample_to_lane_condition[['GSA_ID']]) == '')
  , 'Sample'])) {
  print(id)
}

# create a new final column
sample_to_lane_condition[['sample']] <- NA
# update the GSA ID to have the actual GSA format
sample_to_lane_condition[!is.na(sample_to_lane_condition[['GSA_ID']]) & sample_to_lane_condition[['GSA_ID']] != '', 'GSA_ID'] <- paste('MO', sample_to_lane_condition[!is.na(sample_to_lane_condition[['GSA_ID']]) & sample_to_lane_condition[['GSA_ID']] != '', 'GSA_ID'], sep = '')
# use gsa where possible
sample_to_lane_condition[!is.na(sample_to_lane_condition[['GSA_ID']]) & sample_to_lane_condition[['GSA_ID']] != '', 'sample'] <- sample_to_lane_condition[!is.na(sample_to_lane_condition[['GSA_ID']]) & sample_to_lane_condition[['GSA_ID']] != '', 'GSA_ID']
# UGLI where we can
sample_to_lane_condition[!is.na(sample_to_lane_condition[['UGLI_FINAL']]) & sample_to_lane_condition[['UGLI_FINAL']] != '', 'sample'] <- sample_to_lane_condition[!is.na(sample_to_lane_condition[['UGLI_FINAL']]) & sample_to_lane_condition[['UGLI_FINAL']] != '', 'UGLI_FINAL']

# subset to what we need from a sample sheet
sample_sheet <- sample_to_lane_condition[, c('Lane', 'Condition', 'sample')]
# rename columns
colnames(sample_sheet) <- c('lane', 'condition', 'sample')
# write the resulting sample sheet
write.table(sample_sheet, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv', sep = '\t', row.names = F, col.names = T, quote = F)

# get a list of the UGLI participants, to use as a subsetting list for the UGLI genotype data
write.table(
  data.frame(
    x = unique(sample_to_lane_condition[['UGLI_FINAL']][!is.na(sample_to_lane_condition[['UGLI_FINAL']])]),
    y = unique(sample_to_lane_condition[['UGLI_FINAL']][!is.na(sample_to_lane_condition[['UGLI_FINAL']])])
  ), '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_ugli_plinkfilter.tsv', row.names = F, col.names = F, quote = F)

# create an age+sex file for the UGLI participants
age_sex_ugli <- unique(sample_to_lane_condition[!is.na(sample_to_lane_condition[['UGLI_FINAL']]) & sample_to_lane_condition[['UGLI_FINAL']] != '',  c('UGLI_FINAL', 'Seks', 'Age')])
# set other column names
colnames(age_sex_ugli) <- c('sample', 'sex', 'age')
# write the table
write.table(age_sex_ugli, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_age_sex_ugli.tsv', sep = '\t', row.names = F, col.names = T, quote = F)

# also one for the new GSA participants
age_sex_mo <- unique(sample_to_lane_condition[!is.na(sample_to_lane_condition[['GSA_ID']]) & sample_to_lane_condition[['GSA_ID']] != '',  c('GSA_ID', 'Seks', 'Age')])
colnames(age_sex_mo) <- c('sample', 'sex', 'age')
write.table(age_sex_mo, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_age_sex_batch12.tsv', sep = '\t', row.names = F, col.names = T, quote = F)
