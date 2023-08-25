#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_cellranger_summaries.R
# Function: get the cellranger summaries into an Excel file
############################################################################################################################

####################
# libraries        #
####################

library(xlsx)


####################
# Functions        #
####################

#' get the cellranger summary data per modality for one lane
#' 
#' @param mo_summary_loc the location of the summary file
#' @param joint_columns which columns to extract for the joint summary
#' @param atac_columns which columns to extract for the atac summary
#' @param rna_columns which columns to extract for the atac summary
#' @returns a list with the summary per modality, for the given file
#' summary_230105_lane1 <- read_values_summaries("/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/230105_lane1/outs/summary.csv")
read_values_summaries <- function(mo_summary_loc,
                                  joint_columns=c('Sample.ID', 'Feature.linkages.detected', 'Linked.genes', 'Linked.peaks', 'Estimated.number.of.cells'),
                                  atac_columns=c('Sample.ID', 'ATAC.Median.high.quality.fragments.per.cell', 'ATAC.Mean.raw.read.pairs.per.cell', 'ATAC.Fraction.of.high.quality.fragments.in.cells', 'ATAC.Fraction.of.transposition.events.in.peaks.in.cells', 'ATAC.Fraction.of.high.quality.fragments.overlapping.TSS', 'ATAC.Fraction.of.high.quality.fragments.overlapping.peaks', 'ATAC.Percent.duplicates'),
                                  rna_columns=c('Sample.ID', 'GEX.Median.genes.per.cell', 'GEX.Mean.raw.reads.per.cell', 'GEX.Fraction.of.transcriptomic.reads.in.cells', 'GEX.Median.UMI.counts.per.cell', 'GEX.Median.genes.per.cell', 'GEX.Percent.duplicates')) {
  # read the file
  summary_file <- read.table(mo_summary_loc, sep = ',', header = T)
  # get which columns are present
  columns_present <- colnames(summary_file)
  # warn if anything is missing
  if (length(setdiff(c(joint_columns, atac_columns, rna_columns), columns_present)) > 0) {
    warning(paste('requesting columns not present:', paste(setdiff(c(joint_columns, atac_columns, rna_columns), columns_present), collapse = ',')))
  }
  # extract what we want
  joint_data <- summary_file[, intersect(joint_columns, columns_present)]
  atac_data <- summary_file[, intersect(atac_columns, columns_present)]
  rna_data <- summary_file[, intersect(rna_columns, columns_present)]
  # put it into a list
  data_per_set <- list('joint' = joint_data, 'atac' = atac_data, 'rna' = rna_data)
  return(data_per_set)
}


#' get the cellranger summary data per modality for all the lanes
#' 
#' @param lanes_loc the folder containing the cellranger output per lane
#' @param lanes vector of lanes to get summaries for
#' @param lane_append what to add to the lane directory to get the file
#' @param joint_columns which columns to extract for the joint summary
#' @param atac_columns which columns to extract for the atac summary
#' @param rna_columns which columns to extract for the atac summary
#' @returns a list with the summaries per modalities, for all lanes
#' summaries_per_modality <- read_lane_summaries(lanes_loc)
read_lane_summaries <- function(lanes_loc, 
                                lanes=NULL, 
                                lane_append='outs/summary.csv', 
                                joint_columns=c('Sample.ID', 'Feature.linkages.detected', 'Linked.genes', 'Linked.peaks', 'Estimated.number.of.cells'),
                                atac_columns=c('Sample.ID', 'ATAC.Median.high.quality.fragments.per.cell', 'ATAC.Mean.raw.read.pairs.per.cell', 'ATAC.Fraction.of.high.quality.fragments.in.cells', 'ATAC.Fraction.of.transposition.events.in.peaks.in.cells', 'ATAC.Fraction.of.high.quality.fragments.overlapping.TSS', 'ATAC.Fraction.of.high.quality.fragments.overlapping.peaks', 'ATAC.Percent.duplicates'),
                                rna_columns=c('Sample.ID', 'GEX.Median.genes.per.cell', 'GEX.Mean.raw.reads.per.cell', 'GEX.Fraction.of.transcriptomic.reads.in.cells', 'GEX.Median.UMI.counts.per.cell', 'GEX.Median.genes.per.cell', 'GEX.Percent.duplicates')) {
  # check which lanes to use
  lanes_to_use <- lanes
  # is we didn't specify, we will try all in the directory
  if (is.null(lanes_to_use)) {
    # get the lanes
    lanes_to_use <- list.dirs(lanes_loc, full.names = F, recursive = F)
  }
  # create a list of the modalities
  summaries_per_modality <- list()
  # check each lane
  for (lane in lanes_to_use) {
    # paste together the full path
    full_summary_path <- paste(lanes_loc, '/', lane, '/', lane_append, sep = '')
    # check if the file exists
    if (file.exists(full_summary_path)) {
      # read the file
      full_summary <- read_values_summaries(full_summary_path, joint_columns = joint_columns, atac_columns = atac_columns, rna_columns = rna_columns)
      # paste it all together
      for (modality in names(full_summary)) {
        # check if the modality exists already
        if (!(modality %in% names(summaries_per_modality))) {
          # otherwise we need to create it
          summaries_per_modality[[modality]] <- list()
        }
        # add it to the list
        summaries_per_modality[[modality]][[lane]] <- full_summary[[modality]]
      }
    }
    else{
      warning(paste('summary missing for', lane, ': skipped'))
    }
  }
  # create new list with a joint table for each modality
  table_per_modality <- list()
  # now join the tables for each modality
  for (modality in names(summaries_per_modality)) {
    table_per_modality[[modality]] <- do.call('rbind', summaries_per_modality[[modality]])
  }
  return(table_per_modality)
}

#' write the cellranger summary data to an Excel file
#' 
#' @param lanes_loc the folder containing the cellranger output per lane
#' @param excel_loc where to write the Excel file
#' @param lanes vector of lanes to get summaries for
#' @param lane_append what to add to the lane directory to get the file
#' @param joint_columns which columns to extract for the joint summary
#' @param atac_columns which columns to extract for the atac summary
#' @param rna_columns which columns to extract for the atac summary
#' @returns 0 if succesfull
#' summaries_to_excel(lanes_loc, excel_loc)
summaries_to_excel <- function(lanes_loc, 
                               excel_loc,
                               lanes=NULL, 
                               lane_append='outs/summary.csv', 
                               joint_columns=c('Sample.ID', 'Feature.linkages.detected', 'Linked.genes', 'Linked.peaks', 'Estimated.number.of.cells'),
                               atac_columns=c('Sample.ID', 'ATAC.Median.high.quality.fragments.per.cell', 'ATAC.Mean.raw.read.pairs.per.cell', 'ATAC.Fraction.of.high.quality.fragments.in.cells', 'ATAC.Fraction.of.transposition.events.in.peaks.in.cells', 'ATAC.Fraction.of.high.quality.fragments.overlapping.TSS', 'ATAC.Fraction.of.high.quality.fragments.overlapping.peaks', 'ATAC.Percent.duplicates'),
                               rna_columns=c('Sample.ID','GEX.Median.genes.per.cell', 'GEX.Mean.raw.reads.per.cell', 'GEX.Fraction.of.transcriptomic.reads.in.cells', 'GEX.Median.UMI.counts.per.cell', 'GEX.Median.genes.per.cell', 'GEX.Percent.duplicates')) {
  # get the data per modality
  table_per_modality <- read_lane_summaries(lanes_loc, lanes = lanes,  lane_append = lane_append, joint_columns = joint_columns, atac_columns = atac_columns, rna_columns = rna_columns)
  # create the workbook
  wb <- createWorkbook()
  # add each sheet
  for (modality in names(table_per_modality)) {
    sheet <- createSheet(wb, modality)
    # add to the sheet
    addDataFrame(table_per_modality[[modality]], sheet = sheet, startColumn = 1, row.names = F)
  }
  # write the excel
  saveWorkbook(wb, excel_loc)
  return(0)
}


####################
# Main Code        #
####################

# location of the lanes
lanes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/'

# location for the Excel
excel_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_arc_summary.xlsx'

# write it
summaries_to_excel(lanes_loc, excel_loc)
