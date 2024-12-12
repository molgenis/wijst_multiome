#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_n_cellss_accessible_tables.R
# Function: create tables that have the number of cells per sample that have any reads for a region
############################################################################################################################

####################
# libraries        #
####################

# for reading Seurat object
library(Seurat)
library(Signac)
# for saving the results
library(data.table)


####################
# Functions        #
####################


get_ncell_expressed_sample <- function(seurat_object, sample_present, sample_column='sample_lane', layer='counts', assay='RNA') {
  # subset to that sample
  seurat_object_sample <- seurat_object[, !is.na(seurat_object@meta.data[[sample_column]]) & seurat_object@meta.data[[sample_column]] == sample_present]
  # initalize table
  nonzeroes_datatable <- NULL
  # and check only if we have entries
  if (ncol(seurat_object_sample) > 1 & nrow(seurat_object_sample) > 0) {
    # get the assay data
    assay_data <- GetAssayData(seurat_object_sample, layer = layer, assay = assay)
    # turn into row-wise sparse matrix
    counts_as_row_sparse <- as(assay_data, "RsparseMatrix")
    # extract the features
    features_present <- rownames(counts_as_row_sparse)
    # get the nonzero sum for each gene by doing a rowsum on the T/F value you get from the zero-or-not comparison
    non_zero_cell_nr <- as.vector(unlist(rowSums(counts_as_row_sparse != 0)))
    # make into datatable
    nonzeroes_datatable <- data.table(x = features_present, y = non_zero_cell_nr)
  }
  else if(ncol(seurat_object_sample) == 1 & nrow(seurat_object_sample) > 0) {
    warning(paste('Seurat does not play nice when subsetting using one cell, will slightly erroneously return 0 for', sample_present))
    nonzeroes_datatable <- data.table(x = character(), y = numeric())
  }
  # othewise make a dummy entry
  else {
    nonzeroes_datatable <- data.table(x = character(), y = numeric())
  }
  # set dimension names
  colnames(nonzeroes_datatable) <- c('feature', sample_present)
  # return result
  return(nonzeroes_datatable)
}


get_ncell_expressed_sample_list_multithreaded <- function(seurat_object, sample_column='sample_lane', layer='counts', assay='RNA', nthreads=4) {
  # load parallel library
  library(doParallel)
  # start cluster with threads
  clusterThreads = makeCluster(nthreads)
  # register the cluster
  registerDoParallel(clusterThreads)
  # get unique samples
  unique_samples <- unique(seurat_object@meta.data[[sample_column]])
  # remove NA
  unique_samples <- unique_samples[!is.na(unique_samples)]
  # store results per sample
  ncell_exp_per_sample = foreach(sample_present=unique_samples, .packages=c('Seurat', 'data.table')) %dopar% {
    source('mo_create_n_cellss_expressed_tables_functions.R')
    nonzeroes_datatable <- get_ncell_expressed_sample(seurat_object = seurat_object, sample_present = sample_present, sample_column = sample_column, layer = layer, assay = assay)
  }
  return(ncell_exp_per_sample)
}


get_ncell_expressed_sample_list <- function(seurat_object, sample_column='sample_lane', layer='counts', assay='RNA') {
  # store results per sample
  ncell_exp_per_sample <- list()
  # check each sample
  for (sample_present in unique(seurat_object@meta.data[[sample_column]])) {
    # subset to that sample
    seurat_object_sample <- seurat_object[, !is.na(seurat_object@meta.data[[sample_column]]) & seurat_object@meta.data[[sample_column]] == sample_present]
    # initalize table
    nonzeroes_datatable <- NULL
    # and check only if we have entries
    if (ncol(seurat_object_sample) > 1 & nrow(seurat_object_sample) > 0) {
      # get the assay data
      assay_data <- GetAssayData(seurat_object_sample, layer = layer, assay = assay)
      # turn into row-wise sparse matrix
      counts_as_row_sparse <- as(assay_data, "RsparseMatrix")
      # extract the features
      features_present <- rownames(counts_as_row_sparse)
      # get the nonzero sum for each gene by doing a rowsum on the T/F value you get from the zero-or-not comparison
      non_zero_cell_nr <- as.vector(unlist(rowSums(counts_as_row_sparse != 0)))
      # make into datatable
      nonzeroes_datatable <- data.table(x = features_present, y = non_zero_cell_nr)
    }
    else if(ncol(seurat_object_sample) == 1 & nrow(seurat_object_sample) > 0) {
      warning(paste('Seurat does not play nice when subsetting using one cell, will slightly erroneously return 0 for', sample_present))
      nonzeroes_datatable <- data.table(x = character(), y = numeric())
    }
    # othewise make a dummy entry
    else {
      nonzeroes_datatable <- data.table(x = character(), y = numeric())
    }
    # set dimension names
    colnames(nonzeroes_datatable) <- c('feature', sample_present)
    # put into list
    ncell_exp_per_sample[[sample_present]] <- nonzeroes_datatable
  }
  return(ncell_exp_per_sample)
}

get_ncell_expressed_matrix <- function(seurat_object, sample_column='sample_lane', layer='counts', assay='RNA', multithread=T, nthreads=4) {
  # init list
  ncell_exp_per_sample <- NULL
  # do single or multithreaded
  if (multithread) {
    ncell_exp_per_sample <- get_ncell_expressed_sample_list_multithreaded(seurat_object = seurat_object, sample_column = sample_column, layer = layer, assay = assay, nthreads = nthreads)
  }
  else {
    ncell_exp_per_sample <- get_ncell_expressed_sample_list(seurat_object = seurat_object, sample_column = sample_column, layer = layer, assay = assay)
  }
  # join all the datatables, all needs to be set to TRUE because features with zero across all cells are dropped, without 'all' we would only have genes present in all samples
  ncell_exp_per_sample_dt = Reduce(function(...) merge(..., all = TRUE, by = 'feature'), ncell_exp_per_sample)
  # set NA to zero. When subsetting in Seurat, features with zero across all cells are dropped. It is thus safe to assume that NAs in the joined datatable were zero for that sample
  setnafill(ncell_exp_per_sample_dt, cols = setdiff(colnames(ncell_exp_per_sample_dt), 'feature'), fill = 0)
  # return the result
  return(ncell_exp_per_sample_dt)
}

combine_lanes <- function(lanes) {
  # we'll make a mapping for the lists
  lane_remapping <- list()
  # we'll check each lane
  for (lane in lanes) {
    # we'll extract the last number from the lane
    lane_nr <- regmatches(lane,regexpr("\\d+$",lane))
    # check if it is even
    if (as.numeric(lane_nr) %% 2 == 0) {
      # remove the last number
      lane_no_nr <- substr(lane, 1, nchar(lane) -1)
      # get the previous number
      previous_lane_nr <- as.numeric(lane_nr) - 1
      # get the combined lane
      combined_lane <- paste(lane_no_nr, previous_lane_nr, lane_nr, sep = '')
      # put in list
      lane_remapping[[lane]] <- combined_lane
    }
    else if (as.numeric(lane_nr) %% 2 == 1) {
      # just add the next number
      combined_lane <- paste(lane, as.numeric(lane_nr) + 1, sep = '')
      # put in list
      lane_remapping[[lane]] <- combined_lane
    }
  }
  return(lane_remapping)
}

read_barcode_and_lane <- function(seurat_object) {
  # do the split first
  seurat_object_rowsnames_split <- strsplit(colnames(seurat_object), split = '_')
  # now do a list apply
  df_per_barcode <- lapply(seurat_object_rowsnames_split, FUN = function(x) {
    data.frame(lane = paste(x[2], x[3], sep = '_'), barcode = x[1])
  })
  # merge all of them together
  extra_metadata <- do.call('rbind', df_per_barcode)
  # set the rownames to be the original ones
  rownames(extra_metadata) <- colnames(seurat_object)
  # now add the extra data we have
  seurat_object <- AddMetaData(seurat_object, extra_metadata)
  return(seurat_object)
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


add_inflammation_status_each_object <- function(seurat_object_list, sample_sheet, seurat_lane_column='lane', sheet_lane_column='lane', seurat_participant_column='soup_final_sample_assignment', sheet_participants_column='genoid', seurat_inflammation_column='inflammation_status', sheet_inflammation_column='inflammation_status') {
  # let's go over each object
  for (object_name in names(seurat_object_list)) {
    # do the condition assignment
    seurat_object_list[[object_name]] <- add_inflammation_status(
      seurat_object_list[[object_name]],
      sample_sheet=sample_sheet, 
      seurat_lane_column=seurat_lane_column,
      sheet_lane_column=sheet_lane_column, 
      seurat_participant_column=seurat_participant_column, 
      sheet_participants_column=sheet_participants_column, 
      seurat_inflammation_column=seurat_inflammation_column, 
      sheet_inflammation_column=sheet_inflammation_column
    )
  }
  return(seurat_object_list)
}


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# where to place the outputs
out_tables_folder <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/ncells_by_regions/L1/combined/'

# location of the condition assignment
condition_assignment_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_monocyte_based_condition_numbers.tsv'
# read the conditions
condition_assignments <- read.table(condition_assignment_loc, header = T, sep = '\t')

# location of the cell type objects
cell_type_objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltypemajor_1_80.rds'

# read the object
cell_type_objects <- readRDS(cell_type_objects_loc)

# add barcodes back
for(cell_type in names(cell_type_objects)) {
  cell_type_objects[[cell_type]] <- read_barcode_and_lane(cell_type_objects[[cell_type]])
}
# get the assignment matrices
correlation_mapping_per_barcode_all <- read.table('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched_vs_all.tsv', header = T, sep = '\t')
# set barcodes and remove data we already have
rownames(correlation_mapping_per_barcode_all) <- correlation_mapping_per_barcode_all[['barcode_lane']]
correlation_mapping_per_barcode_all[, c('lane', 'barcode_lane', 'barcode', 'barcode_original')] <- NULL
# now let's get the souporcell data specifically, which would be the same for all and per-lane
soup_only <- correlation_mapping_per_barcode_all[, setdiff(colnames(correlation_mapping_per_barcode_all), c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation'))]
# and the correlation data
correlations_unconfined <- correlation_mapping_per_barcode_all[, c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')]
# add the confined sample
for(cell_type in names(cell_type_objects)) {
  cell_type_objects[[cell_type]] <- AddMetaData(cell_type_objects[[cell_type]], correlations_unconfined[, colnames(correlations_unconfined)])
}

# add the conditions of the original sheet
cell_type_objects <- add_inflammation_status_each_object(cell_type_objects, condition_assignments, seurat_participant_column='best_match_sample', sheet_participants_column = 'sample_final', seurat_inflammation_column = 'inflammation_sheet', sheet_inflammation_column = 'condition')
cell_type_objects <- add_inflammation_status_each_object(cell_type_objects, condition_assignments, seurat_participant_column='best_match_sample', sheet_participants_column = 'sample', seurat_inflammation_column = 'inflammation_prev', sheet_inflammation_column = 'cond_prev')
# rename CA in the original assignment, and set a consensus one
for (cell_type in names(cell_type_objects)) {
  cell_type_object <- cell_type_objects[[cell_type]]
  cell_type_object@meta.data[!is.na(cell_type_object@meta.data[['inflammation_sheet']]) &
                               cell_type_object@meta.data[['inflammation_sheet']] == '24hCa', 'inflammation_sheet'] <- '24hCA'
  # now also set the final inflammation assignment
  cell_type_object@meta.data[['inflammation_final']] <- cell_type_object@meta.data[['inflammation_sheet']]
  cell_type_object@meta.data[is.na(cell_type_object@meta.data[['inflammation_final']]), 'inflammation_final'] <- cell_type_object@meta.data[is.na(cell_type_object@meta.data[['inflammation_final']]), 'inflammation_prev']
  # add back to list
  cell_type_objects[[cell_type]] <- cell_type_object
}

# add pool column
lane_remapping <- combine_lanes(unique(cell_type_objects[['monocyte']]@meta.data$lane))
# add to the object
cell_type_objects[['monocyte']]@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[cell_type_objects[['monocyte']]@meta.data[['lane']]]))
cell_type_objects[['monocyte']]@meta.data[['cell_type']] <- 'monocyte'
cell_type_objects[['monocyte']]@meta.data[['sample_lane']] <- paste(cell_type_objects[['monocyte']]@meta.data[['sample_final']], cell_type_objects[['monocyte']]@meta.data[['lane']], sep = ';;')
# get the table
acc_mono <- get_ncell_expressed_matrix(cell_type_objects[['monocyte']], 'sample_lane', assay = 'peaks')
  # write the result
write.table(acc_mono, gzfile(paste(out_tables_folder, '/', 'monocyte', '.tsv.gz', sep = '')), row.names = F, col.names = T, sep = '\t', quote = F)