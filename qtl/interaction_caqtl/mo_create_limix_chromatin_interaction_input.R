#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Marc-Jan Bonder
# Name: mo_create_limix_chromatin_input.R
# Function: create the limix-QTL compatible input files from the Seurat object
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(Seurat)
library(Signac)
library(matrixStats)
library(textTinyR) # NOT IN CONTAINER
library(pbapply)
library(qsmooth) # NOT IN CONTAINER
library(doParallel)


####################
# Functions        #
####################


#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object_metadata the metadata (from Seurat) to use to create a metadata annotation file
#' @param participant_column the seurat metadata column that denotes the participant
#' @param pool_column the seurat metadata column that denotes the pool that the samples were run in
#' @param condition_column the seurat metadata column that denotes the condition of the sample (inflamed, non-inflamed)
#' @param celltype_column rthe seurat metadata column that denotes the celltype of the cell
#' @returns a list per cell type, with for each cell type a metadata summary table
#' metadata_per_celltype <- create_metadata_tables(seurat_object@meta.data, donor_annotation_psam, participant_column = 'soup_final_sample_assignment')
create_metadata_tables <- function(seurat_object_metadata, psam, participant_column='donor_final', pool_column='day', condition_column='inflammation_status', celltype_column='cell_type_safe', join_pools=T) {
  # we need a separate psam for each cell type
  psam_per_celltype <- list()
  # add the sequencing platform to the psam
  psam[['sequencing_platform']] <- 'BGIseq'
  # some variables are more complicated, as they are a combination of multiple variables
  psam[['sequencing_run']] <- NA # lane
  psam[['sequencing_lane']] <- NA # lane
  psam[['scrna_platform']] <- '10x_v3.1'
  psam[['plate_base']] <- 'N'
  psam[['umi_based']] <- 'Y'
  psam[['biomaterial']] <- 'gut_biopsy'
  psam[['sorting']] <- 'N'
  psam[['cell_treatment']] <- 'UT'
  psam[['sample_condition']] <- NA # inflamed/non-inflamed
  psam[['Donor_Pool']] <- NA # participant;;lane
  # now check the values that vary between individuals
  for (row_i in 1:nrow(psam)) {
    # extract the participant
    participant <- psam[row_i, 'IID']
    # check which lanes this participant is in
    lanes <- unique(seurat_object_metadata[seurat_object_metadata[[participant_column]] == participant, pool_column])
    # check which conditions this sample was in
    conditions <- unique(seurat_object_metadata[seurat_object_metadata[[participant_column]] == participant, condition_column])
    # check if the lanes are not empty
    if (length(lanes) == 0) {
      # set as unknown if this is the case
      lanes <- c(NA)
    }
    # and for conditions
    if (length(conditions) == 0) {
      conditions <- c(NA)
    }
    # now set these values
    psam[row_i, 'sequencing_run'] <- paste(lanes, collapse=',')
    psam[row_i, 'sequencing_lane'] <- paste(lanes, collapse = ',')
    psam[row_i, 'sample_condition'] <- paste(conditions, collapse = ',')
    psam[row_i, 'Donor_Pool'] <- paste(participant, lanes, sep=';;', collapse = ',')
  }
  if (!join_pools) {
    # if we shouldn't join pools, let's get those specific ones
    indices_multiple_pools <- grep(',', psam[['Donor_Pool']])
    rows_multiple_pools <- psam[indices_multiple_pools, ]
    # and remove them from the original
    rows_single_pools <- psam[-indices_multiple_pools, ]
    # create a new list
    rows_split <- list()
    for (row_i in 1:nrow(rows_multiple_pools)){
      # extract the donor pools
      donor_pools <- strsplit(rows_multiple_pools[row_i, 'Donor_Pool'], ',')[[1]]
      # and the sample condition
      sample_conditions <- strsplit(rows_multiple_pools[row_i, 'sample_condition'], ',')[[1]]
      # and the sequencing lane
      sequencing_lanes <- strsplit(rows_multiple_pools[row_i, 'sequencing_lane'], ',')[[1]]
      # and the sequencing run
      sequencing_runs <- strsplit(rows_multiple_pools[row_i, 'sequencing_run'], ',')[[1]]
      # now add those as split column
      for (i in 1:length(donor_pools)) {
        # get individual pool and sample condition
        donor_pool <- donor_pools[i]
        condition <- sample_conditions[i]
        lane <- sequencing_lanes[i]
        run <- sequencing_runs[i]
        # create a new row
        new_row <- rows_multiple_pools[row_i, ]
        new_row[, 'Donor_Pool'] <- donor_pool
        new_row[, 'sample_condition'] <- condition
        new_row[, 'sequencing_lane'] <- lane
        new_row[, 'sequencing_run'] <- run
        # add it to the new split rows
        rows_split[[donor_pool]] <- new_row
      }
    }
    # merge the new rows
    df_rows_split <- do.call(rbind, rows_split)
    # then merge the old and the new
    psam <- rbind(rows_single_pools, df_rows_split)
  }
  # now check each cell type
  for (cell_type in unique(seurat_object_metadata[[celltype_column]])) {
    # subset the metadata to that celltype
    metadata_celltype <- seurat_object_metadata[seurat_object_metadata[[celltype_column]] == cell_type, ]
    # create a new psam for this cell type
    psam_celltype <- psam
    # if we joined the pools, we can just use the donor to get the number of cells
    if (join_pools) {
      # now use table to get the cell type numbers
      cell_numbers_per_donor <- data.frame(table(metadata_celltype[[participant_column]]))
      # now add the cell counts
      psam_celltype[['CellCount']] <- cell_numbers_per_donor[match(psam_celltype[['IID']], cell_numbers_per_donor[['Var1']]), 'Freq']
    }
    # if we have the donors entered multiple times with different pools, we need to get the cells for each specific combination
    else{
      # paste together the donor and the pool
      metadata_celltype[['Donor_Pool']] <- paste(metadata_celltype[[participant_column]], metadata_celltype[[pool_column]], sep = ';;')
      # and instead get the cell numbers like that
      cell_numbers_per_donor_pool <- data.frame(table(metadata_celltype[['Donor_Pool']]))
      # and add that
      psam_celltype[['CellCount']] <- cell_numbers_per_donor_pool[match(psam_celltype[['Donor_Pool']], cell_numbers_per_donor_pool[['Var1']]), 'Freq']
    }
    # turn NA into zero
    psam_celltype[is.na(psam_celltype[['CellCount']]), 'CellCount'] <- 0
    # add to list of psams
    psam_per_celltype[[cell_type]] <- psam_celltype
  }
  return(psam_per_celltype)
}


#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object the metadata (from Seurat) to use to create a metadata annotation file
#' @param participant_column the seurat metadata column that denotes the participant
#' @param celltype_column the seurat metadata column that denotes the celltype of the cell
#' @param batch_column the batch the sample was processed in (optional)
#' @param min_cell_number the minimal number of cells to need to build a pseudobulk, pseudobulks with less cells are removed
#' @param min_peaks the minimal number of peaks to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param verbose print progress or not
#' @returns a list per cell type, each cell type has a list with the raw pseudobulk expression, the filtered pseudobulk expression, and the pcs
#' expression_per_celltype <- create_aggregated_expression_matrices(seurat_object, participant_column = 'soup_final_sample_assignment')
create_aggregated_expression_matrices_quantilemethod <- function(seurat_object, participant_column='donor_final', celltype_column='cell_type_safe', condition_column='inflammation_final', batch_column=NULL, min_cell_number=5, min_peaks=200, npcs=10, sample_cor_column='best_match_correlation', min_sample_cor=0, verbose=T) {
  # subset object if min_peaks parameter is given
  if (!is.null(min_peaks) & !is.na(min_peaks) & min_peaks > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nCount_peaks']] >= min_peaks]
  }
  # if we have a minimal sample correlation metric, we will use that as well
  if (min_sample_cor > 0) {
    seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[[sample_cor_column]]) & seurat_object@meta.data[[sample_cor_column]] >= min_sample_cor]
  }
  # if the batch column is non-empty, we need to make a new column that is the join of the batch column and the participant ID
  if (!is.null(batch_column)) {
    seurat_object@meta.data[['participant_pool']] <- paste(seurat_object@meta.data[[participant_column]], seurat_object@meta.data[[batch_column]], sep = ';;')
    participant_column <- 'participant_pool'
  }
  # we'll store the aggregated count matrices per cell type
  aggregation_per_celltype <- list()
  # get the unique cell types
  cell_types <- unique(seurat_object@meta.data[[celltype_column]])
  # excluding NA of course
  cell_types <- cell_types[!is.na(cell_types)]
  # check each cell type
  for (cell_type in cell_types) {
    if (verbose) {
      message(paste('calculating for', cell_type))
    }
    # get indices of cells that are of this cell type
    indices_cell_type <- which(seurat_object@meta.data[[celltype_column]] == cell_type)
    # get the count matrix where we have the correct cell type
    seurat_object <- seurat_object[, indices_cell_type]
    # ignore genes that are never expressed
    seurat_object <-  seurat_object[which(rowSums(seurat_object@assays$peaks@counts) != 0), ]
    
    if (verbose) {
      message('calculating aggregated expression')
    }
    
    # back up the participant names
    parts <- unique(seurat_object@meta.data[[participant_column]])
    # Seurat will replace underscores with dashes, so we will do the same
    parts_dashes <- gsub('_', '-', parts)
    # now we'll create a mapping
    parts_mapping <- as.list(parts)
    names(parts_mapping) <- parts_dashes
    
    # now create the mean expression matrix
    aggregate_norm_count_matrix <- AggregateExpression(
      seurat_object,
      asssays = c('peaks'),
      group.by = participant_column,
      return.seurat = F
    )[['peaks']]
    # and change back the names by using that mapping we had
    colnames(aggregate_norm_count_matrix) <- as.vector(unlist(parts_mapping[colnames(aggregate_norm_count_matrix)]))
    
    # save the unfiltered mean expression
    aggregate_norm_count_matrix_unfiltered <- aggregate_norm_count_matrix
    
    # calculate the number of cells per donor
    cell_numbers_per_donor <- data.frame(table(seurat_object@meta.data[[participant_column]]))
    # set more descriptive column names
    colnames(cell_numbers_per_donor) <- c('participant', 'number')
    
    # filter on aggregates that are based on a certain number of cells
    samples_with_min_cells <- cell_numbers_per_donor[cell_numbers_per_donor[['number']] >= min_cell_number, 'participant']
    aggregate_norm_count_matrix = aggregate_norm_count_matrix[, which(colnames(aggregate_norm_count_matrix) %in% samples_with_min_cells), drop = F]
    # check if we have any cells left
    if (length(aggregate_norm_count_matrix) != 0 & ncol(aggregate_norm_count_matrix) > 0) {
      # remove genes without any variation
      row_var_info = rowVars(as.matrix(aggregate_norm_count_matrix))
      aggregate_norm_count_matrix = aggregate_norm_count_matrix[which(row_var_info != 0), , drop = F]
      # check if we have any genes left
      if (length(aggregate_norm_count_matrix) != 0 & nrow(aggregate_norm_count_matrix) > 0) {
        if (verbose) {
          message('doing mean expression normalization')
        }
        
        # get grouping factor
        group_factor <- seurat_object@meta.data[match(colnames(aggregate_norm_count_matrix), seurat_object@meta.data[[participant_column]]), condition_column]
        # drop where we have no groups
        aggregate_norm_count_matrix <- aggregate_norm_count_matrix[, !is.na(group_factor)]
        group_factor <- group_factor[!is.na(group_factor)]
        
        # then normalize
        aggregate_norm_count_matrix <- qsmooth(as.matrix(aggregate_norm_count_matrix), group_factor = group_factor)
        # and extract result
        aggregate_norm_count_matrix <- qsmoothData(aggregate_norm_count_matrix)
        
        if (verbose) {
          message('performing PCA')
        }
        # Do PCA, and select first x components.
        pc_out = prcomp(t(aggregate_norm_count_matrix))
        # check how many pcs we have
        npcs_present <- ncol(pc_out$x)
        # now subset to that number of pcs if we can
        cov_out <- NULL
        if (npcs <= npcs_present) {
          cov_out = pc_out$x[, 1:npcs]
        }
        else{
          # if we have less pcs we should warn
          warning(paste('requested', as.character(npcs), 'pcs, but only', as.character(npcs_present), 'are present for', cell_type))
          cov_out <- pc_out$x
        }
        if (verbose) {
          message(paste('finished', cell_type))
        }
        # put the results in a list
        aggregate_summary <- list('cell_type' = cell_type, 'expression' = aggregate_norm_count_matrix, 'expression_unfiltered' = aggregate_norm_count_matrix_unfiltered, 'pc' = cov_out)
        # which in turn is put into another list
        aggregation_per_celltype[[cell_type]] <- aggregate_summary
      }
      else {
        message(paste('no genes left after checking for variation between genes for', cell_type, 'skipping cell type'))
      }
    }
    else {
      message(paste('no cells left after min_cells filter for', cell_type, ', skipping cell type'))
    }
  }
  return(aggregation_per_celltype)
}

inverse_normalize <- function(norm_count_matrix, verbose = T) {
  # get the number of rows
  nrow_matrix <- nrow(norm_count_matrix)
  # get the number of columns
  ncol_matrix <- ncol(norm_count_matrix)
  # calculate how many rows a chunk would be
  nrow_chunk <- chunk_size / ncol_matrix
  nrow_chunk <- ceiling(nrow_chunk)
  # calculate how many chunks we need
  n_chunks <- nrow_matrix / nrow_chunk
  # round up because we need that last chunk
  n_chunks <- ceiling(n_chunks)
  if (verbose) {
    message(paste('work has been split into', as.character(n_chunks), 'chunks, of', as.character(nrow_chunk), 'rows each'))
  }
  # now do this per chunk
  res_per_chunk <- foreach(i = 1:n_chunks) %dopar% {
    # calculate the chunk start and stop
    chunk_start <- (i - 1) * nrow_chunk + 1
    chunk_stop <- i * nrow_chunk
    # check if we are exceeding the number of rows
    if (chunk_stop > nrow_matrix) {
      # because then we will stop there
      chunk_stop <- nrow_matrix
    }
    # extract these rows
    norm_count_matrix_chunk <- norm_count_matrix[chunk_start : chunk_stop, ]
    # do inverse normal transform per gene.
    for (r_i in 1:nrow(norm_count_matrix_chunk)) {
      norm_count_matrix_chunk[r_i, ] = qnorm((rank(norm_count_matrix_chunk[r_i, ], na.last = 'keep')-0.5) / sum(!is.na(norm_count_matrix_chunk[r_i, ])))
      if (verbose & r_i %% 10000 == 0) {
        message(paste('chunk', as.character(i), 'processed', as.character(r_i), 'rows'))
      }
    }
    return(norm_count_matrix_chunk)
  }
  # merge chuncks
  norm_count_matrix <- do.call('rbind', res_per_chunk) # dopar should keep order of input (not order of execution), and as such I should not have to manually set the order to be the same
  return(norm_count_matrix)
}


#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object the metadata (from Seurat) to use to create a metadata annotation file
#' @param participant_column the seurat metadata column that denotes the participant
#' @param celltype_column the seurat metadata column that denotes the celltype of the cell
#' @param batch_column the batch the sample was processed in (optional)
#' @param min_cell_number the minimal number of cells to need to build a pseudobulk, pseudobulks with less cells are removed
#' @param min_peaks the minimal number of peaks to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param verbose print progress or not
#' @returns a list per cell type, each cell type has a list with the raw pseudobulk expression, the filtered pseudobulk expression, and the pcs
create_aggregated_expression_matrices_rnamethod <- function(seurat_object, participant_column='donor_final', celltype_column='cell_type_safe', batch_column=NULL, min_cell_number=5, min_peaks=200, npcs=10, sample_cor_column='best_match_correlation', min_sample_cor=0, verbose=T, single_thread=F) {
  # subset object if min_peaks parameter is given
  if (!is.null(min_peaks) & !is.na(min_peaks) & min_peaks > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nCount_peaks']] >= min_peaks]
  }
  # if we have a minimal sample correlation metric, we will use that as well
  if (min_sample_cor > 0) {
    seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[[sample_cor_column]]) & seurat_object@meta.data[[sample_cor_column]] >= min_sample_cor]
  }
  # if the batch column is non-empty, we need to make a new column that is the join of the batch column and the participant ID
  if (!is.null(batch_column)) {
    seurat_object@meta.data[['participant_pool']] <- paste(seurat_object@meta.data[[participant_column]], seurat_object@meta.data[[batch_column]], sep = ';;')
    participant_column <- 'participant_pool'
  }
  # we'll store the aggregated count matrices per cell type
  aggregation_per_celltype <- list()
  # get the unique cell types
  cell_types <- unique(seurat_object@meta.data[[celltype_column]])
  # excluding NA of course
  cell_types <- cell_types[!is.na(cell_types)]
  # check each cell type
  for (cell_type in cell_types) {
    if (verbose) {
      message(paste('calculating for', cell_type))
    }
    # get indices of cells that are of this cell type
    indices_cell_type <- which(seurat_object@meta.data[[celltype_column]] == cell_type)
    # get the count matrix where we have the correct cell type
    count_matrix_full <- GetAssayData(seurat_object[, indices_cell_type], slot = "counts")
    # ignore genes that are never expressed
    count_matrix <-  count_matrix_full[which(rowSums(count_matrix_full) != 0), ]
    
    # also get the metadata for this cell type
    metadata <- seurat_object@meta.data[indices_cell_type, ]
    
    # extract IDs
    ids <- metadata[[participant_column]]
    unique_id_list <- unique(ids)
    
    # create new object to store the counts in
    norm_count_matrix <- count_matrix
    rm(count_matrix)
    gc()
    
    if (verbose) {
      message('doing single-cell normalization')
    }
    # do mean sample-sum normalization
    sample_sum_info = colSums(norm_count_matrix)
    mean_sample_sum = mean(sample_sum_info)
    sample_scale = sample_sum_info / mean_sample_sum
    
    # divide each column by sample_scale
    norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
    
    if (verbose) {
      message('calculating mean expression')
    }
    
    # now create the mean expression matrix
    aggregate_norm_count_matrix <- as.data.frame(
      pblapply(
        # go through each ID
        unique_id_list, FUN = function(x){
          # get the sparse means over the cells of a participant
          sparse_Means(norm_count_matrix[, ids == x, drop = FALSE], rowMeans = TRUE)
        }
      )
    )
    # set the colnames to be the participants
    colnames(aggregate_norm_count_matrix) <- unique_id_list
    # and the genes as the rows
    rownames(aggregate_norm_count_matrix) <- rownames(norm_count_matrix)
    
    # save the unfiltered mean expression
    aggregate_norm_count_matrix_unfiltered <- aggregate_norm_count_matrix
    
    # calculate the number of cells per donor
    cell_numbers_per_donor <- data.frame(table(metadata[[participant_column]]))
    # set more descriptive column names
    colnames(cell_numbers_per_donor) <- c('participant', 'number')
    
    # filter on aggregates that are based on a certain number of cells
    samples_with_min_cells <- cell_numbers_per_donor[cell_numbers_per_donor[['number']] >= min_cell_number, 'participant']
    aggregate_norm_count_matrix = aggregate_norm_count_matrix[, which(colnames(aggregate_norm_count_matrix) %in% samples_with_min_cells), drop = F]
    # check if we have any cells left
    if (length(aggregate_norm_count_matrix) != 0 & ncol(aggregate_norm_count_matrix) > 0) {
      # remove genes without any variation
      row_var_info = rowVars(as.matrix(aggregate_norm_count_matrix))
      aggregate_norm_count_matrix = aggregate_norm_count_matrix[which(row_var_info != 0), , drop = F]
      # check if we have any genes left
      if (length(aggregate_norm_count_matrix) != 0 & nrow(aggregate_norm_count_matrix) > 0) {
        if (verbose) {
          message(paste('doing mean expression normalization across', as.character(nrow(aggregate_norm_count_matrix)), 'rows'))
        }
        if (single_thread) {
          # do inverse normal transform per gene.
          for (r_i in 1:nrow(aggregate_norm_count_matrix)) {
            aggregate_norm_count_matrix[r_i, ] = qnorm((rank(aggregate_norm_count_matrix[r_i, ], na.last = 'keep')-0.5) / sum(!is.na(aggregate_norm_count_matrix[r_i, ])))
            if (verbose & r_i %% 10000 == 0) {
              message(paste('processed', as.character(r_i), 'rows'))
            }
          }
        }
        else {
          aggregate_norm_count_matrix <- inverse_normalize(aggregate_norm_count_matrix)
        }
        if (verbose) {
          message('performing PCA')
        }
        # Do PCA, and select first x components.
        pc_out = prcomp(t(aggregate_norm_count_matrix))
        # check how many pcs we have
        npcs_present <- ncol(pc_out$x)
        # now subset to that number of pcs if we can
        cov_out <- NULL
        if (npcs <= npcs_present) {
          cov_out = pc_out$x[, 1:npcs]
        }
        else{
          # if we have less pcs we should warn
          warning(paste('requested', as.character(npcs), 'pcs, but only', as.character(npcs_present), 'are present for', cell_type))
          cov_out <- pc_out$x
        }
        if (verbose) {
          message(paste('finished', cell_type))
        }
        # put the results in a list
        aggregate_summary <- list('cell_type' = cell_type, 'expression' = aggregate_norm_count_matrix, 'expression_unfiltered' = aggregate_norm_count_matrix_unfiltered, 'pc' = cov_out)
        # which in turn is put into another list
        aggregation_per_celltype[[cell_type]] <- aggregate_summary
      }
      else {
        message(paste('no genes left after checking for variation between genes for', cell_type, 'skipping cell type'))
      }
    }
    else {
      message(paste('no cells left after min_cells filter for', cell_type, ', skipping cell type'))
    }
  }
  return(aggregation_per_celltype)
}

#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object the metadata (from Seurat) to use to create a metadata annotation file
#' @param participant_column the seurat metadata column that denotes the participant
#' @param celltype_column the seurat metadata column that denotes the celltype of the cell
#' @param batch_column the batch the sample was processed in (optional)
#' @param min_cell_number the minimal number of cells to need to build a pseudobulk, pseudobulks with less cells are removed
#' @param min_peaks the minimal number of peaks to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param verbose print progress or not
#' @param quantile use the smooth quantile normalization methd
#' @returns a list per cell type, each cell type has a list with the raw pseudobulk expression, the filtered pseudobulk expression, and the pcs
#' expression_per_celltype <- create_aggregated_expression_matrices(seurat_object, participant_column = 'soup_final_sample_assignment')
create_aggregated_expression_matrices <- function(seurat_object, participant_column='donor_final', celltype_column='cell_type_safe', condition_column='inflammation_final', batch_column=NULL, min_cell_number=5, min_peaks=200, npcs=10, sample_cor_column='best_match_correlation', min_sample_cor=0, verbose=T, quantile=T) {
  if (quantile) {
    expression_per_celltype <- create_aggregated_expression_matrices_quantilemethod(
      seurat_object = seurat_object, 
      participant_column = participant_column, 
      celltype_column = celltype_column, 
      batch_column = batch_column,
      min_cell_number = min_cell_number, 
      npcs = npcs, 
      sample_cor_column = sample_cor_column,
      min_sample_cor = min_sample_cor,
      min_peaks = min_peaks,
      verbose = verbose
    )
  }
  else {
    expression_per_celltype <- create_aggregated_expression_matrices_rnamethod(
      seurat_object = seurat_object, 
      participant_column = participant_column, 
      celltype_column = celltype_column, 
      batch_column = batch_column,
      min_cell_number = min_cell_number, 
      npcs = npcs, 
      sample_cor_column = sample_cor_column,
      min_sample_cor = min_sample_cor,
      min_peaks = min_peaks,
      verbose = verbose
    )
  }
}


#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param expression_per_celltype list with celltypes as key, with the filtered expression under the expression key, unfiltered expression under the expression_unfiltered key, and pcs under the pc key
#' @param metadata_per_celltype metadata annation per cell type in a list, where the keys are the cell types
#' @param output_loc where to place the output files
#' @param merge_pcs_into_covariates whether to merge the PCs into the covariates file
#' @param rename_samples_to_samplepools rename the sample names 'TEST_81' to the sample+pool 'TEST_81;;190102_lane1'. Required if pools were joined when creating the expression matrices
#' @returns 0 if succesfull
#' write_limix_input(expression_per_celltype, metadata_per_celltype, '/groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/eqtl/sc-eqtlgen/input/elmentaite_adult_martin_immune/cell_type_safe/')
write_limix_input <- function(expression_per_celltype, metadata_per_celltype, output_loc='./', merge_pcs_into_covariates=F, rename_samples_to_samplepools=T) {
  # we can only do the data that we have expression and metadata for
  cell_types <- intersect(names(expression_per_celltype), names(metadata_per_celltype))
  # create the directory to place the files in if it does not exist yet
  dir.create(output_loc, recursive = T)
  # now check each cell type
  for (cell_type in cell_types) {
    # extract the data
    metadata <- metadata_per_celltype[[cell_type]]
    mean_expression <- expression_per_celltype[[cell_type]][['expression_unfiltered']]
    qtl_expression <- expression_per_celltype[[cell_type]][['expression']]
    pcs <- expression_per_celltype[[cell_type]][['pc']]
    # rename the samples with sample+pool, if not done already
    if (rename_samples_to_samplepools) {
      # subset the expression with entries we have metadata for
      mean_expression <- mean_expression[, colnames(mean_expression) %in% metadata[['IID']]]
      qtl_expression <- qtl_expression[, colnames(qtl_expression) %in% metadata[['IID']]]
      # replace the sample names with sample+lane
      colnames(mean_expression) <- metadata[
        match(colnames(mean_expression), metadata[['IID']]), 'Donor_Pool'
      ]
      colnames(qtl_expression) <- metadata[
        match(colnames(qtl_expression), metadata[['IID']]), 'Donor_Pool'
      ]
      rownames(pcs) <- metadata[
        match(rownames(pcs), metadata[['IID']]), 'Donor_Pool'
      ]
    }
    else{
      # subset the expression with entries we have metadata for
      mean_expression <- mean_expression[, colnames(mean_expression) %in% metadata[['Donor_Pool']]]
      qtl_expression <- qtl_expression[, colnames(qtl_expression) %in% metadata[['Donor_Pool']]]
    }
    # paste together the output names
    # qtl_output_loc <- paste(output_loc, '/', cell_type, '.qtlInput.txt', sep = '')
    # pcs_output_loc <- paste(output_loc, '/', cell_type, '.qtlInput.Pcs.txt', sep = '')
    # exp_output_loc <- paste(output_loc, '/', cell_type, '.Exp.txt', sep = '')
    # metadata_output_loc <- paste(output_loc, '/', cell_type, '.covariates.txt', sep = '')
    # gzip the files when pasting together the path
    qtl_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.qtlInput.txt.gz', sep = ''))
    pcs_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.qtlInput.Pcs.txt.gz', sep = ''))
    exp_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.Exp.txt.gz', sep = ''))
    metadata_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.covariates.txt.gz', sep = ''))
    
    # write the files
    write.table(qtl_expression, qtl_output_loc, quote = F, sep = '\t', col.names = NA)
    write.table(mean_expression, exp_output_loc, quote = F, sep = '\t', col.names = NA)
    # change X.FFID back to #FID
    colnames(metadata) <- gsub('X\\.FID', '#FID', colnames(metadata))
    # either write the PCs together or separate from the covariates
    if (merge_pcs_into_covariates) {
      # turn into dataframes
      pcs <- data.frame(pcs)
      metadata <- data.frame(metadata)
      # extract samples from PCs
      pc_samples <- rownames(pcs)
      # extract the original column names
      pc_columns <- colnames(pcs)
      # extract the covariate samples
      cov_samples <- as.character(metadata[['Donor_Pool']])
      # extract the original columns
      cov_columns <- colnames(metadata)
      # check which we have in both cases
      joint_samples <- intersect(pc_samples, cov_samples)
      # subset both
      pcs <- pcs[pc_samples, ]
      metadata[as.character(metadata[['Donor_Pool']]) %in% joint_samples, ]
      # add the donor pool to the pcs as a explicit column
      pcs[['Donor_Pool']] <- rownames(pcs)
      # finally join them
      metadata <- merge(metadata, pcs, by = 'Donor_Pool')
      # the reorder back the columns
      metadata <- metadata[, c(cov_columns, pc_columns)]
      # change X.FFID back to #FID
      colnames(metadata) <- gsub('X\\.FID', '#FID', colnames(metadata))
      # and write the result
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = NA)
    }
    else {
      write.table(pcs, pcs_output_loc, quote = F, sep = '\t', col.names = NA)
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = NA)
    }
  }
  # extract the first metadata
  metadata_first <- metadata_per_celltype[[1]]
  # get the unique combinations of the donor pool and the donor
  smf <-  unique(metadata_first[,c('IID', 'Donor_Pool')])
  colnames(smf) <- c('genotype_id', 'phenotype_id')
  # create the output file
  smf_output_loc <- paste(output_loc, 'smf.txt', sep = '')
  # write the result
  write.table(smf, smf_output_loc, quote=F, sep = '\t', row.names=F)
  return(0)
}

#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object_metadata the metadata (from Seurat) to use to create a metadata annotation file
#' @param output_loc where to place the output files
#' @param participant_column the seurat metadata column that denotes the participant
#' @param pool_column the seurat metadata column that denotes the pool that the samples were run in
#' @param condition_column the seurat metadata column that denotes the condition of the sample (inflamed, non-inflamed)
#' @param celltype_column the seurat metadata column that denotes the celltype of the cell
#' @param join_pools whether to join the pools
#' @param min_cell_number the minimal number of cells to need to build a pseudobulk, pseudobulks with less cells are removed
#' @param min_peaks the minimal number of UMIs to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param merge_pcs_into_covariates whether to merge the PCs into the covariates file
#' @param verbose print progress or not
#' @param quantile use the smooth quantile normalization methd
#' @returns 0 if succesfull
#' do_limix_input_pipeline(seurat_object, psam, output_loc='./', partipant_column='soup_final_sample_assignment')
do_limix_input_pipeline <- function(seurat_object, 
                                    psam, 
                                    output_loc='./',
                                    participant_column='donor_final', 
                                    pool_column='day', 
                                    condition_column='inflammation_status', 
                                    celltype_column='cell_type_safe', 
                                    join_pools=T,
                                    min_cell_number=5, 
                                    min_peaks=200,
                                    npcs=10,
                                    sample_cor_column='soup_final_sample_correlation',
                                    min_sample_cor=0,
                                    merge_pcs_into_covariates=F, 
                                    verbose=T,
                                    quantile=T) {
  if (verbose) {
    message('creating metadata files')
  }
  # create metadata tables
  metadata_per_celltype <- create_metadata_tables(
    seurat_object_metadata = seurat_object@meta.data, 
    psam = psam, 
    participant_column = participant_column, 
    pool_column = pool_column, 
    condition_column = condition_column, 
    celltype_column = celltype_column,
    join_pools = join_pools
  )
  if (verbose) {
    message('creating mean expression matrices')
  }
  # if we joined the pools, we need that as a parameter for the count matrices
  batch_column <- NULL
  # as well as that we don't need to rename the expression matrix columns then
  rename_samples_to_samplepools <- T
  if (!join_pools) {
    batch_column <- pool_column
    rename_samples_to_samplepools <- F
  }
  # create aggregated counts
  expression_per_celltype <- create_aggregated_expression_matrices(
    seurat_object = seurat_object, 
    participant_column = participant_column, 
    celltype_column = celltype_column, 
    batch_column = batch_column,
    min_cell_number = min_cell_number, 
    npcs = npcs, 
    sample_cor_column = sample_cor_column,
    min_sample_cor = min_sample_cor,
    min_peaks = min_peaks,
    verbose = verbose,
    quantile = quantile
  )
  if (verbose) {
    message('writing results')
  }
  # write the results
  write_limix_input(
    expression_per_celltype = expression_per_celltype,
    metadata_per_celltype = metadata_per_celltype,
    output_loc = output_loc,
    merge_pcs_into_covariates = merge_pcs_into_covariates,
    rename_samples_to_samplepools = rename_samples_to_samplepools
  )
  return(0)
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


split_samples_by_metadata <- function(base_input_dir, base_output_dir, metadata_split_column='sample_condition', metadata_sample_column='Donor_Pool') {
  # create a regular expression
  metadata_regex <- paste('*', '.covariates.txt.gz', sep = '')
  # list the files
  metadata_files <- list.files(base_input_dir, pattern = metadata_regex)
  # loop through each file
  for (metadata_file in metadata_files) {
    # extract the base name
    celltype <- gsub('.covariates.txt.gz', '', metadata_file)
    # now read the actual metadata file
    metadata <- read.table(paste(base_input_dir, metadata_file, sep = ''), header = T, check.names = F, comment.char = '')
    # now extract the different splits
    splits <- unique(metadata[[metadata_split_column]])
    # but remove empty entries
    splits <- splits[!is.na(splits)]
    # now read the qtlInput
    qtl_input <- read.table(paste(base_input_dir, '/', celltype, '.qtlInput.txt.gz', sep = ''), header = T, check.names = F, comment.char = '', row.names = 1)
    # and the PCs if present
    pcs <- NULL
    pcs_loc <- paste(paste(base_input_dir, '/', celltype, '.qtlInput.Pcs.txt.gz', sep = ''))
    if (file.exists(pcs_loc)) {
      pcs <- read.table(pcs_loc, header = T, check.names = F, comment.char = '', row.names = 1)
    }
    # now check each split
    for (split in splits) {
      # create the folder
      dir.create(paste(base_output_dir, '/', split, '/', sep = ''), recursive = T)
      # get the samples in the split
      samples_split <- metadata[!is.na(metadata[[metadata_split_column]]) & metadata[[metadata_split_column]] == split, metadata_sample_column]
      # now subset the qtl input
      qtl_input_split <- qtl_input[, colnames(qtl_input) %in% samples_split]
      # and write the result
      qtl_input_split_loc <- paste(base_output_dir, '/', split, '/', celltype, '.qtlInput.txt.gz', sep = '')
      write.table(qtl_input_split, gzfile(qtl_input_split_loc), quote = F, sep = '\t', col.names = NA)
      # subset the metadata as well
      metadata_split <- metadata[!is.na(metadata[[metadata_split_column]]) & metadata[[metadata_split_column]] == split, ]
      # and write it
      metadata_split_loc <- paste(base_output_dir, '/', split, '/', celltype, '.covariates.txt.gz', sep = '')
      write.table(metadata_split, gzfile(metadata_split_loc), quote = F, sep = '\t', col.names = NA)
      # finally the pcs if they were there
      if (!is.null(pcs)) {
        # subset to this split
        pcs_split <- pcs[rownames(pcs) %in% samples_split, ]
        # and write result
        pcs_split_loc <- paste(base_output_dir, '/', split, '/', celltype, '.qtlInput.Pcs.txt.gz', sep = '')
        write.table(pcs_split, pcs_split_loc, quote = F, sep = '\t', col.names = NA)
      }
    }
  }
  return(0)
}

merge_split_expression_matrices <- function(split_matrices_dir, output_base_dir, expression_file_prepend='', expression_file_append='.Exp.txt.gz', npcs=10, single_thread=F, verbose=T) {
  # list the directories
  split_dirs <- list.dirs(split_matrices_dir, recursive = F, full.names = F)
  # get just the ones that have the matrices
  split_dirs <- split_dirs[grepl('(\\d+)_(\\d+)', split_dirs)]
  # we'll save per condition
  matrices_per_condition <- list()
  # go through each matrix dir
  for (matrix_dir in split_dirs) {
    # paste the path together
    matrices_condition_dir <- paste(split_matrices_dir, '/', matrix_dir, '/', sep = '')
    # now list all the conditions
    conditions <- list.dirs(matrices_condition_dir, recursive = F, full.names = F)
    # check each condition
    for (condition in conditions) {
      # add to the list if not already there
      if (!(condition %in% names(matrices_per_condition))) {
        matrices_per_condition[[condition]] <- list()
      }
      # paste the path of that condition together
      matrices_celltypes_dir <- paste(matrices_condition_dir, '/', condition, '/', sep = '')
      # now list the celltypes in that directory
      celltype_files <- list.files(matrices_celltypes_dir, full.names = F, include.dirs = F)
      # now subset to make sure we only have the files we need
      matrices_celltypes_pattern <- paste(expression_file_prepend, '*', expression_file_append, sep = '')
      celltype_files <- celltype_files[grepl(matrices_celltypes_pattern, celltype_files)]
      # check each celltype file
      for (celltype_file in celltype_files) {
        # get just the celltype name
        celltype <- gsub(expression_file_prepend, '', celltype_file)
        celltype <- gsub(expression_file_append, '', celltype)
        # add the celltype in the list if it is not already there
        if (!(celltype %in% names(matrices_per_condition[[condition]]))) {
          matrices_per_condition[[condition]][[celltype]] <- list()
        }
        # read the file
        expression_input <- read.table(paste(matrices_celltypes_dir, '/', celltype_file, sep = ''), header = T, sep = '\t', check.names = F, comment.char = '', row.names = 1)
        # add to the list
        matrices_per_condition[[condition]][[celltype]][[matrix_dir]] <- expression_input
      }
    }
  }
  # check each condition
  for (condition in names(matrices_per_condition)) {
    # check each celltype
    for (celltype in names(matrices_per_condition[[condition]])) {
      # now merge all of them
      merged_matrices_celltype <- NULL
      for (matrix_name in names(matrices_per_condition[[condition]][[celltype]])) {
        # get the matrix
        matrix_celltype <- matrices_per_condition[[condition]][[celltype]][[matrix_name]]
        # if it was the first part, it will be that one
        if (is.null(merged_matrices_celltype)) {
          merged_matrices_celltype <- matrix_celltype
        }
        # otherwise we need to merge
        else {
          # get sample names in both
          samples_both <- intersect(colnames(merged_matrices_celltype), colnames(matrix_celltype))
          # and merge them
          merged_matrices_celltype <- rbind(merged_matrices_celltype[, samples_both], matrix_celltype[, samples_both])
        }
      }
      # check if we have any cells left
      if (length(merged_matrices_celltype) != 0 & ncol(merged_matrices_celltype) > 0) {
        # backup original matrix
        aggregate_norm_count_matrix <- merged_matrices_celltype
        # remove genes without any variation
        row_var_info = rowVars(as.matrix(aggregate_norm_count_matrix))
        aggregate_norm_count_matrix = aggregate_norm_count_matrix[which(row_var_info != 0), , drop = F]
        # check if we have any genes left
        if (length(aggregate_norm_count_matrix) != 0 & nrow(aggregate_norm_count_matrix) > 0) {
          if (verbose) {
            message(paste('doing mean expression normalization across', as.character(nrow(aggregate_norm_count_matrix)), 'rows'))
          }
          if (single_thread) {
            # do inverse normal transform per gene.
            for (r_i in 1:nrow(aggregate_norm_count_matrix)) {
              aggregate_norm_count_matrix[r_i, ] = qnorm((rank(aggregate_norm_count_matrix[r_i, ], na.last = 'keep')-0.5) / sum(!is.na(aggregate_norm_count_matrix[r_i, ])))
              if (verbose & r_i %% 10000 == 0) {
                message(paste('processed', as.character(r_i), 'rows'))
              }
            }
          }
          else {
            aggregate_norm_count_matrix <- inverse_normalize(aggregate_norm_count_matrix)
          }
          if (verbose) {
            message('performing PCA')
          }
          # Do PCA, and select first x components.
          pc_out = prcomp(t(aggregate_norm_count_matrix))
          # check how many pcs we have
          npcs_present <- ncol(pc_out$x)
          # now subset to that number of pcs if we can
          cov_out <- NULL
          if (npcs <= npcs_present) {
            cov_out = pc_out$x[, 1:npcs]
          }
          else{
            # if we have less pcs we should warn
            warning(paste('requested', as.character(npcs), 'pcs, but only', as.character(npcs_present), 'are present for', celltype))
            cov_out <- pc_out$x
          }
          if (verbose) {
            message(paste('finished', celltype))
          }
          # make the output directories
          qtl_output_loc <- gzfile(paste(output_base_dir, '/', condition,  '/', celltype, '.qtlInput.txt.gz', sep = ''))
          pcs_output_loc <- gzfile(paste(output_base_dir, '/', condition, '/', celltype, '.qtlInput.Pcs.txt.gz', sep = ''))
          exp_output_loc <- gzfile(paste(output_base_dir, '/', condition, '/', celltype, '.Exp.txt.gz', sep = ''))

          # write the files
          write.table(aggregate_norm_count_matrix, qtl_output_loc, quote = F, sep = '\t', col.names = NA)
          write.table(merged_matrices_celltype, exp_output_loc, quote = F, sep = '\t', col.names = NA)
          write.table(cov_out, pcs_output_loc, quote = F, sep = '\t', col.names = NA)
        }
        else {
          message(paste('no genes left after checking for variation between genes for', celltype, 'skipping cell type'))
        }
      }
      else {
        message(paste('no cells left after min_cells filter for', celltype, ', skipping cell type'))
      }
    }
  }
}

####################
# Main Code        #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 500 * 1000 * 1024^2)

# set seed
set.seed(7777)

# size of chunks to normalize
chunk_size <- 5000000
registerDoParallel(cores = 8)

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

# donor annotation psam
donor_annotation_psam_batch1_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/processed/genotype/GSA2023_1044_025_V3/unimputed/GSA2022_1044_025_V3.psam'
donor_annotation_psam_batch1 <- read.delim(donor_annotation_psam_batch1_loc, as.is = T, check.names = F)
donor_annotation_psam_batch2_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/processed/genotype/GSA2023_1009/unimputed/mo_individuals.psam'
donor_annotation_psam_batch2 <- read.delim(donor_annotation_psam_batch2_loc, as.is = T, check.names = F)
donor_annotation_psam_batch3_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/processed/genotype/ugli/unimputed/chr_all.psam'
donor_annotation_psam_batch3 <- read.delim(donor_annotation_psam_batch3_loc, as.is = T, check.names = F)
# merge all
donor_annotation_psam <- rbind(donor_annotation_psam_batch1, donor_annotation_psam_batch2)
donor_annotation_psam <- rbind(donor_annotation_psam, donor_annotation_psam_batch3)


# add pool column
lane_remapping <- combine_lanes(unique(cell_type_objects[['monocyte']]@meta.data$lane))
# add to the object
cell_type_objects[['monocyte']]@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[cell_type_objects[['monocyte']]@meta.data[['lane']]]))
cell_type_objects[['monocyte']]@meta.data[['cell_type']] <- 'monocyte'
# create input matrices
do_limix_input_pipeline(seurat_object = cell_type_objects[['monocyte']], 
                        psam = donor_annotation_psam, 
                        output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/UT/',
                        participant_column='best_match_sample', 
                        pool_column='lane', 
                        condition_column='inflammation_final', 
                        celltype_column='cell_type',
                        join_pools=F,
                        min_cell_number=5, 
                        min_peaks=200,
                        npcs=10,
                        sample_cor_column='best_match_correlation', 
                        min_sample_cor=0,
                        merge_pcs_into_covariates=F, 
                        verbose=T,
                        quantile=F)


# add to the object
cell_type_objects[['NK']]@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[cell_type_objects[['NK']]@meta.data[['lane']]]))
cell_type_objects[['NK']]@meta.data[['cell_type']] <- 'NK'
# create input matrices
do_limix_input_pipeline(seurat_object = cell_type_objects[['NK']], 
                        psam = donor_annotation_psam, 
                        output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/UT/',
                        participant_column='best_match_sample', 
                        pool_column='lane', 
                        condition_column='inflammation_final', 
                        celltype_column='cell_type',
                        join_pools=F,
                        min_cell_number=5, 
                        min_peaks=200,
                        npcs=10,
                        sample_cor_column='best_match_correlation', 
                        min_sample_cor=0,
                        merge_pcs_into_covariates=F, 
                        verbose=T,
                        quantile=F)
# CD4T bow
cell_type_objects[['CD4T']]@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[cell_type_objects[['CD4T']]@meta.data[['lane']]]))
cell_type_objects[['CD4T']]@meta.data[['cell_type']] <- 'CD4T'
do_limix_input_pipeline(seurat_object = cell_type_objects[['CD4T']], 
                        psam = donor_annotation_psam, 
                        output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/UT/',
                        participant_column='best_match_sample', 
                        pool_column='lane', 
                        condition_column='inflammation_final', 
                        celltype_column='cell_type',
                        join_pools=F,
                        min_cell_number=5, 
                        min_peaks=200,
                        npcs=10,
                        sample_cor_column='best_match_correlation', 
                        min_sample_cor=0,
                        merge_pcs_into_covariates=F, 
                        verbose=T,
                        quantile=F)
# B
cell_type_objects[['B']]@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[cell_type_objects[['B']]@meta.data[['lane']]]))
cell_type_objects[['B']]@meta.data[['cell_type']] <- 'B'
do_limix_input_pipeline(seurat_object = cell_type_objects[['B']][, cell_type_objects[['B']][['inflammation_final']] == 'UT'], 
                        psam = donor_annotation_psam, 
                        output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/UT/',
                        participant_column='best_match_sample', 
                        pool_column='lane', 
                        condition_column='inflammation_final', 
                        celltype_column='cell_type',
                        join_pools=F,
                        min_cell_number=5, 
                        min_peaks=200,
                        npcs=10,
                        sample_cor_column='best_match_correlation', 
                        min_sample_cor=0,
                        merge_pcs_into_covariates=F, 
                        verbose=T,
                        quantile=F)

# DC
cell_type_objects[['DC']]@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[cell_type_objects[['DC']]@meta.data[['lane']]]))
cell_type_objects[['DC']]@meta.data[['cell_type']] <- 'DC'
do_limix_input_pipeline(seurat_object = cell_type_objects[['DC']][, cell_type_objects[['DC']][['inflammation_final']] == 'UT'], 
                        psam = donor_annotation_psam, 
                        output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/UT/',
                        participant_column='best_match_sample', 
                        pool_column='lane', 
                        condition_column='inflammation_final', 
                        celltype_column='cell_type',
                        join_pools=F,
                        min_cell_number=5, 
                        min_peaks=200,
                        npcs=10,
                        sample_cor_column='best_match_correlation', 
                        min_sample_cor=0,
                        merge_pcs_into_covariates=F, 
                        verbose=T,
                        quantile=F)
do_limix_input_pipeline(seurat_object = cell_type_objects[['DC']][, cell_type_objects[['DC']][['inflammation_final']] == '24hCA'], 
                        psam = donor_annotation_psam, 
                        output_loc='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/24hCA/',
                        participant_column='best_match_sample', 
                        pool_column='lane', 
                        condition_column='inflammation_final', 
                        celltype_column='cell_type',
                        join_pools=F,
                        min_cell_number=5, 
                        min_peaks=200,
                        npcs=10,
                        sample_cor_column='best_match_correlation', 
                        min_sample_cor=0,
                        merge_pcs_into_covariates=F, 
                        verbose=T,
                        quantile=F)

