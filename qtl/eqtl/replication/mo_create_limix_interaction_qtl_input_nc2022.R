#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Marc-Jan Bonder
# Name: mo_create_limix_interaction_qtl_input_nc2022.R
# Function: create the limix-QTL compatible input files from the Seurat object
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(Seurat)
library(matrixStats)
library(textTinyR) # NOT IN CONTAINER
library(pbapply)
library(doParallel)
library(ggplot2)
library(roycols)
library(mdfiver)


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
create_metadata_tables <- function(seurat_object_metadata, psam, participant_column='donor_final', pool_column='lane', condition_column='inflammation_status', celltype_column='cell_type_safe', join_pools=T) {
  # we need a separate psam for each cell type
  psam_per_celltype <- list()
  # add the sequencing platform to the psam
  psam[['sequencing_platform']] <- 'BGIseq'
  # some variables are more complicated, as they are a combination of multiple variables
  psam[['sequencing_run']] <- NA # lane
  psam[['sequencing_lane']] <- NA # lane
  psam[['scrna_platform']] <- '10x_multiome_gex'
  psam[['plate_base']] <- 'N'
  psam[['umi_based']] <- 'Y'
  psam[['biomaterial']] <- 'PBMC'
  psam[['sorting']] <- 'N'
  psam[['cell_treatment']] <- NA
  psam[['sample_condition']] <- 'control' # inflamed/non-inflamed
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
    psam[row_i, 'cell_treatment'] <- paste(conditions, collapse = ',')
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
        new_row[, 'cell_treatment'] <- condition
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
#' @param condition_column the seurat metadata column that denotes the condition of the sample (inflamed, non-inflamed)
#' @param min_cell_number the minimal number of cells to need to build a pseudobulk, pseudobulks with less cells are removed
#' @param min_numi the minimal number of UMIs to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param verbose print progress or not
#' @param single_thread do all work on a single thread (for debugging purposes)
#' @returns a list per cell type, each cell type has a list with the raw pseudobulk expression, the filtered pseudobulk expression, and the pcs
#' expression_per_celltype <- create_aggregated_expression_matrices(seurat_object, participant_column = 'soup_final_sample_assignment')
create_aggregated_expression_matrices <- function(seurat_object, participant_column='donor_final', celltype_column='cell_type_safe', batch_column=NULL, condition_column='inflammation_final', min_cell_number=5, min_numi=200, npcs=10, sample_cor_column='best_match_correlation', min_sample_cor=0, verbose=T, single_thread=F) {
  # subset object if min_numi parameter is given
  if (!is.null(min_numi) & !is.na(min_numi) & min_numi > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nFeature_RNA']] >= min_numi]
  }
  # if we have a minimal sample correlation metric, we will use that as well
  if (!is.null(min_sample_cor) & min_sample_cor > 0) {
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
    
    # do pseudobulking for the non-normalized counts
    aggregate_raw_count_matrix <- as.data.frame(
      pblapply(
        # go through each ID
        unique_id_list, FUN = function(x){
          # get the sparse means over the cells of a participant
          sparse_Means(count_matrix[, ids == x, drop = FALSE], rowMeans = TRUE)
        }
      )
    )
    # set the colnames to be the participants
    colnames(aggregate_raw_count_matrix) <- unique_id_list
    # and the genes as the rows
    rownames(aggregate_raw_count_matrix) <- rownames(count_matrix)
    
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
          message('doing mean expression normalization')
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
        # get the interaction term status for each sample
        part_to_interaction <- unique(seurat_object@meta.data[, c(participant_column, condition_column)])
        # first get the original rownames and colnames
        cov_rownames <- rownames(cov_out)
        cov_colnames <- colnames(cov_out)
        # get the conditions for the samples
        conditions <- part_to_interaction[match(rownames(cov_out), part_to_interaction[[participant_column]]), condition_column]
        # add the condition info
        cov_out <- cbind(conditions, cov_out)
        # set the dimension names again
        rownames(cov_out) <- cov_rownames
        colnames(cov_out) <- c(condition_column, cov_colnames)
        # put the results in a list
        aggregate_summary <- list('cell_type' = cell_type, 'expression' = aggregate_norm_count_matrix, 'expression_unfiltered' = aggregate_norm_count_matrix_unfiltered, 'expression_nonorm' = aggregate_raw_count_matrix, 'pc' = cov_out)
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
    qtl_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.qtlInput.txt.gz', sep = ''))
    pcs_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.qtlInput.Pcs.txt.gz', sep = ''))
    exp_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.Exp.txt.gz', sep = ''))
    metadata_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.covariates.txt.gz', sep = ''))
    #qtl_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.qtlInput.txt', sep = ''))
    #pcs_output_loc <- paste(output_loc, '/', cell_type, '.qtlInput.Pcs.txt', sep = '')
    #exp_output_loc <- paste(output_loc, '/', cell_type, '.Exp.txt', sep = '')
    #metadata_output_loc <- paste(output_loc, '/', cell_type, '.covariates.txt', sep = '')
    # write the files
    write.table(qtl_expression, qtl_output_loc, quote = F, sep = '\t', col.names = T)
    mdfiver::create_sha256_for_file(paste(output_loc, '/', cell_type, '.qtlInput.txt.gz', sep = ''))
    write.table(mean_expression, exp_output_loc, quote = F, sep = '\t', col.names = T)
    mdfiver::create_sha256_for_file(paste(output_loc, '/', cell_type, '.Exp.txt.gz', sep = ''))
    # in the case of pflogpf we also have a pseudobulk without pflogpf
    if ('expression_nonorm' %in% names(expression_per_celltype[[cell_type]])) {
      raw_output_loc <- gzfile(paste(output_loc, '/', cell_type, '.Exp.raw.txt.gz', sep = ''))
      write.table(expression_per_celltype[[cell_type]][['expression_nonorm']], raw_output_loc, quote = F, sep = '\t', col.names = NA)
      mdfiver::create_sha256_for_file(paste(output_loc, '/', cell_type, '.Exp.raw.txt.gz', sep = ''))
    }
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
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = T, row.names = F)
      mdfiver::create_sha256_for_file(paste(output_loc, '/', cell_type, '.covariates.txt.gz', sep = ''))
    }
    else {
      write.table(pcs, pcs_output_loc, quote = F, sep = '\t', col.names = T, row.names = T)
      mdfiver::create_sha256_for_file(paste(output_loc, '/', cell_type, '.qtlInput.Pcs.txt.gz', sep = ''))
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = T, row.names = F)
      mdfiver::create_sha256_for_file(paste(output_loc, '/', cell_type, '.covariates.txt.gz', sep = ''))
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
  write.table(smf, smf_output_loc, quote = F, sep = '\t', row.names = F)
  mdfiver::create_sha256_for_file(smf_output_loc)
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
#' @param min_numi the minimal number of UMIs to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param merge_pcs_into_covariates whether to merge the PCs into the covariates file
#' @param verbose print progress or not
#' @param single_thread do all work on a single thread (for debugging purposes)
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
                                    min_numi=200,
                                    npcs=10,
                                    sample_cor_column='soup_final_sample_correlation',
                                    min_sample_cor=0,
                                    merge_pcs_into_covariates=F, 
                                    verbose=T, 
                                    single_thread=F) {
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
    condition_column = condition_column, 
    min_cell_number = min_cell_number, 
    npcs = npcs, 
    sample_cor_column = sample_cor_column,
    min_sample_cor = min_sample_cor,
    min_numi = min_numi,
    verbose = verbose, 
    single_thread = single_thread
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


####################
# Settings        #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)

# size of chunks to normalize
chunk_size <- 5000
registerDoParallel(cores = 8)


####################
# Main Code        #
####################

# location of the v2 data
nc2022_v2_loc <- '/groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/demultiplexing/output/wg1_final_assignments_singlets_v2.rds'
# read the object
nc2022_v2 <- readRDS(nc2022_v2_loc)

# location of the condition annotations
condition_annotation_loc <- '/groups/umcg-franke-scrna/tmp02/releases/wijst-2020-hg19/v1/metadata/1M_timepoints_v2.tsv'
# and the original annotations for the donors
donor_annotation_loc <- '/groups/umcg-franke-scrna/tmp02/releases/wijst-2020-hg19/v1/metadata/1M_assignments_v2.tsv'
# read these two tables
condition_annotation <- read.table(condition_annotation_loc, header = T, sep = '\t')
donor_annotation <- read.table(donor_annotation_loc, header = T, sep = '\t')
# set the barcodes as the rownames
rownames(condition_annotation) <- condition_annotation[['barcode']]
rownames(donor_annotation) <- donor_annotation[['barcode']]

# add this information to the object
nc2022_v2 <- AddMetaData(nc2022_v2, donor_annotation['assignment'], 'assignment_previous')
nc2022_v2 <- AddMetaData(nc2022_v2, condition_annotation['timepoint'], 'condition_previous')

# check the combinations of original to new sample
nc2022_v2_old_donor_to_new <- data.frame(table(nc2022_v2@meta.data[, c('Assignment', 'assignment_previous', 'lane')]))
# rename the last column
colnames(nc2022_v2_old_donor_to_new)[4] <- 'overlapping'
# remove all the zeroes
nc2022_v2_old_donor_to_new <- nc2022_v2_old_donor_to_new[nc2022_v2_old_donor_to_new$overlapping > 0, ]
# get the totals for each donor as well
nc2022_v2_new_donor_totals <- data.frame(table(nc2022_v2@meta.data[, c('Assignment', 'lane')]))
# add that information
nc2022_v2_old_donor_to_new[['Assignment_total']] <- nc2022_v2_new_donor_totals[match(paste(nc2022_v2_old_donor_to_new[['Assignment']], nc2022_v2_old_donor_to_new[['lane']]), paste(nc2022_v2_new_donor_totals[['Assignment']], nc2022_v2_new_donor_totals[['lane']])), 'Freq']
# so we can calculate the fractions
nc2022_v2_old_donor_to_new[['Assignment_frac']] <- nc2022_v2_old_donor_to_new[['overlapping']] / nc2022_v2_old_donor_to_new[['Assignment_total']]

# quickly check to see how this looks at the donor level
nc2022_v2_old_donor_to_new[['Assignment_lane']] <- paste(nc2022_v2_old_donor_to_new[['Assignment']], nc2022_v2_old_donor_to_new[['lane']])
ggplot(data = nc2022_v2_old_donor_to_new, mapping = aes(x = Assignment_lane, y = Assignment_frac, fill = assignment_previous)) +
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('New donor/sample') +
  ylab('Fraction of old donor/sample') +
  ggtitle('Fraction overlap of sample/donor betwen NC2022 and sc-eQTLgen') +
  scale_fill_manual(values = roycols::get_color_list(nc2022_v2_old_donor_to_new[['assignment_previous']])) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  theme(legend.position="none") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# check the combinations of original to new sample
nc2022_v2_stim_to_new_donor <- data.frame(table(nc2022_v2@meta.data[, c('Assignment', 'condition_previous', 'lane')]))
# rename the last column
colnames(nc2022_v2_stim_to_new_donor)[4] <- 'overlapping'
# remove all the zeroes
nc2022_v2_stim_to_new_donor <- nc2022_v2_stim_to_new_donor[nc2022_v2_stim_to_new_donor$overlapping > 0, ]
# get the totals for each donor as well
nc2022_v2_new_donor_totals <- data.frame(table(nc2022_v2@meta.data[, c('Assignment', 'lane')]))
# add that information
nc2022_v2_stim_to_new_donor[['Assignment_total']] <- nc2022_v2_new_donor_totals[match(paste(nc2022_v2_stim_to_new_donor[['Assignment']], nc2022_v2_stim_to_new_donor[['lane']]), paste(nc2022_v2_new_donor_totals[['Assignment']], nc2022_v2_new_donor_totals[['lane']])), 'Freq']
# so we can calculate the fractions
nc2022_v2_stim_to_new_donor[['Assignment_frac']] <- nc2022_v2_stim_to_new_donor[['overlapping']] / nc2022_v2_stim_to_new_donor[['Assignment_total']]

# quickly check to see how this looks at the donor level
nc2022_v2_stim_to_new_donor[['Assignment_lane']] <- paste(nc2022_v2_stim_to_new_donor[['Assignment']], nc2022_v2_stim_to_new_donor[['lane']])
ggplot(data = nc2022_v2_stim_to_new_donor, mapping = aes(x = Assignment_lane, y = Assignment_frac, fill = condition_previous)) +
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('New donor/sample') +
  ylab('Fraction of old stim assignment') +
  ggtitle('Fraction overlap of old stim to new sc-eQTLgen donor assignment') +
  scale_fill_manual(values = roycols::get_color_list(nc2022_v2_stim_to_new_donor[['condition_previous']])) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  theme(legend.position="none") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# those look good enough, so now make a mapping
stim_assignment <- nc2022_v2_stim_to_new_donor[, c('Assignment', 'lane', 'condition_previous', 'Assignment_frac')]
# remove any NA
stim_assignment <- stim_assignment[!is.na(stim_assignment[['Assignment_frac']]), ]
# and sort by the fractions
stim_assignment <- stim_assignment[order(stim_assignment[['Assignment_frac']], decreasing = T), ]
# and keep only the best one for each sample-lane
stim_assignment <- stim_assignment[!duplicated(paste(stim_assignment[['Assignment']], stim_assignment[['lane']])), ]

# now add this information to the Seurat object
nc2022_v2@meta.data[, c('condition_previous', 'Assignment_frac')] <- stim_assignment[match(
  paste(nc2022_v2@meta.data[['Assignment']], nc2022_v2@meta.data[['lane']]), 
  paste(stim_assignment[['Assignment']], stim_assignment[['lane']]), 
), c('condition_previous', 'Assignment_frac')]


# location of the celltype annotation
celltype_annotation_loc <- '/groups/umcg-franke-scrna/tmp02/releases/wijst-2020-hg19/v1/metadata/1M_cell_types_v2.tsv'
celltype_annotation <- read.table(celltype_annotation_loc, header = T, sep = '\t')
rownames(celltype_annotation) <- celltype_annotation[['barcode']]

# add the celltype annotation to the Seurat object
nc2022_v2 <- AddMetaData(nc2022_v2, celltype_annotation['cell_type'])
# make anything that is NA, into unknown
nc2022_v2@meta.data[is.na(nc2022_v2@meta.data[['cell_type']]), 'cell_type'] <- 'unknown'
# make a mapping table
lowres_mapping <- list(
  "memory CD8T" = 'CD8T',
  "mono 2" = 'monocyte', 
  "NK" = 'NK', 
  "naive CD4T" = 'CD4T', 
  "mono 1" = 'monocyte', 
  "naive CD8T" = 'CD8T', 
  "th2 CD4T" = 'CD4T', 
  "NKdim" = 'NK', 
  "reg CD4T" = 'CD4T', 
  "pDC" = 'DC', 
  "unknown" = 'unknown', 
  "mDC" = 'DC', 
  "megakaryocyte" = 'megakaryocyte', 
  "plasma B" = 'plasmablast',  
  "mono 4" = 'monocyte', 
  "th1 CD4T" = 'CD4T', 
  "B" = 'B', 
  "mono2" = 'monocyte', 
  "hemapoietic stem" = 'hemapoietic_stem', 
  "NKbright" = 'NK', 
  "memory CD4T" = 'CD4T'
)
# add this lower resolution annotation
nc2022_v2@meta.data[['cell_type_lowerres']] <- as.vector(unlist(lowres_mapping[nc2022_v2@meta.data[['cell_type']]]))
# remove anywhere where we are unsure of the stim condition
nc2022_v2 <- nc2022_v2[, !is.na(nc2022_v2@meta.data[['condition_previous']]) & !is.na(nc2022_v2@meta.data[['Assignment_frac']]) & nc2022_v2@meta.data[['Assignment_frac']] > .7]

# make a mapping for the conditions
condition_mapping <- list(
  'UT' = 0,
  '3hCA' = 1,
  '3hMTB' = 1, 
  '3hPA' = 1,
  '24hCA' = 2, 
  '24hMTB' = 2,
  '24hPA' = 2
)
# add coded condition
nc2022_v2@meta.data[['condition_coded']] <- as.vector(unlist(condition_mapping[nc2022_v2@meta.data[['condition_previous']]]))

# get the psam location
psam_loc <- '/groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/unimputed/all/cytosnp_all_1m/cytosnp_all_1m.psam'
psam <- read.delim(psam_loc, as.is = T, check.names = F)
# write the matrices
do_limix_input_pipeline(
  seurat_object = nc2022_v2[, !is.na(nc2022_v2@meta.data[['condition_previous']]) & nc2022_v2@meta.data[['condition_previous']] %in% c('UT', '24hCA')], 
  psam = psam, 
  output_loc='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/nc2022_replication/24hCA/L1/',
  participant_column='Assignment',
  celltype_column='cell_type_lowerres',
  merge_pcs_into_covariates=F,
  join_pools=F,
  min_numi=200,
  min_sample_cor=0,
  condition_column='condition_coded',
  pool_column = 'lane',
  sample_cor_column=NULL
)
# write the matrices
do_limix_input_pipeline(
  seurat_object = nc2022_v2[, !is.na(nc2022_v2@meta.data[['condition_previous']]) & nc2022_v2@meta.data[['condition_previous']] %in% c('UT', '24hCA', '24hMTB', '24hPA')], 
  psam = psam, 
  output_loc='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/nc2022_replication/24h/L1/',
  participant_column='Assignment',
  celltype_column='cell_type_lowerres',
  merge_pcs_into_covariates=F,
  join_pools=F,
  min_numi=200,
  min_sample_cor=0,
  condition_column='condition_coded',
  pool_column = 'lane',
  sample_cor_column=NULL
)
# write the matrices
do_limix_input_pipeline(
  seurat_object = nc2022_v2[, !is.na(nc2022_v2@meta.data[['condition_previous']])], 
  psam = psam, 
  output_loc='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/nc2022_replication/timepoint/L1/',
  participant_column='Assignment',
  celltype_column='cell_type_lowerres',
  merge_pcs_into_covariates=F,
  join_pools=F,
  min_numi=200,
  min_sample_cor=0,
  condition_column='condition_coded',
  pool_column = 'lane',
  sample_cor_column=NULL
)
