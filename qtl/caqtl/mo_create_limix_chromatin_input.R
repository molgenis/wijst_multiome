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
#' @param min_numi the minimal number of UMIs to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param verbose print progress or not
#' @returns a list per cell type, each cell type has a list with the raw pseudobulk expression, the filtered pseudobulk expression, and the pcs
#' expression_per_celltype <- create_aggregated_expression_matrices(seurat_object, participant_column = 'soup_final_sample_assignment')
create_aggregated_expression_matrices <- function(seurat_object, participant_column='donor_final', celltype_column='cell_type_safe', condition_column='inflammation_final', batch_column=NULL, min_cell_number=5, min_numi=200, npcs=10, sample_cor_column='best_match_correlation', min_sample_cor=0, verbose=T) {
  # subset object if min_numi parameter is given
  if (!is.null(min_numi) & !is.na(min_numi) & min_numi > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nCount_peaks']] >= min_numi]
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
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = T, row.names = F)
    }
    else {
      write.table(pcs, pcs_output_loc, quote = F, sep = '\t', col.names = NA, row.names = F)
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = NA, row.names = F)
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
#' @param min_numi the minimal number of UMIs to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param sample_cor_column the column with the sample correlations
#' @param min_sample_cor the minimal sample correlation to keep a cell (leave at zero to do no filtering)
#' @param merge_pcs_into_covariates whether to merge the PCs into the covariates file
#' @param verbose print progress or not
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
                                    verbose=T) {
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
    min_numi = min_numi,
    verbose = verbose
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
# Main Code        #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 500 * 1000 * 1024^2)

# set seed
set.seed(7777)

# location of the condition assignment
condition_assignment_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_monocyte_based_condition_numbers.tsv'
# read the conditions
condition_assignments <- read.table(condition_assignment_loc, header = T, sep = '\t')

# location of the cell type objects
cell_type_objects_wstatus_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltypemajor_wstatus_1_64.rds'

# read the object
cell_type_objects <- readRDS(cell_type_objects_wstatus_loc)

# donor annotation psam
donor_annotation_psam_batch1_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/genotype/GSA2023_1044_025_V3/unimputed/GSA2022_1044_025_V3.psam'
donor_annotation_psam_batch1 <- read.delim(donor_annotation_psam_batch1_loc, as.is = T, check.names = F)
donor_annotation_psam_batch2_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/genotype/GSA2023_1009/unimputed/mo_individuals.psam'
donor_annotation_psam_batch2 <- read.delim(donor_annotation_psam_batch2_loc, as.is = T, check.names = F)
donor_annotation_psam_batch3_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/genotype/ugli/unimputed/chr_all.psam'
donor_annotation_psam_batch3 <- read.delim(donor_annotation_psam_batch3_loc, as.is = T, check.names = F)
# merge all
donor_annotation_psam <- rbind(donor_annotation_psam_batch1, donor_annotation_psam_batch2)
donor_annotation_psam <- rbind(donor_annotation_psam, donor_annotation_psam_batch3)


# add pool column
lane_remapping <- combine_lanes(unique(cell_type_objects[['monocyte']]@meta.data$lane))
# add to the object
cell_type_objects[['monocyte']]@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[cell_type_objects[['monocyte']]@meta.data[['lane']]]))
cell_type_objects[['monocyte']]@meta.data[['cell_type']] <- 'monocyte'

create_aggregated_expression_matrices(seurat_object = cell_type_objects[['monocyte']], 
                                                  participant_column='best_match_sample', 
                                                  celltype_column='cell_type', 
                                                  condition_column='inflammation_final', 
                                                  batch_column='lane_both', 
                                                  min_cell_number=5, 
                                                  min_numi=200, 
                                                  npcs=10, 
                                                  sample_cor_column='best_match_correlation', 
                                                  min_sample_cor=0, 
                                                  verbose=T)
  

do_limix_input_pipeline(seurat_object = cell_type_objects[['monocyte']], 
                                    psam = donor_annotation_psam, 
                                    output_loc='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/',
                                    participant_column='best_match_sample', 
                                    pool_column='lane_both', 
                                    condition_column='inflammation_final', 
                                    celltype_column='cell_type',
                                    join_pools=F,
                                    min_cell_number=5, 
                                    min_numi=200,
                                    npcs=10,
                                    sample_cor_column='best_match_correlation', 
                                    min_sample_cor=0,
                                    merge_pcs_into_covariates=T, 
                                    verbose=T)
