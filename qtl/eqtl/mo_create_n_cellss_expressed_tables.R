#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_n_cellss_expressed_tables.R
# Function: create tables that have the number of cells per sample that express each gene
############################################################################################################################

####################
# libraries        #
####################

# for reading Seurat object
library(Seurat)
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


write_ncell_expressed_matrices <- function(seurat_objects_folder, output_table_folder, seurat_object_prepend='mo_all_20240619_seuratv5_annotated_agesexcovid_', seurat_object_append='.rds', cell_types=NULL, sample_column='sample_lane', layer='counts', assay='RNA', multithread=T, nthreads=4) {
  # get the files to do
  files_to_process <- NULL
  # and the cell types they belong to
  file_celltypes <- NULL
  # either get them using the supplied cell types
  if (!is.null(cell_types)) {
    files_to_process <- paste(seurat_objects_folder, seurat_object_prepend, cell_types, seurat_object_append, sep = '')
    # and those will be the cell types
    file_celltypes <- cell_types
  }
  # or do based on the names and a wildcard
  else {
    # all files
    files_to_process <- list.files(seurat_objects_folder)
    # then make a regex
    files_regex <- paste('^', seurat_object_prepend, '.*', seurat_object_append, '$', sep = '')
    # and subset
    files_to_process <- files_to_process[grepl(files_regex, files_to_process)]
    # then get the celltypes from those files
    file_celltypes <- gsub(paste('^', seurat_object_prepend, sep = ''), '', files_to_process)
    file_celltypes <- gsub(paste(seurat_object_append, '$', sep = ''), '', file_celltypes)
  }
  # do each of the files
  for (i in 1:length(files_to_process)) {
    # get the file
    file_to_process <- files_to_process[i]
    # and the cell type
    celltype_to_process <- file_celltypes[i]
    # paste together the full paths
    file_to_process_full <- paste(seurat_objects_folder, '/', file_to_process, sep = '')
    output_table_loc <- paste(output_table_folder, '/', celltype_to_process, '.tsv.gz', sep = '')
    # read the Seurat object
    seurat_object <- readRDS(file_to_process_full)
    # get the results for this cell type
    celltype_result <- get_ncell_expressed_matrix(seurat_object, sample_column = sample_column, layer = layer, assay = assy, multithread = multithread, nthreads = nthreads)
    # write the result to a table
    write.table(celltype_result, gzfile(output_table_loc), row.names = F, col.names = T, sep = '\t', quote = F)
  }
  return(0)
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

# locations of objects
objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
mo_object_loc <- paste(objects_loc, 'mo_all_20240223_seuratv5_normalized.rds', sep = '')
# the seurat object prepend
seurat_object_prepend <- 'mo_all_20240619_seuratv5_annotated_agesexcovid_'
seurat_object_append <- '.rds'
# where to place the outputs
out_tables_folder <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/ncells_by_genes/L1/combined/'

# check the cell types
for (cell_type in c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')) {
  # read the file
  seurat_object_ct <- readRDS(paste(objects_loc, '/', seurat_object_prepend, cell_type, seurat_object_append, sep = ''))
  # add sample + lane
  seurat_object_ct@meta.data[['sample_lane']] <- paste(seurat_object_ct@meta.data[['sample_final']], seurat_object_ct@meta.data[['lane']], sep = ';;')
  # get the table
  expr_celltype <- get_ncell_expressed_matrix(seurat_object_ct, 'sample_lane')
  # write the result
  write.table(expr_celltype, gzfile(paste(out_tables_folder, '/', cell_type, '.tsv.gz', sep = '')), row.names = F, col.names = T, sep = '\t', quote = F)
}
