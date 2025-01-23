#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_deconstruct_signac_objects.R
# Function: deconstruct signac object to import into pycistopic
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(mdfiver) # not in container yet
library(Matrix)
library(spam)
library(spam64) # not in container yet

####################
# Functions        #
####################

export_atac_counts_object <- function(signac_object, out_folder) {
  # features
  features_gz <- gzfile(paste(out_folder, 'features.tsv.gz', sep = ''))
  write.table(data.frame(x = rownames(signac_object@assays$peaks@counts)), features_gz, row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(paste(out_folder, 'features.tsv.gz', sep = ''))
  # barcodes
  barcodes_gz <- gzfile(paste(out_folder, 'barcodes.tsv.gz', sep = ''))
  write.table(data.frame(x = colnames(signac_object@assays$peaks@counts)), barcodes_gz, row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(paste(out_folder, 'barcodes.tsv.gz', sep = ''))
  # metadata
  metadata_gz <- gzfile(paste(out_folder, 'metadata.tsv.gz', sep = ''))
  write.table(cbind(data.frame(bc = rownames(signac_object@meta.data)), signac_object@meta.data), metadata_gz, row.names = F, col.names = T, quote = F, sep = '\t')
  mdfiver::create_md5_for_file(paste(out_folder, 'metadata.tsv.gz', sep = ''))
  # and finally the count matrix
  counts_gz <- paste(out_folder, 'matrix.mtx', sep = '')
  writeMM(signac_object@assays$peaks@counts, counts_gz)
  mdfiver::create_md5_for_file(paste(out_folder, 'matrix.mtx', sep = ''))
  # upon success
  return(0)
}


####################
# Settings        #
####################

options(spam.force64 = TRUE)    # forcing 64-bit structure


####################
# Debugging        #
####################


####################
# Main Code        #
####################

# location to put the deconstructed objects
deconstructed_folders_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/'

# location of the Signac objects
signac_objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/'
# prepend
signac_prepend <- 'mo_cpeaks_filtered_'
# append of each file
signac_append <- '_wstatus_1_80_20240709.rds'

# check each cell type
for (cell_type in c('b', 'cd4t', 'cd8t', 'dc', 'nk', 'monocyte')) {
  # full location to the object
  signac_full_loc <- paste(signac_objects_loc, '/', signac_prepend, cell_type, signac_append, sep = '')
  signac_object <- readRDS(signac_full_loc)
  # folder to the deconstructed output
  out_folder <- paste(deconstructed_folders_loc, '/', cell_type, '/', sep = '')
  # create the folder
  dir.create(out_folder)
  # do the run
  export_atac_counts_object(
    signac_object = signac_object, 
    out_folder = out_folder
  )
  # clear up the memory
  rm(signac_object)
}

# initialize a value in which we will store the combined atac matrix
atac_full_matrix <- NULL
atac_full_columns <- NULL
atac_full_features <- NULL
# check each cell type again
for (cell_type in c('b', 'cd4t', 'cd8t', 'dc', 'nk', 'monocyte')) {
  print(paste('reading', cell_type))
  # read accessibility
  accessibility <- spam::read.MM(paste(deconstructed_folders_loc, '/', cell_type, '/matrix.mtx', sep = ''))
  # read cells
  cell_names <- read.table(paste(deconstructed_folders_loc, '/', cell_type, '/barcodes.tsv.gz', sep = ''), header = F)$V1
  # read features
  features_names <- read.table(paste(deconstructed_folders_loc, '/', cell_type, '/features.tsv.gz', sep = ''), header = F)$V1
  # merge onto the matrix we already have
  if (is.null(atac_full_matrix)) {
    atac_full_matrix <- accessibility
    atac_full_columns <- cell_names
    atac_full_features <- features_names
  }else {
    # get features in either, so outer join
    all_features <- union(atac_full_features, atac_full_features)
    # order them, not required but just easier later on
    all_features <- all_features[order(all_features)]
    
    # get missing features in new matrix
    features_missing_accessibility <- setdiff(all_features, features_names)
    # make those as empty entries
    entries_missing_accessibility <- spam::spam(0, nrow = length(features_missing_accessibility), ncol = ncol(accessibility))
    # add those to the matrix
    #accessibility <- rbind(accessibility, entries_missing_accessibility)
    accessibility <- t(cbind(t(accessibility), t(entries_missing_accessibility)))
    # update the feature names
    features_names <- c(features_names, features_missing_accessibility)
    
    # get features missing in existing matrix
    features_missing_full <- setdiff(all_features, atac_full_features)
    # make those as empty entries
    entries_missing_full <- spam::spam(0, nrow = length(features_missing_full), ncol = ncol(atac_full_matrix))
    # add those to the matrix
    #atac_full_matrix <- rbind(atac_full_matrix, entries_missing_full)
    atac_full_matrix <- t(cbind(t(atac_full_matrix), t(entries_missing_full)))
    # update the feature names
    atac_full_features <- c(atac_full_features, features_missing_full)
    
    # order them to be the same
    atac_full_matrix <- atac_full_matrix[match(all_features, atac_full_features), ]
    accessibility <- accessibility[match(all_features, features_names), ]
    
    # and finally merge them
    atac_full_matrix <- cbind(atac_full_matrix, accessibility)
    # update the feature names
    atac_full_features <- all_features
    # and cell names
    atac_full_columns <- c(atac_full_columns, cell_names)
  }
}

# check each cell type for metadata
metadata_all_list <- list()
for (cell_type in c('b', 'cd4t', 'cd8t', 'dc', 'nk', 'monocyte')) {
  # read the metadata
  metadata_ct <- read.table(paste(deconstructed_folders_loc, '/', cell_type, '/metadata.tsv.gz', sep = ''), header = T, sep = '\t')
  # put in list
  metadata_all_list[[cell_type]] <- metadata_ct
}
# combine all metadata
metadata_all <- do.call('rbind', metadata_all_list)
# order it to be the same as the fragment matrix
metadata_all <- metadata_all[match(atac_full_columns, metadata_all[['bc']]), ]

# setup paths to output
merged_out_folder <- paste(deconstructed_folders_loc, '/', 'merged_major_celltypes/', sep = '')
dir.create(merged_out_folder)
merged_mtx_loc <- paste(merged_out_folder, '/matrix.mtx', sep = '')
merged_rds_loc <- paste(merged_out_folder, '/matrix.rds', sep = '')
merged_barcodes_loc <- paste(merged_out_folder, '/barcodes.tsv.gz', sep = '')
merged_features_loc <- paste(merged_out_folder, '/features.tsv.gz', sep = '')
merged_metadata_loc <- paste(merged_out_folder, '/metadata.tsv.gz', sep = '')
# write the results
write.table(data.frame(x = atac_full_features), gzfile(merged_features_loc), row.names = F, col.names = F, quote = F)
write.table(data.frame(x = atac_full_columns), gzfile(merged_barcodes_loc), row.names = F, col.names = F, quote = F)
write.table(metadata_all, gzfile(merged_metadata_loc), row.names = F, col.names = T, quote = F, sep = '\t')
saveRDS(atac_full_matrix, merged_rds_loc)
# generate checksums
mdfiver::create_md5_for_file(merged_rds_loc)
mdfiver::create_md5_for_file(merged_barcodes_loc)
mdfiver::create_md5_for_file(merged_features_loc)
mdfiver::create_md5_for_file(merged_metadata_loc)


# reload matrix
atac_full_matrix <- readRDS(merged_rds_loc)
# do a chunked save of the matrix
chunk_row_size <- 100000
# get how many chunks we need
chunks_needed <- ceiling(nrow(atac_full_matrix) / chunk_row_size)
# set 32 bit structure for conversion
options(spam.force64 = FALSE)
# check each chunk
for (chunk in 1 : chunks_needed) {
  # start of the chunk
  chunk_start <- ((chunk - 1) * chunk_row_size) + 1
  # end of the chunk
  chunk_end <- chunk * chunk_row_size
  # if the end is bigger than the total, we need to just take the end
  if (chunk_end > nrow(atac_full_matrix)) {
    chunk_end <- nrow(atac_full_matrix)
  }
  # subset the matrix
  atac_full_matrix_chunk <- atac_full_matrix[chunk_start : chunk_end, ]
  # covert to DGR
  atac_full_matrix_chunk_dgr <- spam::as.dgRMatrix.spam(atac_full_matrix_chunk)
  # save
  chunk_name <- chunk
  if (chunk_name < 10) {
    chunk_name <- paste('00', chunk_name, sep = '')
  }
  else if (chunk_name < 100) {
    chunk_name <- paste('0', chunk_name, sep = '')
  }
  Matrix::writeMM(atac_full_matrix_chunk_dgr, paste(merged_out_folder, 'matrix_chunk_', chunk_name, '.mtx', sep = ''))
  mdfiver::create_md5_for_file(paste(merged_out_folder, 'matrix_chunk_', chunk_name, '.mtx', sep = ''))
}
