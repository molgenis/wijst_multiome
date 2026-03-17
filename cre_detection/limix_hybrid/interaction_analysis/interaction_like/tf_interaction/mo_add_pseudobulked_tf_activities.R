#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_add_pseudobulked_tf_activities.R
# Function: add pseudobulked TF activity data to covariate matrices
############################################################################################################################

####################
# libraries        #
####################

library(Matrix)
library(data.table)
library(parallel)
library(foreach)
library(doParallel)

####################
# Functions        #
####################

read_mtx_barcodes_and_features <- function(parts_path, barcodes_file='eregulon_gene_auc_barcodes.txt.gz', features_file='eregulon_gene_auc_eregnames.txt.gz', mtx_file='eregulon_gene_auc.mtx.gz', barcodes_index=1, features_index=1) {
  # paste the filepaths together
  mtx_path <- paste(parts_path, mtx_file, sep = '/')
  features_path <- paste(parts_path, features_file, sep = '/')
  barcodes_path <- paste(parts_path, barcodes_file, sep = '/')
  # read features and barcodes first, they are the smallest
  barcodes <- read.table(barcodes_path, header = F, sep = '\t')[[barcodes_index]]
  features <- read.table(features_path, header = F, sep = '\t')[[features_index]]
  # now read the matrix file
  mtx <- Matrix::readMM(mtx_path)
  # set the column and row names
  colnames(mtx) <- barcodes
  rownames(mtx) <- features
  return(mtx)
}


####################
# Settings         #
####################


####################
# Main Code        #
####################

# location of the metadata
metadata_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
# location of the TF activity matrix
activity_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/'
# where to place the per-sample TF activity matrices
activity_ps_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/pseudobulked/L1/'
# columns to pseudobulk over
ps_columns <- c('sample_final', 'lane')
sample_column <- 'sample_final'
ct_column <- 'celltype_imputed_lowerres'
ps_columns_sep <- ';;'
# column to denote the cell in the metadata
cell_column='barcode_lane'

# get the TF matrix
activity <- read_mtx_barcodes_and_features(activity_loc)
# read the metadata
metadata <- fread(metadata_loc, header = T, sep = '\t')

# subset to barcodes present in both
barcodes_both <- intersect(metadata[[cell_column]], colnames(activity))
# and order both matrices in the same order
metadata <- metadata[match(barcodes_both, metadata[[cell_column]]), ]
activity <- activity[, barcodes_both]

# add the aggregating column
for (ps_column in c(ps_columns)) {
  # initialize column
  if (!('ps_column' %in% colnames(metadata))) {
    metadata[['ps_column']] <- metadata[[ps_column]]
  } else {
    # otherwise add to that column
    metadata[['ps_column']] <- paste(metadata[['ps_column']], metadata[[ps_column]], sep = ps_columns_sep)
  }
}
# make a list of things we keep track of
cols_to_keep <- unique(c('ps_column', ct_column, sample_column, ps_columns))
# store the ps column with the cell type column
ps_to_ct <- unique(metadata[, ..cols_to_keep])
# then add the ct column to the ps column
ps_to_ct[['ps_ct']] <- paste(ps_to_ct[['ps_column']], ps_to_ct[[ct_column]], sep = ps_columns_sep)
metadata[['ps_column']] <- paste(metadata[['ps_column']], metadata[[ct_column]], sep = ps_columns_sep)

# extract unique ps column entries
ps_column_entries <- unique(metadata[['ps_column']])

# we'll do this in portions so we can multithread this
chunk_size <- 10
# check how many chunks we need
chunks_required <- ceiling(nrow(activity) / chunk_size)
# check each chunk
activity_ps <- foreach(chunki = 1 : chunks_required, .combine=cbind) %dopar% {
# activity_ps_l <- list()
# activity_ps_l <- for(chunki in 1 : chunks_required) {
  # get the end index
  end_index <- chunki * chunk_size
  # get the starting index
  start_index <- end_index - chunk_size + 1
  # if the end index is larger than the number of rows, we were actually using the last chunk, which can be smaller
  if (end_index > nrow(activity)) {
    end_index <- nrow(activity)
  }
  # subset the matrix to these indices
  mtx_ss <- activity[start_index : end_index, ]
  # extract eregulon names
  ereg_names <- rownames(mtx_ss)
  # do the aggregation
  mtx_ss_agg <- aggregate(data.frame(as.matrix(t(mtx_ss))), by = list('ps' = metadata[['ps_column']]), FUN = mean)
  # order the same way as the ps column entries
  mtx_ss_agg <- mtx_ss_agg[match(mtx_ss_agg[['ps']], ps_column_entries), ]
  # drop the ps column
  mtx_ss_agg[['ps']] <- NULL
  # make sure the ereg names are set
  colnames(mtx_ss_agg) <- ereg_names
  return(mtx_ss_agg)
  # activity_ps_l[[chunki]] <- mtx_ss_agg
}
# activity_ps <- do.call('cbind', activity_ps_l)

# transpose because that is the format we expect
activity_ps_t <- t(activity_ps)
# set the samples as columns again
colnames(activity_ps_t) <- ps_column_entries

# check each cell type
for (cell_type in unique(ps_to_ct[[ct_column]])) {
  # get the entries for this cell type from the mapping
  mapping_ct <- ps_to_ct[ps_to_ct[[ct_column]] == cell_type, ]
  # keep only what we were able to aggregate
  mapping_ct_aggregated <- mapping_ct[mapping_ct[['ps_ct']] %in% colnames(activity_ps_t), ]
  # extract those entries from the bigger matrix
  activity_ps_t_ct <- activity_ps_t[, mapping_ct_aggregated[['ps_ct']]]
  # then replace the column names
  colnames(activity_ps_t_ct) <- mapping_ct_aggregated[['ps_column']]
  # extract the TF names
  tf_names <- rownames(activity_ps_t_ct)
  # extract the sample names
  sample_names <- colnames(activity_ps_t_ct)
  # convert to a dataframe
  activity_ps_t_ct <- data.frame(activity_ps_t_ct)
  # make sure we keep the right sample naems
  colnames(activity_ps_t_ct) <- sample_names
  # set the TF as first column
  activity_ps_t_ct <- cbind(data.frame('eregulon' = tf_names), activity_ps_t_ct)
  # make path to write
  ereg_dir <- paste(activity_ps_loc, cell_type, sep = '/')
  # create directory
  dir.create(ereg_dir, recursive = T, showWarnings = F)
  # make the path to the file
  ereg_file_loc <- paste(ereg_dir, 'eregulons.tsv.gz', sep = '/')
  # write the file
  write.table(activity_ps_t_ct, gzfile(ereg_file_loc), sep = '\t', quote = F, row.names = F, col.names = T)
  # make a checksum
  mdfiver::create_sha256_for_file(ereg_file_loc)
  # also create a covariates file
  covariates_columns <- c('ps_column', setdiff(colnames(mapping_ct_aggregated), 'ps_column'))
  covariates <- mapping_ct_aggregated[, ..covariates_columns]
  # create a location for it
  covar_file_loc <- paste(ereg_dir, 'covariates.tsv.gz', sep = '/')
  write.table(covariates, gzfile(covar_file_loc), sep = '\t', quote = F, row.names = F, col.names = T)
  mdfiver::create_sha256_for_file(covar_file_loc)
}

# we can also generate an smf from this
# smf <- unique(ps_to_ct[ , c('ps_column', ..sample_column)])
smf <- unique(ps_to_ct[ , c(..sample_column, 'ps_column')])
# set column names
colnames(smf) <- c('genotype_id', 'sample_id')

# save the smf as well
smf_loc <- paste(activity_ps_loc, 'smf.tsv.gz', sep = '/')
write.table(smf, gzfile(smf_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(smf_loc)
