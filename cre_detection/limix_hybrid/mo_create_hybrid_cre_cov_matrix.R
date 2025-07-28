#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_hybrid_cre_cov_matrix.R
# Function: create single-cell LIMIX covariate matrix files
# Example: 
# Rscript mo_create_hybrid_cre_cov_matrix.R \
# --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_monocyte_1_80_20240521.rds \
# --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/cov_matrix.txt \
# --mcolumns lane,condition_final
# --column sample_final
#
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(data.table)
library(optparse)


####################
# Functions        #
####################


####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- T


####################
# Main code        #
####################

# make command line options
option_list <- list(
  make_option(c("-i", "--in"), type="character", default=NULL, 
              help="input Seurat object to generate matrices for", metavar="character"),
  make_option(c("-o", "--out_cov"), type="character", default=NULL, 
              help="output covariates file", metavar="character"), 
  make_option(c("-m", "--mcolumns"), type="character", default=NULL, 
              help="comma separated string of metadata columns to include in the covariate file", metavar="character"), 
  make_option(c("-d", "--dcolumn"), type="character", default=NULL, 
              help="donor column to use for kinship file", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize some variables we will used
ct_object_loc <- NULL
out_cov <- NULL
out_ks <- NULL
mcolumns_string <- NULL
dcolumn <- NULL

if (debug) {
  # location of the cell type object
  ct_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_monocyte_1_80_20240521.rds'
  
  # location of the output
  out_cov <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/cov_matrix.tsv.gz'
  out_ks <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/low_rank_kinship.tsv.gz'
  
  # covariate annotation columns
  mcolumns_string <- 'lane,condition_final'
  
  # donor annotation
  dcolumn <- 'sample_final'
} else {
  # there are some things we cannot allow
  if (is.null(opt[['in']])) {
    stop('seurat file must be supplied')
  } else {
    # location of the cell type object
    ct_object_loc <- opt[['in']]
  }
  if (is.null(opt[['out_cov']])) {
    stop('output file must be supplied')
  } else {
    # location of the output
    out_cov <- opt[['out_cov']]
  }
  if (is.null(opt[['out_ks']])) {
    stop('output file must be supplied')
  } else {
    # location of the output
    out_ks <- opt[['out_ks']]
  }
  if (is.null(opt[['mcolumns']])) {
    stop('metadata columns string')
  } else {
    # location of the metadata columns to create covariates for
    columns_string <- opt[['mcolumns']]
  }
  if (is.null(opt[['dcolumn']])) {
    stop('metadata columns string')
  } else {
    # location of the output
    dcolumn <- opt[['dcolumn']]
  }
}

# get the metadata columns
metadata_columns <- strsplit(mcolumns_string, ',')[[1]]

# read the seurat object
seurat_object <- readRDS(ct_object_loc)

# extract the metadata
metadata <- seurat_object@meta.data

# clear memory
rm(seurat_object)

# extract the barcodes
barcodes_present <- rownames(metadata)
# save the in a list
binary_metadata_list <- list('barcode' = barcodes_present)

# check each metadata column
for (metadata_column in metadata_columns) {
  # get the unique values for that metadata column
  metadata_unique_values <- unique(metadata[[metadata_column]])
  # only include it if there is more than one level
  if (length(metadata_unique_values) > 1) {
    # order them
    metadata_unique_values <- metadata_unique_values[order(metadata_unique_values)]
    # we'll remove the last one, as it will automatically be that one if it is none of the others
    metadata_unique_values <- metadata_unique_values[-(length(metadata_unique_values))]
    # now check each of the unique values
    for (metadata_unique in metadata_unique_values) {
      # get the barcodes belonging to that unique metadata variable
      barcodes_this_value <- rownames(metadata[!is.na(metadata[[metadata_column]]) & metadata[[metadata_column]] == metadata_unique, ])
      # and get the 0/1 assignment for the barcodes
      binary_metadata_list[[metadata_unique]] <- as.numeric(barcodes_present %in% barcodes_this_value)
    }
  }
}

# make list of memberships into table
binary_metadata <- do.call('cbind', binary_metadata_list)
# with gzfile filehandle if needed
out_cov_out_fh <- out_cov
# and write
if (grepl('.gz$', out_cov)) {
  # gzip if ends with .gz
  out_fh <- gzfile(out_cov_out_fh)
}
write.table(binary_metadata, out_cov_out_fh, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum
mdfiver::create_md5_for_file(out_cov)


# repeat for sample assignment
binary_kinship_list <- list('barcode' = barcodes_present)
# get the unique values for that metadata column
donor_unique_values <- unique(metadata[[dcolumn]])
# order them
donor_unique_values <- donor_unique_values[order(donor_unique_values)]
# now check each of the unique values
for (donor_unique in donor_unique_values) {
  # get the barcodes belonging to that unique metadata variable
  barcodes_this_value <- rownames(metadata[!is.na(metadata[[dcolumn]]) & metadata[[dcolumn]] == donor_unique, ])
  # and get the 0/1 assignment for the barcodes
  binary_kinship_list[[donor_unique]] <- as.numeric(barcodes_present %in% barcodes_this_value)
}

# make list of memberships into table
binary_kinship <- do.call('cbind', binary_kinship_list)
# with gzfile filehandle if needed
kin_out_fh <- out_ks
# and write
if (grepl('.gz$', out_ks)) {
  # gzip if ends with .gz
  kin_out_fh <- gzfile(out_ks)
}
write.table(binary_kinship, kin_out_fh, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum
mdfiver::create_md5_for_file(out_ks)
