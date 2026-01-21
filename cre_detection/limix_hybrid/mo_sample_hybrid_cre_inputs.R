#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_sample_hybrid_cre_inputs.R
# Function: sample single-cell LIMIX input files
# Example: 
# Rscript mo_sample_hybrid_cre_inputs.R \
# --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD8T/ \
# --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input_sampled/sampling2/L1/CD8T/ \
# --seed 1337 \
# --fraction 0.5
#
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(optparse)


####################
# Functions        #
####################


subset_matrices <- function(input_directory, output_directory, barcodes_to_keep, matrices=c('accessibility.tsv.gz', 'expression.tsv.gz')) {
  # read each of the files
  for (matrix_file in matrices) {
    # paste the file location together
    matrix_input_loc <- paste(input_directory, matrix_file, sep = '/')
    # read that file
    matrix_input <- fread(matrix_input_loc, sep = '\t')
    # set the first column as 'X'
    colnames(matrix_input)[1] <- 'X'
    # keep these columns
    barcodes_kept_and_present <- colnames(matrix_input)[colnames(matrix_input) %in% barcodes_to_keep]
    # add the x back
    columns_to_keep <- c('X', barcodes_kept_and_present)
    # now subset the matrix to the 'X' and the barcodes
    matrix_input <- matrix_input[, ..columns_to_keep]
    # paste the output loc
    matrix_out_loc <- paste(output_directory, matrix_file, sep = '/')
    # with gzfile filehandle if needed
    matrix_out_fh <- matrix_out_loc
    # and write
    if (grepl('.gz$', matrix_out_loc)) {
      # gzip if ends with .gz
      matrix_out_fh <- gzfile(matrix_out_loc)
    }
    write.table(matrix_input, matrix_out_fh, sep = '\t', row.names = F, col.names = T, quote = F)
    # also make a checksum
    mdfiver::create_md5_for_file(matrix_out_loc)
  }
  return(0)
}


####################
# Settings         #
####################

# whether we are in debug mode
debug <- F


####################
# Main code        #
####################

# these we'll keep default for now
smf_filename <- 'smf.tsv.gz'
cov_filename <- 'cov_matrix.txt'
ks_filename <- 'low_rank_kinship.txt'


# make command line options
option_list <- list(
  make_option(c("-i", "--in"), type="character", default=NULL, 
              help="input directory of chunks", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory of chunks", metavar="character"), 
  make_option(c("-s", "--seed"), type="numeric", default=NULL, 
              help="the seed to use for the sampling (will create one if not supplied)", metavar="numeric"),
  make_option(c("-f", "--fraction"), type="numeric", default=0.5, 
              help="fraction of the total data to sample from", metavar="numeric")
)


# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize our variables
input_dir <- NULL
output_dir <- NULL
seed <- NULL
sample_fraction <- NULL

if (debug) {
  # location of the region-to-gene files
  input_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD8T/'
  
  # read the location of the genes
  output_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input_sampled/sampling1/L1/CD8T/'
  
  # location of the cell type object
  seed <- 1337
  
  # location of the output
  sample_fraction <- 0.5
  
} else {
  # the input location
  if (is.null(opt[['in']])) {
    stop('input must be supplied')
  } else {
    input_dir <- opt[['in']]
  }
  # the output location
  if (is.null(opt[['out']])) {
    stop('output must be supplied')
  } else {
    output_dir <- opt[['out']]
  }
  # get the seed
  seed <- opt[['seed']]
  # and the fraction of output
  sample_fraction <- opt[['fraction']]
}

# the fraction has to be less than 1
if (!(sample_fraction < 1)) {
  stop('sample fraction must be smaller than 1')
}


# get the seed if null
if (is.null(seed)) {
  seed <- sample.int(.Machine$integer.max, 1)[1]
}
# set the seed
set.seed(seed)

# write the seed
seed_loc <- paste(output_dir, 'seed.txt.gz', sep = '/')
write.table(data.frame(x = c(seed)), gzfile(seed_loc), row.names = F, col.names = F, quote = F)
mdfiver::create_md5_for_file(seed_loc)

# get smf loc
smf_loc <- paste(input_dir, smf_filename, sep = '/')
# read that smf
smf <- fread(smf_loc, header = T, sep = '\t')
# grab the barcodes
all_barcodes <- smf[['phenotype_id']]
# get the number of barcodes
n_barcodes <- length(all_barcodes)
# get how many we would like to get from there
n_barcodes_to_keep <- round(sample_fraction * n_barcodes, digits = 0)
# randomly sample indices
barcode_indices <- sample(n_barcodes, n_barcodes_to_keep)
# and get those barcodes
barcodes_to_keep <- all_barcodes[barcode_indices]

# subset the smf
smf <- smf[smf[['phenotype_id']] %in% barcodes_to_keep, ]
# paste the output loc
smf_out_loc <- paste(output_dir, smf_filename, sep = '/')
# with gzfile filehandle if needed
smf_out_fh <- smf_out_loc
# and write
if (grepl('.gz$', smf_out_loc)) {
  # gzip if ends with .gz
  smf_out_fh <- gzfile(smf_out_loc)
}
write.table(smf, smf_out_fh, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum
mdfiver::create_md5_for_file(smf_out_loc)

# get covariate matrix
covariate_loc <- paste(input_dir, cov_filename, sep = '/')
# read that covariate file
covariates <- fread(covariate_loc, header = T, sep = '\t')
# set the first column as 'X'
colnames(covariates)[1] <- 'X'
# keep only the entries in the first column, that is in our barcodes list
covariates <- covariates[covariates[['X']] %in% barcodes_to_keep, ]
# paste the output loc
covariates_out_loc <- paste(output_dir, cov_filename, sep = '/')
# with gzfile filehandle if needed
covariates_out_fh <- covariates_out_loc
# and write
if (grepl('.gz$', covariates_out_loc)) {
  # gzip if ends with .gz
  covariates_out_fh <- gzfile(covariates_out_loc)
}
write.table(covariates, covariates_out_fh, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum
mdfiver::create_md5_for_file(covariates_out_loc)

# get kinship matrix
ks_loc <- paste(input_dir, ks_filename, sep = '/')
# read that kinship file
ks <- fread(ks_loc, header = T, sep = '\t')
# set the first column as 'X'
colnames(ks)[1] <- 'X'
# keep only the entries in the first column, that is in our barcodes list
ks <- ks[ks[['X']] %in% barcodes_to_keep, ]
# paste the output loc
ks_out_loc <- paste(output_dir, ks_filename, sep = '/')
# with gzfile filehandle if needed
ks_out_fh <- ks_out_loc
# and write
if (grepl('.gz$', ks_out_loc)) {
  # gzip if ends with .gz
  ks_out_fh <- gzfile(ks_out_loc)
}
write.table(covariates, ks_out_fh, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum
mdfiver::create_md5_for_file(ks_out_loc)

# check each of the folders
chunk_folders <- list.dirs(input_dir, recursive = F, full.names = F)
# and filter based on the pattern
chunk_folders <- chunk_folders[grepl('chr\\d+-\\d+-\\d+', chunk_folders)]
# do each chunk
for (chunk_folder in chunk_folders) {
  # get to and from folders
  chunk_folder_from <- paste(input_dir, chunk_folder, sep = '/')
  chunk_folder_to <- paste(output_dir, chunk_folder, sep = '/')
  # make the 'to' folder
  dir.create(chunk_folder_to, recursive = T, showWarnings = F)
  # now do the run
  subset_matrices(
    input_directory = chunk_folder_from, 
    output_directory = chunk_folder_to, 
    barcodes_to_keep = barcodes_to_keep
  )
}
