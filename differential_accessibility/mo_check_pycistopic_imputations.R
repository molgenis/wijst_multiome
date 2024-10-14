#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_pycistopic_imputations.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# read the object
library(Seurat)
library(Signac)
# plotting
library(ggplot2)
library(cowplot)
# convert count matrices
library(Matrix)


####################
# Functions        #
####################


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 900 * 1000 * 1024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# the cell type
cell_type <- 'monocyte'

# location of the object
signac_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_monocyte_wstatus_1_80_20240709.rds'
signac_celltype <- readRDS(signac_object_loc)

# the location of the input matrices
imputed_matrix_loc <- paste('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/imputed_pycistopic_matrices/', cell_type, '/', sep = '')
# the prepend
imputed_matrix_prepend <- 'matrix_'
# the append
imputed_matrix_append <- '.mtx.gz'

# list the files in the directory
matrix_files <- list.files(imputed_matrix_loc, pattern = paste(imputed_matrix_prepend, '(\\d+)_(\\d+)', imputed_matrix_append, sep = ''))

# we'll store rowsums per chunk
rowsums_per_chunk <- list()

# go through each matrix
for (matrix_file in matrix_files) {
  # extract the region
  regions_string <- stringr::str_extract(matrix_file, '(\\d+)_(\\d+)')
  # split by underscore
  regions_vector <- regions_string[[1]]
  # read the accompanying features file
  features_matrix_loc <- paste(imputed_matrix_loc, 'features_', regions_string, '.tsv.gz', sep = '')
  features_matrix <- read.table(features_matrix_loc)$V1
  # read the matrix
  matrix_regions <- Matrix::readMM(paste(imputed_matrix_loc, matrix_file, sep = ''))
  # get the row sums
  sums_regions <- Matrix::rowSums(matrix_regions)
  # get the non-zero fraction
  nnz_fraction <- tabulate(matrix_regions@i + 1) / ncol(matrix_regions)
  # put in a nice dataframe with the region names
  df_sums_regions <- data.frame(region = features_matrix, total_count = sums_regions, nnz_fraction = nnz_fraction)
  # add that to the list
  rowsums_per_chunk[[matrix_file]] <- df_sums_regions
}

# merge all the regions
rowsums_imputed <- do.call('rbind', rowsums_per_chunk)
# now get the unimputed regions as well
rowsums_unimputed_vector <- Matrix::rowSums(signac_celltype@assays$peaks@counts)
# get the fraction of non-zeroes
nnz_fraction_unimputed <- tabulate(signac_celltype@assays$peaks@counts@i + 1, nrow(signac_celltype@assays$peaks@counts)) / ncol(signac_celltype@assays$peaks@counts)
# and make into a dataframe
rowsums_unimputed <- data.frame(region = names(rowsums_unimputed_vector), total_count = as.vector(unlist(rowsums_unimputed_vector)), nnz_fraction = nnz_fraction_unimputed)
# add whether or not this was also imputed
rowsums_unimputed[['imputation_included']] <- rowsums_unimputed[['region']] %in% gsub(':', '-', rowsums_imputed[['region']])
# make that into a character
rowsums_unimputed[['imputation_included']] <- as.character(rowsums_unimputed[['imputation_included']])

# plot the the counts
ggplot(rowsums_unimputed, aes(log10(total_count), fill = imputation_included)) + 
  geom_density(alpha = 0.2) + 
  scale_fill_manual(name = 'included in imputation', values = list('TRUE' = 'blue', 'FALSE' = 'red')) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('log10 sum of region counts') + 
  ylab('Density')
# and the non-zero fraction
ggplot(rowsums_unimputed, aes(-log10(nnz_fraction), fill = imputation_included)) + 
  geom_density(alpha = 0.2) + 
  scale_fill_manual(name = 'included in imputation', values = list('TRUE' = 'blue', 'FALSE' = 'red')) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('-log10 of non-zero fraction') + 
  ylab('Density')

# read the cpeaks annotation file
cpeaks_annotation_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_info.tsv'
cpeaks_annotation <- read.table(cpeaks_annotation_loc, header = T)
# add this annotation to the unimputed stats
rowsums_unimputed[['housekeeping']] <- cpeaks_annotation[match(rowsums_unimputed[['region']], paste(cpeaks_annotation[['chr_hg38']], cpeaks_annotation[['start_hg38']], cpeaks_annotation[['end_hg38']], sep = '-')), 'housekeeping']
rowsums_unimputed[['inferredElements']] <- cpeaks_annotation[match(rowsums_unimputed[['region']], paste(cpeaks_annotation[['chr_hg38']], cpeaks_annotation[['start_hg38']], cpeaks_annotation[['end_hg38']], sep = '-')), 'inferredElements']
