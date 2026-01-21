#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: aging_subset_atac_objects.R
# Function: subset the ATAC objects of the multiomics study to the subset to test topic modelling for age/sex
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(data.table)
library(mdfiver)
library(Matrix)
library(R.utils)


####################
# Functions        #
####################

#' Write Data Sample
#'
#' This function writes expression data, chromatin data, and metadata to specified folders, ensuring that the sample name is POSIX-safe. It also creates MD5 checksums for the written files.
#'
#' @param chromatin_data A matrix containing the chromatin data.
#' @param metadata A data frame containing the metadata.
#' @param output_folder A character string specifying the output folder path.
#' @param chromatin_name A character string for the chromatin data folder name. Default is 'peaks'.
#' @param binarize_chromatin_matrix A logical value indicating whether to binarize the chromatin matrix. Default is TRUE.
#' @return An integer value of 0 upon successful completion.
#' @examples
#' \dontrun{
#' write_data_sample(expression_data, chromatin_data, metadata, "sample1", "/path/to/output")
#' }
#' @export
write_data_sample <- function(expression_data, chromatin_data, metadata, output_folder, chromatin_name='peaks', binarize_chromatin_matrix=T) {
  # the chromatin folder as well
  chromatin_data_folder <- output_folder
  dir.create(chromatin_data_folder, recursive = T)
  
  # add barcode as explicit column
  metadata <- cbind(data.frame('barcode' = rownames(metadata)), metadata)
  # write the metadata
  metadata_loc <- paste0(output_folder, '/metadata.tsv.gz')
  # write it
  write.table(metadata, gzfile(metadata_loc), row.names = F, col.names = T, quote = F, sep = '\t')
  # and make an md5
  mdfiver::create_md5_for_file(metadata_loc)
  
  # binarize the chromatin data if requested
  if (binarize_chromatin_matrix) {
    chromatin_data@x[chromatin_data@x > 1] <- 1
  }
  
  # and write the chromatin data
  chromatin_barcodes_loc <- paste0(chromatin_data_folder, 'barcodes.tsv.gz')
  write.table(data.frame(x = colnames(chromatin_data)), gzfile(chromatin_barcodes_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(chromatin_barcodes_loc)
  chromatin_features_loc <- paste0(chromatin_data_folder, 'features.tsv.gz')
  write.table(data.frame(x = rownames(chromatin_data)), gzfile(chromatin_features_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(chromatin_features_loc)
  # write matrix
  chromatin_matrix_loc <- paste0(chromatin_data_folder, 'matrix.mtx')
  writeMM(chromatin_data, chromatin_matrix_loc)
  # Compress the MTX file using gzip
  chromatin_matrix_gz_loc <- paste0(chromatin_data_folder, 'matrix.mtx.gz')
  gzip(chromatin_matrix_loc, chromatin_matrix_gz_loc, overwrite = TRUE)
  # and make md5 checksum
  mdfiver::create_md5_for_file(chromatin_matrix_gz_loc)
  return(0)
}


#' Deconstruct Seurat Object
#'
#' This function deconstructs a Seurat object by extracting RNA and chromatin data, subsetting to specified regions and genes, and writing the data for each sample to the specified output directory.
#'
#' @param seurat_object A Seurat object containing the data.
#' @param regions A character vector of regions to subset from the chromatin data.
#' @param output_dir A character string specifying the output directory path.
#' @param chromatin_assay A character string for the chromatin assay name. Default is 'peaks'.
#' @param chromatin_layer A character string for the chromatin data layer. Default is 'counts'.
#' @param binarize_chromatin_matrix A logical value indicating whether to binarize the chromatin matrix. Default is TRUE.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#' @return An integer value of 0 upon successful completion.
#' @examples
#' \dontrun{
#' deconstruct_seurat_object(seurat_object, regions, genes, "/path/to/output")
#' }
#' @export
deconstruct_seurat_object <- function(seurat_object, regions, output_dir, chromatin_assay='peaks', chromatin_layer='counts', binarize_chromatin_matrix=T, verbose=T) {
  # get the metadata
  metadata <- seurat_object@meta.data
  # get the chromatin data
  chromatin_data <- GetAssayData(seurat_object, assay = chromatin_assay, layer = chromatin_layer)
  # check which regions we have
  regions_have <- intersect(rownames(chromatin_data), regions)
  # report on what is missing if present
  if (verbose & length(regions_have) < length(regions)) {
    message(paste('only', as.character(length(regions_have)), 'out of ', length(regions), 'peaks present in data'))
  }
  # subset to those regions and genes
  chromatin_data <- chromatin_data[regions_have, ]
  # write the data
  write_data_sample(chromatin_data = chromatin_data, 
                        metadata = metadata, 
                        output_folder = output_dir, 
                        chromatin_name = chromatin_assay, 
                        binarize_chromatin_matrix = binarize_chromatin_matrix)
  return(0)
}




####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

####################
# Main Code        #
####################

# object location for Seurat
mo_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/'
# the prepend of the file
mo_object_prepend <- 'mo_cpeaks_filtered_'
# the append of the file
mo_object_append <- '_wstatus_1_80_20240709.rds'

# set where we'll output the matrices
mo_disassembled_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/aging/objects/atac/disassembled/'

# the location of the subset file
aging_subset_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/aging_subset_samples.tsv.gz'
# read the subset table
aging_subset <- fread(aging_subset_loc, header = F, sep = '\t')
# set column names
colnames(aging_subset) <- c('lane', 'participant', 'timepoint')

# make a mapping of what the cell type is called, and how it is named in the filesystem
ct_to_file <- list('B' = 'b', 'CD4T' = 'cd4t', 'CD8T' = 'cd8t', 'DC' = 'dc', 'monocyte' = 'monocyte', 'NK' = 'nk', 'plasmablast' = 'plasmablast', 'T_other' = 't_other')

# check each of the cell types
for (cell_type in names(ct_to_file)) {
  # get the file location
  object_loc <- paste0(mo_object_loc, '/', mo_object_prepend, ct_to_file[[cell_type]], mo_object_append)
  # let them know where we are
  message(paste('reading', cell_type, 'at', object_loc))
  # read the object
  object_ct <- readRDS(object_loc)
  # subset the object
  object_ct <- object_ct[, 
           paste(object_ct@meta.data[['lane']], object_ct@meta.data[['sample_final']], object_ct@meta.data[['condition_final']]) 
           %in%
           paste(aging_subset[['lane']], aging_subset[['participant']], aging_subset[['timepoint']])]
  # paste together the location of where to store the subset
  object_aging_loc <- paste0(mo_disassembled_loc, '/', cell_type, '/', 'peaks', '/')
  # make the directory
  dir.create(object_aging_loc, recursive = T)
  # write the object
  deconstruct_seurat_object(
    seurat_object = object_ct, regions = rownames(object_ct), output_dir = object_aging_loc, chromatin_assay='peaks', chromatin_layer='counts', binarize_chromatin_matrix=F, verbose=T
  )
}
