#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_split_sample_and_celltype.R
# Function: take a multimodal Seurat object and split it per donor into matrices for the RNA and accessibility
# Example
# Rscript mo_split_sample_and_celltype.R \
#   --cell_type DC \
#   --cell_type_column celltype_imputed_lowerres \
#   --seurat_object_path /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_dc_1_80_20240521.rds \
#   --cre_pairs_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/coeqtl/trial_run/cre_lists/monocyte_eregulon_pairs.tsv.gz \
#   --seurat_assignment_column sample_final,lane \
#   --output_folder /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/coeqtl/trial_run/matrices/ \
#   --no_binarize_atac \
#   --quietly
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(Matrix)
library(R.utils)
library(optparse)
library(data.table)


####################
# Functions        #
####################

#' Make a POSIX-safe String
#'
#' This function takes an input string and replaces any character that is not a letter, digit, hyphen, or underscore with an underscore.
#'
#' @param input_string A character string that needs to be converted to a POSIX-safe format.
#' @return A character string where all non-POSIX-safe characters are replaced with underscores.
#' @examples
#' make_posix_safe("example string!") # Returns "example_string_"
#' make_posix_safe("another@string#") # Returns "another_string_"
#' @export
make_posix_safe <- function(input_string) {
  # replace any character that is not a letter, digit, hyphen, or underscore with an underscore
  posix_safe_string <- gsub("[^a-zA-Z0-9_-]", "_", input_string)
  return(posix_safe_string)
}

#' Write Data Sample
#'
#' This function writes expression data, chromatin data, and metadata to specified folders, ensuring that the sample name is POSIX-safe. It also creates MD5 checksums for the written files.
#'
#' @param expression_data A matrix containing the expression data.
#' @param chromatin_data A matrix containing the chromatin data.
#' @param metadata A data frame containing the metadata.
#' @param sample_name A character string representing the sample name.
#' @param output_folder A character string specifying the output folder path.
#' @param expression_name A character string for the expression data folder name. Default is 'RNA'.
#' @param chromatin_name A character string for the chromatin data folder name. Default is 'peaks'.
#' @param binarize_chromatin_matrix A logical value indicating whether to binarize the chromatin matrix. Default is TRUE.
#' @return An integer value of 0 upon successful completion.
#' @examples
#' \dontrun{
#' write_data_sample(expression_data, chromatin_data, metadata, "sample1", "/path/to/output")
#' }
#' @export
write_data_sample <- function(expression_data, chromatin_data, metadata, sample_name, output_folder, expression_name='RNA', chromatin_name='peaks', binarize_chromatin_matrix=T) {
  # make the sample name posix safe
  sample_name_posix <- make_posix_safe(sample_name)
  # warn if this make the name different
  if (sample_name != sample_name_posix) {
    warning(paste('sample name was not POSIX safe,', sample_name, 'was renamed to', sample_name_posix))
  }
  # create folder for the expression data
  expression_data_folder <- paste0(output_folder, '/', sample_name_posix, '/', expression_name, '/')
  # the chromatin folder as well
  chromatin_data_folder <- paste0(output_folder, '/', sample_name_posix, '/', chromatin_name, '/')
  # create the directories
  dir.create(expression_data_folder, recursive = T)
  dir.create(chromatin_data_folder, recursive = T)
  
  # add barcode as explicit column
  metadata <- cbind(data.frame('barcode' = rownames(metadata)), metadata)
  # write the metadata
  metadata_loc <- paste0(output_folder, '/', sample_name_posix, '/metadata.tsv.gz')
  # write it
  write.table(metadata, gzfile(metadata_loc), row.names = F, col.names = T, quote = F, sep = '\t')
  # and make an md5
  mdfiver::create_md5_for_file(metadata_loc)
  
  # write the expression data
  expression_barcodes_loc <- paste0(expression_data_folder, 'barcodes.tsv.gz')
  write.table(data.frame(x = colnames(expression_data)), gzfile(expression_barcodes_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(expression_barcodes_loc)
  expression_features_loc <- paste0(expression_data_folder, 'features.tsv.gz')
  write.table(data.frame(x = rownames(expression_data)), gzfile(expression_features_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(expression_features_loc)
  # write matrix
  expression_matrix_loc <- paste0(expression_data_folder, 'matrix.mtx')
  writeMM(expression_data, expression_matrix_loc)
  # Compress the MTX file using gzip
  expression_matrix_gz_loc <- paste0(expression_data_folder, 'matrix.mtx.gz')
  gzip(expression_matrix_loc, expression_matrix_gz_loc, overwrite = TRUE)
  # and make md5 checksum
  mdfiver::create_md5_for_file(expression_matrix_gz_loc)
  
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
#' @param genes A character vector of genes to subset from the RNA data.
#' @param output_dir A character string specifying the output directory path.
#' @param rna_assay A character string for the RNA assay name. Default is 'RNA'.
#' @param chromatin_assay A character string for the chromatin assay name. Default is 'peaks'.
#' @param rna_layer A character string for the RNA data layer. Default is 'counts'.
#' @param chromatin_layer A character string for the chromatin data layer. Default is 'counts'.
#' @param seurat_assignment_column A character string for the column in metadata that contains sample assignments. Default is 'sample_final'.
#' @param binarize_chromatin_matrix A logical value indicating whether to binarize the chromatin matrix. Default is TRUE.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#' @return An integer value of 0 upon successful completion.
#' @examples
#' \dontrun{
#' deconstruct_seurat_object(seurat_object, regions, genes, "/path/to/output")
#' }
#' @export
deconstruct_seurat_object <- function(seurat_object, regions, genes, output_dir, rna_assay='RNA', chromatin_assay='peaks', rna_layer='counts', chromatin_layer='counts', seurat_assignment_column='sample_final', binarize_chromatin_matrix=T, verbose=T) {
  # get the metadata
  metadata <- seurat_object@meta.data
  # get the RNA data
  expression_data <- GetAssayData(seurat_object, assay = rna_assay, layer = rna_layer)
  # get the chromatin data
  chromatin_data <- GetAssayData(seurat_object, assay = chromatin_assay, layer = chromatin_layer)
  # check which regions we have
  regions_have <- intersect(rownames(chromatin_data), regions)
  # report on what is missing if present
  if (verbose & length(regions_have) < length(regions)) {
    message(paste('only', as.character(length(regions_have)), 'out of ', length(regions), 'peaks present in data'))
  }
  # check which genes we have
  genes_have <- intersect(rownames(expression_data), genes)
  # report on what is missing if present
  if (verbose & length(genes_have) < length(genes)) {
    message(paste('only', as.character(length(genes_have)), 'out of ', length(genes), 'genes present in data'))
  }
  # subset to those regions and genes
  expression_data <- expression_data[genes_have, ]
  chromatin_data <- chromatin_data[regions_have, ]
  # check each sample
  samples_present <- metadata[[seurat_assignment_column]]
  # and remove entries where we have no data
  samples_present <- samples_present[!is.na(samples_present)]
  # check each sample
  for (sample_name in samples_present) {
    if (verbose) {
      message(paste('starting sample', sample_name))
    }
    # get indices for sample
    indices_sample <- !is.na(metadata[[seurat_assignment_column]]) & metadata[[seurat_assignment_column]] == sample_name
    # subset the data to that sample
    expression_data_sample <- expression_data[, indices_sample]
    chromatin_data_sample <- chromatin_data[, indices_sample]
    metadata_sample <- metadata[indices_sample, ]
    # check if we have enough data for this sample
    if (is.null(nrow(metadata_sample)) | nrow(metadata_sample) < 2) {
      warning(paste('skipping', sample_name, 'due to too few cells'))
    }
    else if(is.null(nrow(expression_data_sample)) | nrow(expression_data_sample) < 2) {
      warning(paste('skipping', sample_name, 'due to too few expressing genes'))
    }
    else if(is.null(nrow(chromatin_data_sample)) | nrow(chromatin_data_sample) < 2) {
      warning(paste('skipping', sample_name, 'due to too few regions with observations'))
    }
    else {
      # write for this sample
      write_data_sample(expression_data = expression_data_sample, 
                        chromatin_data = chromatin_data_sample, 
                        metadata = metadata_sample, 
                        sample_name = sample_name, 
                        output_folder = output_dir, 
                        expression_name = rna_assay, 
                        chromatin_name = chromatin_assay, 
                        binarize_chromatin_matrix = binarize_chromatin_matrix)
    }
    if(verbose) {
      message(paste('finished sample', sample_name))
    }
  }
  return(0)
}


####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-c", "--cell_type"), type="character", default=NULL,
              help="cell type working on", metavar="character"),
  make_option(c("-t", "--cell_type_column"), type="character", default=NULL,
              help="cell type in metadata that has the cell type", metavar="character"),
  make_option(c("-s", "--seurat_object_path"), type="character", default=NULL,
              help="location to Seurat object (rds)", metavar="character"), 
  make_option(c("-e", "--cre_pairs_loc"), type="character", default=NULL,
              help="location of tab separated table, that has the region-gene pairs to test", metavar="character"), 
  make_option(c("-a", "--seurat_assignment_column"), type="character", default=NULL,
              help="name of column or comma separated list of columns in Seurat metadata that has the assignment of the cell to a sample", metavar="character"), 
  make_option(c("-o", "--output_folder"), type="character", default=NULL,
              help="path to the folder to store the sample-specific matrices ", metavar="character"), 
  make_option(c("-n", "--chromatin_assay"), type="character", default='peaks',
              help="name of chromatin assay in Seurat object", metavar="character"), 
  make_option(c("-r", "--rna_assay"), type="character", default='RNA',
              help="name of RNA assay in Seurat object", metavar="character"), 
  make_option(c("-i", "--chromatin_layer"), type="character", default='counts',
              help="name of chromatin layer in Seurat object to use for given assay", metavar="character"), 
  make_option(c("-l", "--rna_layer"), type="character", default='counts',
              help="name of RNA layer in Seurat object to use for given assay", metavar="character"), 
  make_option(c("-q", "--quietly"), action="store_true", default=FALSE, 
              help="do not print progress messages"), 
  make_option(c("-b", "--no_binarize_atac"), action="store_true", default=FALSE,
              help="do not binarize the ATAC data to 0/1")
)

# initialize optparser
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# check all arguments
if (is.null(opt[['cell_type']])) {
  stop('-c/--cell_type not supplied, this is a mandatory parameter')
}
if (is.null(opt[['cell_type_column']])) {
  stop('-t/--cell_type_column not supplied, this is a mandatory parameter')
}
if (is.null(opt[['seurat_object_path']])) {
  stop('-s/--seurat_object_path not supplied, this is a mandatory parameter')
}
if (is.null(opt[['cre_pairs_loc']])) {
  stop('-e/--cre_pairs_loc not supplied, this is a mandatory parameter')
}
if (is.null(opt[['seurat_assignment_column']])) {
  stop('-ac/--seurat_assignment_column not supplied, this is a mandatory parameter')
}
if (is.null(opt[['output_folder']])) {
  stop('-o/--output_folder not supplied, this is a mandatory parameter')
}

# fetch options
cell_type <- opt[['cell_type']]
cell_type_column <- opt[['cell_type_column']]
seurat_object_path <- opt[['seurat_object_path']]
cre_pairs_loc <- opt[['cre_pairs_loc']]
seurat_assignment_column <- opt[['seurat_assignment_column']]
output_folder <- opt[['output_folder']]
chromatin_assay <- opt[['chromatin_assay']]
rna_assay <- opt[['rna_assay']]
chromatin_layer <- opt[['chromatin_layer']]
rna_layer <- opt[['rna_layer']]
quietly <- opt[['quietly']]
no_binarize_atac <- opt[['no_binarize_atac']]
# the last two are easier if we flip them
verbose <- !quietly
binarize_atac <- !no_binarize_atac

# check if the output folder exists
if (!dir.exists(output_folder)) {
  stop(paste('output directory', output_folder, 'does not exist'))
}
# try to load the CREs to test
cre_pairs <- NULL
if (file.exists(cre_pairs_loc)) {
  cre_pairs <- fread(cre_pairs_loc, header = F, sep = '\t')
} else {
  stop(paste('cre-pair location file', cre_pairs_loc, 'does not exist'))
}
# extract from those the regions and genes
regions_to_use <- unique(cre_pairs[['V1']])
genes_to_use <- unique(cre_pairs[['V2']])
# if there are none, there is nothing to do
if (length(genes_to_use) == 0) {
  stop(paste('no genes present in', cre_pairs_loc))
}
if (length(regions_to_use) == 0) {
  stop(paste('no regions present in', cre_pairs_loc))
}

# try to load the Seurat object next
seurat_object <- NULL
if (file.exists(seurat_object_path)) {
  seurat_object <- readRDS(seurat_object_path)
} else {
  stop(paste(seurat_object_path, 'does not exist'))
}

# check if the Seurat object has both modalities
if (!rna_assay %in% names(seurat_object@assays)) {
  stop(paste0('RNA assay\'', rna_assay, '\', not present in Seurat object'))
}
if (!chromatin_assay %in% names(seurat_object@assays)) {
  stop(paste0('chromatin assay\'', chromatin_assay, '\', not present in Seurat object'))
}

# next, check if the cell type and assignment columns are present
if (!cell_type_column %in% colnames(seurat_object@meta.data)) {
  stop(paste0('cell type column \'', cell_type_column, '\', not present in Seurat object'))
}
# check if there is a comma
if (grepl(',', seurat_assignment_column)) {
  # if there is, we split by that comma
  seurat_assignment_columns <- strsplit(seurat_assignment_column, ',')[[1]]
  # make a new column that is all of the other columns, by first taking the first column
  seurat_object@meta.data[['aggregate_columns']] <- seurat_object@meta.data[[seurat_assignment_columns[[1]]]]
  # then check all other columns
  for (seurat_column_i in 2:length(seurat_assignment_columns)) {
    # get the column
    seurat_column <- seurat_assignment_columns[seurat_column_i]
    # check if the column exists
    if (seurat_column %in% colnames(seurat_object@meta.data)) {
      seurat_object@meta.data[['aggregate_columns']] <- paste(seurat_object@meta.data[['aggregate_columns']], seurat_object@meta.data[[seurat_column]], sep = '_')
    }
    else {
      stop(paste0('sample assignment column \'', seurat_column, '\', not present in Seurat object'))
    }
    # set that new column as the new assignment column
    seurat_assignment_column <- 'aggregate_columns'
  }
} else {
  # if there is no comma, we use the column as-is
  if (!seurat_assignment_column %in% colnames(seurat_object@meta.data)) {
    stop(paste0('sample assignment column \'', seurat_assignment_column, '\', not present in Seurat object'))
  }
}


# check if we have the cell type requested
if (!(cell_type %in% seurat_object@meta.data[[cell_type_column]])) {
  stop(paste0('cell type\'', cell_type, '\' not present in metadata using cell type column\'', cell_type_column, '\''))
}

# first step, subset the data to this cell type
seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[[cell_type_column]]) & seurat_object@meta.data[[cell_type_column]] == cell_type]

# set with the cell type as output folder
output_folder_full <- paste0(output_folder, '/', cell_type, '/')
dir.create(output_folder_full, recursive = T)

# start the procedure
deconstruct_seurat_object(
  seurat_object = seurat_object, 
  regions = regions_to_use, 
  genes = genes_to_use, 
  output_dir = output_folder_full, 
  rna_assay = rna_assay, 
  chromatin_assay = chromatin_assay, 
  rna_layer = rna_layer, 
  chromatin_layer = chromatin_layer, 
  seurat_assignment_column = seurat_assignment_column, 
  binarize_chromatin_matrix = binarize_atac, 
  verbose = verbose)

if (verbose) {
  message(paste(seurat_object_path, 'finished'))
}