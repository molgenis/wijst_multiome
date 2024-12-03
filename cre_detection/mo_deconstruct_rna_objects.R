#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_deconstruct_rna_objects.R
# Function: deconstruct Seurat RNA objects so that they can be used to create scanpy objects for the scenic pipeline
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Matrix)


####################
# Functions        #
####################

#' Deconstruct Seurat Object
#'
#' This function deconstructs a Seurat object, extracting metadata, cell barcodes,
#' and count data, and writes them to specified output folders.
#'
#' @param seurat_object_loc Character. The file path to the Seurat object (.rds file).
#' @param output_folder Character. The directory where the deconstructed files will be saved. Default is the current directory.
#' @return Integer. Returns 0 upon successful completion.
#' @examples
#' deconstruct_seurat_object(seurat_object_loc = "path/to/seurat_object.rds", 
deconstruct_seurat_object <- function(seurat_object_loc, output_folder='./') {
    # load the Seurat object
    seurat_object <- readRDS(seurat_object_loc)

    # remove the empty cells
    seurat_object <- seurat_object[, colSums(seurat_object) > 0]
    # remove empty cells
    #seurat_object <- seurat_object[rowSums(seurat_object) > 0, , drop = F]

    # create the root folder
    dir.create(output_folder, recursive = T)

    # extract the metadata
    metadata <- seurat_object@meta.data
    # add the rownames as cell IDs
    metadata <- cbind(data.frame(cell_id = rownames(metadata)), metadata)
    # write the metadata
    write.table(metadata, gzfile(paste(output_folder, '/metadata.tsv.gz', sep = '')), row.names = F, col.names = T, quote = F, sep = '\t')

    # get the features for this assay
    features_any_raw <- data.frame(seurat_object@assays$RNA@features)
    # and get the ones present in the counts layer, as these should be interpreted as the row names or features.tsv.gz
    features_present_raw <- rownames(features_any_raw[features_any_raw[['counts']] == T, , drop = F])
    # make them unique
    features_present_raw <- make.unique(features_present_raw, sep = ".")
    # create the directory
    dir.create(paste(output_folder, '/RNA/counts/', sep = ''), recursive = T)
    # write these
    write.table(data.frame(x = features_present_raw), gzfile(paste(output_folder, '/RNA/counts/features.tsv.gz', sep = '')), row.names = F, col.names = F, quote = F)
    # extract the raw count data
    counts_raw <- seurat_object@assays$RNA@layers$counts
    # extract the cell barcodes
    cell_barcodes_raw <- colnames(seurat_object)
    # write those barcodes
    write.table(data.frame(x = cell_barcodes_raw), gzfile(paste(output_folder, '/RNA/counts/barcodes.tsv.gz', sep = '')), row.names = F, col.names = F, quote = F)
    # write the matrix
    Matrix::writeMM(obj = counts_raw, file = paste(output_folder, '/RNA/counts/matrix.mtx', sep = ''))
    
    # extract the normalized counts as well, with the same aproach
    features_present_sct <- rownames(seurat_object@assays$SCT@counts)
    features_present_sct <- make.unique(features_present_sct, sep = ".")
    dir.create(paste(output_folder, '/SCT/counts/', sep = ''), recursive = T)
    write.table(data.frame(x = features_present_raw), gzfile(paste(output_folder, '/SCT/counts/features.tsv.gz', sep = '')), row.names = F, col.names = F, quote = F)
    counts_sct <- seurat_object@assays$SCT@counts
    cell_barcodes_sct <- colnames(counts_sct)
    write.table(data.frame(x = cell_barcodes_sct), gzfile(paste(output_folder, '/SCT/counts/barcodes.tsv.gz', sep = '')), row.names = F, col.names = F, quote = F)
    Matrix::writeMM(obj = counts_sct, file = paste(output_folder, '/SCT/counts/matrix.mtx', sep = ''))

    # return 0 upon success
    return(0)
}

#' Write Deconstructed Seurat Object
#'
#' This function processes Seurat object files from a specified input directory,
#' deconstructs them, and writes the results to an output directory.
#'
#' @param input_dir Character. The directory containing the Seurat object files.
#' @param output_dir Character. The directory where the deconstructed files will be saved. Default is the current directory.
#' @param seurat_objects_prepend Character. The prefix pattern for Seurat object filenames. Default is '^mo_all_20240619_seuratv5_annotated_agesexcovid_'.
#' @param seurat_objects_append Character. The suffix pattern for Seurat object filenames. Default is '\\.rds$'.
#' @return Integer. Returns 0 upon successful completion.
#' @examples
#' write_deconstructed_seurat_object(input_dir = "path/to/input", output_dir = "path/to/output")
write_deconstructed_seurat_object <- function(input_dir, output_dir='./', seurat_objects_prepend='^mo_all_20240619_seuratv5_annotated_agesexcovid_', seurat_objects_append='\\.rds$', verbose = F) {
    # create the file list
    seurat_files <- list.files(input_dir, full.names = FALSE, recursive = FALSE)
    # create the matching regex
    seurat_file_regex <- paste(seurat_objects_prepend, '.*', seurat_objects_append, sep = '')
    # filter the seurat file list
    seurat_files <- seurat_files[grepl(seurat_file_regex, seurat_files)]
    # now check each file
    for (seurat_file in seurat_files) {
        # extract the 'cell type'
        seurat_file_name <- gsub(seurat_objects_prepend, '', seurat_file)
        seurat_file_name <- gsub(seurat_objects_append, '', seurat_file_name)
        # use that to make the output directory
        output_deconstruction_folder <- paste(output_dir, '/', seurat_file_name, '/', sep = '')
        # get the full path to the input as well
        input_seurat_file <- paste(input_dir, '/', seurat_file, sep = '')
        # verbosity
        if (verbose) {
            message(paste('deconstructing', input_seurat_file))
        }
        # then do the deconstruction
        deconstruct_seurat_object(input_seurat_file, output_folder = output_deconstruction_folder)
    }
    return(0)
}

####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
set.seed(7777)


####################
# Main Code        #
####################

# location of Seurat files
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
# start and end of each file
seurat_objects_prepend <- '^mo_all_20240619_seuratv5_annotated_agesexcovid_'
seurat_objects_append <- '\\.rds$'
# where to place the output
output_deconstructed_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/'
# now do the deconstruction
write_deconstructed_seurat_object(
    input_dir = seurat_objects_loc, 
    output_dir = output_deconstructed_loc, 
    seurat_objects_prepend = seurat_objects_prepend, 
    seurat_objects_append = seurat_objects_append, 
    verbose = T)
