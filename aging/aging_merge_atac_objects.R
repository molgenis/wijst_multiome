#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: aging_merge_atac_objects.R
# Function: merge atac objects
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
# this database is needed for annotations
library(EnsDb.Hsapiens.v86)


####################
# Functions        #
####################



####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# get annotations from ensemble database
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# we use UCSC gencode
seqlevelsStyle(annotations) <- "UCSC"
# and this was aligned on build 38
genome(annotations) <- "hg38"


####################
# Main Code        #
####################


# set where we'll output the matrices
mo_disassembled_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/aging/objects/atac/disassembled/'

# make a mapping of what the cell type is called, and how it is named in the filesystem
ct_to_file <- list('B' = 'b', 'CD4T' = 'cd4t', 'CD8T' = 'cd8t', 'DC' = 'dc', 'monocyte' = 'monocyte', 'NK' = 'nk', 'plasmablast' = 'plasmablast', 'T_other' = 't_other')

# keep track of current matrix
current_matrix <- NULL
current_features <- NULL
current_barcodes <- NULL
current_metadata <- NULL

# check each of the cell types
for (cell_type in names(ct_to_file)) {
  # paste together the location of where to store the subset
  object_aging_loc <- paste0(mo_disassembled_loc, '/', cell_type, '/', 'peaks', '/')
  # show where we are
  message(paste('reading', cell_type, 'at', object_aging_loc))
  # get the location of the metadata
  metadata_loc <- paste(object_aging_loc, 'metadata.tsv.gz', sep = '/')
  # read the metadata
  metadata <- read.table(metadata_loc, header = T, sep = '\t', row.names = 1)
  # the barcodes
  barcodes_loc <- paste(object_aging_loc, 'barcodes.tsv.gz', sep = '/')
  barcodes <- read.table(barcodes_loc, header = F)$V1
  # the features
  features_loc <- paste(object_aging_loc, 'features.tsv.gz', sep = '/')
  features <- read.table(features_loc, header = F)$V1
  # get the matrix
  matrix_loc <- paste(object_aging_loc, 'matrix.mtx.gz', sep = '/')
  matrix <- Matrix::readMM(matrix_loc)
  # set the current matrix, features and metadata
  if (is.null(current_matrix)) {
    current_matrix <- matrix
    current_features <- features
    current_barcodes <- barcodes
    current_metadata <- metadata
  } else {
    # get which new features we dont have in the old matrix
    new_features_not_in_old <- setdiff(features, current_features)
    print(length(new_features_not_in_old))
    # and which old features are not in the new matrix
    old_features_not_in_new <- setdiff(current_features, features)
    # create zeroes for those
    new_features_expresssion_not_in_old <- Matrix(0, nrow = length(new_features_not_in_old), ncol = length(current_barcodes), sparse = TRUE)
    old_features_expression_not_in_new <- Matrix(0, nrow = length(old_features_not_in_new), ncol = length(barcodes), sparse = TRUE)
    # now add those to the two matrices
    current_matrix <- rbind(current_matrix, new_features_expresssion_not_in_old)
    matrix <- rbind(matrix, old_features_expression_not_in_new)
    # and the current order of the features
    current_features <- c(current_features, new_features_not_in_old)
    features <- c(features, old_features_not_in_new)
    # get the same order for features in both of them the same
    current_features_order <- order(current_features)
    features_order <- order(features)
    # now make sure we have that order
    current_features <- current_features[current_features_order]
    features <- features[features_order]
    current_matrix <- current_matrix[current_features_order, ]
    matrix <- matrix[features_order, ]
    # now merge the matrices, which can be done with a cbind, as we have added the missing features on both sides and ordered them the same
    current_matrix <- cbind(current_matrix, matrix)
    # add the barcodes, as we did a cbind
    current_barcodes <- c(current_barcodes, barcodes)
    # with the same for the metadata
    current_metadata <- rbind(current_metadata, metadata)
  }
}

# and write the chromatin data
chromatin_barcodes_loc <- paste0(mo_disassembled_loc, '/merged/barcodes.tsv.gz')
write.table(data.frame(x = current_barcodes), gzfile(chromatin_barcodes_loc), row.names = F, col.names = F, quote = F)
mdfiver::create_sha256_for_file(chromatin_barcodes_loc)
chromatin_features_loc <- paste0(mo_disassembled_loc, '/merged/features.tsv.gz')
write.table(data.frame(x = current_features), gzfile(chromatin_features_loc), row.names = F, col.names = F, quote = F)
mdfiver::create_sha256_for_file(chromatin_features_loc)
# write matrix
chromatin_matrix_loc <- paste0(mo_disassembled_loc, '/merged/matrix.mtx')
writeMM(current_matrix, chromatin_matrix_loc)
# Compress the MTX file using gzip
chromatin_matrix_gz_loc <- paste0(mo_disassembled_loc, '/merged/matrix.mtx.gz')
gzip(chromatin_matrix_loc, chromatin_matrix_gz_loc, overwrite = TRUE)
# and make md5 checksum
mdfiver::create_sha256_for_file(chromatin_matrix_gz_loc)
# add barcode as explicit column
current_metadata <- cbind(data.frame('barcode' = rownames(metadata)), current_metadata)
# write the metadata
metadata_loc <- paste0(mo_disassembled_loc, '/merged/metadata.tsv.gz')
# write it
write.table(current_metadata, gzfile(metadata_loc), row.names = F, col.names = T, quote = F, sep = '\t')
# and make an md5
mdfiver::create_sha256_for_file(metadata_loc)

# set row and column names for the chromatin assay
rownames(current_matrix) <- current_features
colnames(current_matrix) <- current_barcodes
# create the chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = current_matrix,
  sep = c("-", "-"),
  fragments = NULL,
  min.cells = 10,
  min.features = 200
)
# create object
seurat_object <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = current_metadata,
  project = 'wijst_multiome'
)
# set the annotations to the object now
Annotation(seurat_object) <- annotations
# make the location of the seurat object
signac_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/aging/objects/atac/assembled/aging_mo_atac_subset_20251022.rds'
# and write the object
saveRDS(seurat_object, signac_object_loc)
# make a checksum
mdfiver::create_sha256_for_file(signac_object_loc)
