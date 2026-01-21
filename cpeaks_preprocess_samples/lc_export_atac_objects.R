#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lc_export_atac_objects.R
# Function: create an ATAC Seurat object to use in LONG-COVID data
############################################################################################################################


####################
# libraries        #
####################

# object loading
library(Seurat)
library(Signac)
# needed for checksum
library(mdfiver)
# need to export sparse matrix
library(Matrix)
# needed to zip
library(R.utils)
# this database is needed for annotations
library(EnsDb.Hsapiens.v86)

####################
# functions        #
####################


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 120 * 2000 * 2024^2)

# set seed
set.seed(7777)

# get annotations from ensemble database
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# we use UCSC gencode
seqlevelsStyle(annotations) <- "UCSC"
# and this was aligned on build 38
genome(annotations) <- "hg38"

####################
# Main Code        #
####################

# location of the atac objects
mo_objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/'
# start of ATAC object name
mo_object_prepend <- 'mo_cpeaks_filtered_'
mo_object_append <- '_wstatus_1_80_20240709.rds'

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
  mo_object_loc <- paste0(mo_objects_loc, '/', mo_object_prepend, ct_to_file[[cell_type]], mo_object_append)
  # show where we are
  message(paste('reading', cell_type, 'at', mo_object_loc))
  # read object
  mo_object_ct <- readRDS(mo_object_loc)
  # subset to the long covid samples
  mo_object_ct <- mo_object_ct[, !is.na(mo_object_ct@meta.data[['condition_final']]) & mo_object_ct@meta.data[['condition_final']] == 'UT']
  # and keep only the assigned long covid samples
  mo_object_ct <- mo_object_ct[, !is.na(mo_object_ct@meta.data[['LONG_COVID_method']]) & mo_object_ct@meta.data[['LONG_COVID_method']] == 'assigned']
  # remove some columns
  mo_object_ct@meta.data[, c('condition_imputed', 'celltype_imputed', 'final_condition', 'condition_sheet', 'condition_prev', 'LONG_COVID_method', 'sample_final', 'best_match_sample', 'second_match_sample', 'confined_best_match_sample', 'confined_second_match_sample', 'unconfined_best_match_sample', 'unconfined_second_match_sample')] <- NULL
  # read the metadata
  metadata <- mo_object_ct@meta.data
  # get the barcodes
  barcodes <- colnames(mo_object_ct)
  # get the accessibility
  matrix <- mo_object_ct@assays$peaks@counts
  # get the regions
  features <- rownames(matrix)
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

# load the sample mapping
mo_sample_mapping_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_covid_sample_mapping.tsv.gz'
mo_sample_mapping <- read.table(mo_sample_mapping_loc, header = T, sep = '\t')
# add these new columns
current_metadata[['snumber']] <- mo_sample_mapping[match(current_metadata[['realid']], mo_sample_mapping[['realid']]), 'snumber']
# remove the original annotation
current_metadata[['realid']] <- NULL

# and write the barcodes
chromatin_barcodes_loc <- paste0(mo_objects_loc, '/mo_atac_lc_study_only_20251217_barcodes.tsv.gz')
write.table(data.frame(x = current_barcodes), gzfile(chromatin_barcodes_loc), row.names = F, col.names = F, quote = F)
mdfiver::create_sha256_for_file(chromatin_barcodes_loc)
# the regions
chromatin_features_loc <- paste0(mo_objects_loc, '/mo_atac_lc_study_only_20251217_raw_features.tsv.gz')
write.table(data.frame(x = current_features), gzfile(chromatin_features_loc), row.names = F, col.names = F, quote = F)
mdfiver::create_sha256_for_file(chromatin_features_loc)
# write matrix
chromatin_matrix_loc <- paste0(mo_objects_loc, '/mo_rna_lc_study_only_20251217_matrix.mtx')
writeMM(current_matrix, chromatin_matrix_loc)
# zip matrix
gzip(chromatin_matrix_loc)
mdfiver::create_sha256_for_file(paste0(chromatin_matrix_loc, '.gz'))
# write the metadata
metadata_loc <- paste0(mo_objects_loc, 'mo_atac_lc_study_only_20251217_metadata.tsv.gz')
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
mo <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = current_metadata,
  project = 'wijst_multiome'
)
# set the annotations to the object now
Annotation(mo) <- annotations
# add some info
mo@misc[['processed_by']] <- 'r.oelen@umcg.nl'
mo@misc[['generated_by']] <- 'm.g.p.van.der.wijst@umcg.nl'
mo@misc[['exported_at']] <- '2025-12-17'
# and the object
mo_object_covid_loc <- paste0(mo_objects_loc, 'mo_atac_lc_study_only_20251217.rds')
# export
saveRDS(mo, mo_object_covid_loc)
# and make checksum
mdfiver::create_sha256_for_file(mo_object_covid_loc)
