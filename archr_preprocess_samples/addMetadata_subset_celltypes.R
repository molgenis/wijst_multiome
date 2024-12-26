#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Martijn van der Werf, Roy Oelen
# Name: addMetadata_subset_celltypes.R
# Function: add celltype annotation metadata from the Seurat RNA object
#
############################################################################################################################

setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)
library(Seurat)
library(stringr)

addArchRThreads(threads=4)
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject('Project_all_cells')

# Load Azimuth cell type annotations for all ATAC cells
barcode_celltype_anno <- read.csv('../cell_type_assignment/azimuth/10x_multiome_PBMCs/mo_azimuth_ct_10xmultiome.tsv', sep='\t')
azimuth_celltype_anno <- read.csv('../cell_type_assignment/azimuth/celltypes_10x_ref_to_lowerres.tsv', sep='\t')
# Add lower resolution cell types for all 1 milion cells
barcode_celltype_anno$cell_type_lower <- azimuth_celltype_anno$ten_x_lower[match(barcode_celltype_anno$predicted.mo_10x_cell_type, azimuth_celltype_anno$ten_x)]
barcode_celltype_anno$archr_barcode <- paste0(str_sub(barcode_celltype_anno$barcode, start=-12), "#", substr(barcode_celltype_anno$barcode, start=1, stop=16), "-1")

# Add cell type annotation, both resolutions
proj@cellColData$cell_type <- barcode_celltype_anno$predicted.mo_10x_cell_type[match(rownames(proj@cellColData), barcode_celltype_anno$archr_barcode)]
proj@cellColData$cell_type_lowerres <- barcode_celltype_anno$cell_type_lower[match(rownames(proj@cellColData), barcode_celltype_anno$archr_barcode)]

# Save project 
saveArchRProject(proj, 'Project_all_cells')

#########################################
# Subsetting project
#########################################

# Load seurat object for annotation
seurat_obj <- readRDS('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240223_seuratv5_normalized.rds')
metadata <- seurat_obj@meta.data

metadata$archr_barcode <- paste0(metadata$batch, "#", metadata$barcode_1)
proj@cellColData['barcode_archr'] <- rownames(proj@cellColData)

# 700,000 cells left
metadata_common <- metadata[match(metadata$archr_barcode, rownames(proj@cellColData)),]

# Subset project ArchR for common barcodes
proj_common <- proj[rownames(proj@cellColData) %in% metadata_common$archr_barcode]
saveArchRProject(proj_common, 'Project_common_atac_rna_cells')

# Add cell type annotation
proj_common@cellColData['celltype'] <- metadata_common$celltype_imputed_lowerres[match(rownames(proj_common@cellColData), metadata_common$archr_barcode)] 
proj_common@cellColData['donor'] <- metadata_common$unconfined_best_match_sample[match(rownames(proj_common@cellColData), metadata_common$archr_barcode)]

saveArchRProject(proj_common, 'Project_common_atac_rna_cells')