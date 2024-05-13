#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_multimodal_object_per_celltype.R
# Function: 
############################################################################################################################


####################
# libraries        #
####################

library(SeuratData)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)

####################
# Functions        #
####################

merge_atac_and_rna <- function(rna_object, atac_object) {
  # get joint barcodes
  joint_barcodes <- intersect(colnames(rna_object), colnames(atac_object))
  
  # ATAC analysis add gene annotation information
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  Annotation(atac_object) <- annotations
  
  # repeat QC steps performed in the WNN vignette
  rna_object <- rna_object[, joint_barcodes]
  atac_object <- atac_object[, joint_barcodes]
  
  # get the metadata for each modality
  rna_metadata <- rna_object@meta.data
  atac_metadata <- atac_object@meta.data
  # get the ones that are exclusive to the atac
  atac_only_colnames <- setdiff(colnames(atac_metadata), colnames(atac_metadata))
  # now only use the metadata that is new for a join
  both_metadata <- cbind(rna_metadata, atac_metadata[, atac_only_colnames])
  
  # We exclude the first dimension as this is typically correlated with sequencing depth
  atac_object <- RunTFIDF(atac_object)
  atac_object <- FindTopFeatures(atac_object, min.cutoff = "q0")
  
  # instead just get the chromatin assay data
  counts <- GetAssayData(atac_object, slot = "counts")
  # clear memory
  rm(atac_object)
  # create a new assay, while filtering out what we already removed in the RNA data
  chromatinassay <- CreateChromatinAssay(counts = counts, genome = "hg38")
  # clear memory
  rm(counts)
  
  # put all in one object
  both_object <- CreateSeuratObject(counts = CreateAssayObject(rna_object@assays$RNA$counts))
  # clear memory
  rm(rna_object)
  # add this assay to the previously created Seurat object
  both_object[["peaks"]] <- chromatinassay
  # clear memory
  rm(chromatinassay)
  
  # try to do the peaks again
  DefaultAssay(both_object) <- "peaks"
  # We exclude the first dimension as this is typically correlated with sequencing depth
  both_object <- RunTFIDF(both_object)
  both_object <- FindVariableFeatures(both_object, assay = 'peaks')
  both_object <- RunSVD(both_object)
  # both_object <- RunUMAP(both_object, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_", return.model = T)
  # both_object <- FindNeighbors(both_object, dims = 2:30, graph.name = 'peaks_snn', reduction = 'lsi')
  # both_object <- FindClusters(both_object, resolution = 1.2, verbose = FALSE, graph.name = 'peaks_snn')
  
  # then do RNA
  DefaultAssay(both_object) <- "RNA"
  # perform visualization and clustering steps
  both_object <- NormalizeData(both_object)
  both_object <- FindVariableFeatures(both_object)
  both_object <- ScaleData(both_object)
  both_object <- RunPCA(both_object, verbose = T)
  # both_object <- FindNeighbors(both_object, dims = 1:30)
  # both_object <- FindClusters(both_object, resolution = 1.2, verbose = FALSE)
  # both_object <- RunUMAP(both_object, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", return.model = T)
  
  # # do WNN
  # both_object <- FindMultiModalNeighbors(
  #   both_object, reduction.list = list("pca", "lsi"),
  #   dims.list = list(1:30, 2:30), modality.weight.name = c("RNA.weight", "ATAC.weigth")
  # )
  # both_object <- RunUMAP(both_object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = T)
  # both_object <- FindClusters(both_object, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
  # 
  # # calculate spca
  # both_object <- RunSPCA(both_object, assay = 'RNA', graph = 'wsnn')
  # # and cache neighbourhood index
  # both_object <- FindNeighbors(
  #   object = both_object,
  #   reduction = "spca",
  #   dims = 1:50,
  #   graph.name = "spca.annoy.neighbors",
  #   k.param = 50,
  #   cache.index = TRUE,
  #   return.neighbor = TRUE,
  #   l2.norm = TRUE
  # )
  DefaultAssay(both_object) <- 'RNA'
  
  # add back the metadata
  both_object <- AddMetaData(both_object, both_metadata[, colnames(both_metadata)])
  
  # return the result
  return(both_object)
}


combine_lanes <- function(lanes) {
  # we'll make a mapping for the lists
  lane_remapping <- list()
  # we'll check each lane
  for (lane in lanes) {
    # we'll extract the last number from the lane
    lane_nr <- regmatches(lane,regexpr("\\d+$",lane))
    # check if it is even
    if (as.numeric(lane_nr) %% 2 == 0) {
      # remove the last number
      lane_no_nr <- substr(lane, 1, nchar(lane) -1)
      # get the previous number
      previous_lane_nr <- as.numeric(lane_nr) - 1
      # get the combined lane
      combined_lane <- paste(lane_no_nr, previous_lane_nr, lane_nr, sep = '')
      # put in list
      lane_remapping[[lane]] <- combined_lane
    }
    else if (as.numeric(lane_nr) %% 2 == 1) {
      # just add the next number
      combined_lane <- paste(lane, as.numeric(lane_nr) + 1, sep = '')
      # put in list
      lane_remapping[[lane]] <- combined_lane
    }
  }
  return(lane_remapping)
}


#' add the inflammation assignments  to the Seurat object
#' 
#' @param seurat_object The Seurat object to add the inflammation status to
#' @param sample_sheet The sample sheet containing lanes, participants and inflammation statuses
#' @param seurat_lane_column The column in the Seurat metadata denoting the 10x lane
#' @param sheet_lane_column The column in the sample sheet denoting the 10x lane
#' @param seurat_participant_column The column in the Seurat metadata denoting the participant assignment
#' @param sheet_participants_column The column in the sample sheet containing the participants per lane
#' @param seurat_inflammation_column The column in the Seurat metadata to add the inflammation status in
#' @param sheet_inflammation_column The column in the sample sheet containing the inflammation statuses per lane
#' @returns the Seurat object with the inflammation status added
#' lpmcv2 <- add_inflammation_status(lpmcv2, sample_sheet)
add_inflammation_status <- function(seurat_object, sample_sheet, seurat_lane_column='lane', sheet_lane_column='lane', seurat_participant_column='soup_final_sample_assignment', sheet_participants_column='genoid', seurat_inflammation_column='inflammation_status', sheet_inflammation_column='inflammation_status') {
  # create a mapping of lane+sample to inflammation status
  mapping_per_lane_list <- list()
  for (i in 1:nrow(sample_sheet)) {
    # extract lane
    lane <- sample_sheet[i, sheet_lane_column]
    # extract the participants
    participant <- sample_sheet[i, sheet_participants_column]
    # and the inflammation condition
    condition <- sample_sheet[i, sheet_inflammation_column]
    # if not set, set to unknown
    if (is.null(condition)) {
      condition <- 'unknown'
    }
    # subset to the barcodes which have are this lane and participant
    barcodes_match <- rownames(seurat_object@meta.data[!is.na(seurat_object@meta.data[[seurat_lane_column]]) &
                                                         seurat_object@meta.data[[seurat_lane_column]] == lane &
                                                         !is.na(seurat_object@meta.data[[seurat_participant_column]]) &
                                                         seurat_object@meta.data[[seurat_participant_column]] == participant, ])
    # only add if there are matching barcodes
    if (length(barcodes_match) > 0) {
      # create dataframe
      df_lane_part <- data.frame(barcode = barcodes_match, condition = rep(condition, times = length(barcodes_match)))
      # set the colname to be the one we chose
      colnames(df_lane_part) <- c('barcode', seurat_inflammation_column)
      # then add to the list
      mapping_per_lane_list[[paste(lane, participant, sep = ':')]] <- df_lane_part
    }
  }
  # now merge all together
  mapping_all <- do.call('rbind', mapping_per_lane_list)
  # set the barcode as rownames
  rownames(mapping_all) <- mapping_all[['barcode']]
  # finally add to the object
  seurat_object <- AddMetaData(seurat_object, mapping_all[seurat_inflammation_column])
  return(seurat_object)
}


read_barcode_and_lane <- function(seurat_object) {
  # do the split first
  seurat_object_rowsnames_split <- strsplit(colnames(seurat_object), split = '_')
  # now do a list apply
  df_per_barcode <- lapply(seurat_object_rowsnames_split, FUN = function(x) {
    data.frame(lane = paste(x[2], x[3], sep = '_'), barcode = x[1])
  })
  # merge all of them together
  extra_metadata <- do.call('rbind', df_per_barcode)
  # set the rownames to be the original ones
  rownames(extra_metadata) <- colnames(seurat_object)
  # now add the extra data we have
  seurat_object <- AddMetaData(seurat_object, extra_metadata)
  return(seurat_object)
}

#' add the inflammation assignments  to the Seurat object
#' 
#' @param seurat_object The Seurat object to add the inflammation status to
#' @param sample_sheet The sample sheet containing lanes, participants and inflammation statuses
#' @param seurat_lane_column The column in the Seurat metadata denoting the 10x lane
#' @param sheet_lane_column The column in the sample sheet denoting the 10x lane
#' @param seurat_participant_column The column in the Seurat metadata denoting the participant assignment
#' @param sheet_participants_column The column in the sample sheet containing the participants per lane
#' @param seurat_inflammation_column The column in the Seurat metadata to add the inflammation status in
#' @param sheet_inflammation_column The column in the sample sheet containing the inflammation statuses per lane
#' @returns the Seurat object with the inflammation status added
#' lpmcv2 <- add_inflammation_status(lpmcv2, sample_sheet)
add_inflammation_status <- function(seurat_object, sample_sheet, seurat_lane_column='lane', sheet_lane_column='lane', seurat_participant_column='soup_final_sample_assignment', sheet_participants_column='genoid', seurat_inflammation_column='inflammation_status', sheet_inflammation_column='inflammation_status') {
  # create a mapping of lane+sample to inflammation status
  mapping_per_lane_list <- list()
  for (i in 1:nrow(sample_sheet)) {
    # extract lane
    lane <- sample_sheet[i, sheet_lane_column]
    # extract the participants
    participant <- sample_sheet[i, sheet_participants_column]
    # and the inflammation condition
    condition <- sample_sheet[i, sheet_inflammation_column]
    # if not set, set to unknown
    if (is.null(condition)) {
      condition <- 'unknown'
    }
    # subset to the barcodes which have are this lane and participant
    barcodes_match <- rownames(seurat_object@meta.data[!is.na(seurat_object@meta.data[[seurat_lane_column]]) &
                                                         seurat_object@meta.data[[seurat_lane_column]] == lane &
                                                         !is.na(seurat_object@meta.data[[seurat_participant_column]]) &
                                                         seurat_object@meta.data[[seurat_participant_column]] == participant, ])
    # only add if there are matching barcodes
    if (length(barcodes_match) > 0) {
      # create dataframe
      df_lane_part <- data.frame(barcode = barcodes_match, condition = rep(condition, times = length(barcodes_match)))
      # set the colname to be the one we chose
      colnames(df_lane_part) <- c('barcode', seurat_inflammation_column)
      # then add to the list
      mapping_per_lane_list[[paste(lane, participant, sep = ':')]] <- df_lane_part
    }
  }
  # now merge all together
  mapping_all <- do.call('rbind', mapping_per_lane_list)
  # set the barcode as rownames
  rownames(mapping_all) <- mapping_all[['barcode']]
  # finally add to the object
  seurat_object <- AddMetaData(seurat_object, mapping_all[seurat_inflammation_column])
  return(seurat_object)
}


add_inflammation_status_each_object <- function(seurat_object_list, sample_sheet, seurat_lane_column='lane', sheet_lane_column='lane', seurat_participant_column='soup_final_sample_assignment', sheet_participants_column='genoid', seurat_inflammation_column='inflammation_status', sheet_inflammation_column='inflammation_status') {
  # let's go over each object
  for (object_name in names(seurat_object_list)) {
    # do the condition assignment
    seurat_object_list[[object_name]] <- add_inflammation_status(
      seurat_object_list[[object_name]],
      sample_sheet=sample_sheet, 
      seurat_lane_column=seurat_lane_column,
      sheet_lane_column=sheet_lane_column, 
      seurat_participant_column=seurat_participant_column, 
      sheet_participants_column=sheet_participants_column, 
      seurat_inflammation_column=seurat_inflammation_column, 
      sheet_inflammation_column=sheet_inflammation_column
    )
  }
  return(seurat_object_list)
}



####################
# Options          #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory, we'll take 95 percent of the max memory in GB
mem_gb_requested <- 1900
options(future.globals.maxSize = (0.95 * mem_gb_requested) * 1000 * 1024^2)

# set the max parallel threads
#register(MulticoreParam(4))

# set seed
set.seed(7777)

####################
# Main Code        #
####################

# location of the condition assignment
condition_assignment_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_monocyte_based_condition_numbers.tsv'
# read the conditions
condition_assignments <- read.table(condition_assignment_loc, header = T, sep = '\t')

# location of the cell type objects
cell_type_objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltypemajor_1_80.rds'

# read the object
cell_type_objects <- readRDS(cell_type_objects_loc)


# add barcodes back
for(cell_type in names(cell_type_objects)) {
  cell_type_objects[[cell_type]] <- read_barcode_and_lane(cell_type_objects[[cell_type]])
}
# get the assignment matrices
correlation_mapping_per_barcode_all <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched_vs_all.tsv', header = T, sep = '\t')
# set barcodes and remove data we already have
rownames(correlation_mapping_per_barcode_all) <- correlation_mapping_per_barcode_all[['barcode_lane']]
correlation_mapping_per_barcode_all[, c('lane', 'barcode_lane', 'barcode', 'barcode_original')] <- NULL
# now let's get the souporcell data specifically, which would be the same for all and per-lane
soup_only <- correlation_mapping_per_barcode_all[, setdiff(colnames(correlation_mapping_per_barcode_all), c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation'))]
# and the correlation data
correlations_unconfined <- correlation_mapping_per_barcode_all[, c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')]
# add the confined sample
for(cell_type in names(cell_type_objects)) {
  cell_type_objects[[cell_type]] <- AddMetaData(cell_type_objects[[cell_type]], correlations_unconfined[, colnames(correlations_unconfined)])
}

# add the conditions of the original sheet
cell_type_objects <- add_inflammation_status_each_object(cell_type_objects, condition_assignments, seurat_participant_column='best_match_sample', sheet_participants_column = 'sample_final', seurat_inflammation_column = 'inflammation_sheet', sheet_inflammation_column = 'condition')
cell_type_objects <- add_inflammation_status_each_object(cell_type_objects, condition_assignments, seurat_participant_column='best_match_sample', sheet_participants_column = 'sample', seurat_inflammation_column = 'inflammation_prev', sheet_inflammation_column = 'cond_prev')
# rename CA in the original assignment, and set a consensus one
for (cell_type in names(cell_type_objects)) {
  cell_type_object <- cell_type_objects[[cell_type]]
  cell_type_object@meta.data[!is.na(cell_type_object@meta.data[['inflammation_sheet']]) &
                               cell_type_object@meta.data[['inflammation_sheet']] == '24hCa', 'inflammation_sheet'] <- '24hCA'
  # now also set the final inflammation assignment
  cell_type_object@meta.data[['inflammation_final']] <- cell_type_object@meta.data[['inflammation_sheet']]
  cell_type_object@meta.data[is.na(cell_type_object@meta.data[['inflammation_final']]), 'inflammation_final'] <- cell_type_object@meta.data[is.na(cell_type_object@meta.data[['inflammation_final']]), 'inflammation_prev']
  # add back to list
  cell_type_objects[[cell_type]] <- cell_type_object
}

# locations of objects
objects_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
mo_object_loc <- paste(objects_loc, 'mo_all_20240223_seuratv5_normalized.rds', sep = '')
# read the Seurat object
seurat_object <- readRDS(mo_object_loc)
# get the lane remapping
lane_remapping <- combine_lanes(unique(seurat_object@meta.data$lane))
# add to the object
seurat_object@meta.data[['lane_both']] <- as.vector(unlist(lane_remapping[seurat_object@meta.data[['lane']]]))
# do the condition assignment
seurat_object <- add_inflammation_status(
  seurat_object,
  sample_sheet=condition_assignments, 
  seurat_lane_column='lane',
  sheet_lane_column='lane', 
  seurat_participant_column='sample_final', 
  sheet_participants_column='sample_final', 
  seurat_inflammation_column='inflammation_sheet', 
  sheet_inflammation_column='condition'
)
seurat_object <- add_inflammation_status(
  seurat_object,
  sample_sheet=condition_assignments, 
  seurat_lane_column='lane',
  sheet_lane_column='lane', 
  seurat_participant_column='sample_final', 
  sheet_participants_column='sample', 
  seurat_inflammation_column='inflammation_prev', 
  sheet_inflammation_column='cond_prev'
)
# rename 24hCA
seurat_object@meta.data[!is.na(seurat_object@meta.data[['inflammation_sheet']]) &
                          seurat_object@meta.data[['inflammation_sheet']] == '24hCa', 'inflammation_sheet'] <- '24hCA'
# now also set the final inflammation assignment
seurat_object@meta.data[['inflammation_final']] <- seurat_object@meta.data[['inflammation_sheet']]
seurat_object@meta.data[is.na(seurat_object@meta.data[['inflammation_final']]), 'inflammation_final'] <- seurat_object@meta.data[is.na(seurat_object@meta.data[['inflammation_final']]), 'inflammation_prev']
# save result
mo_annotated_loc <- paste(objects_loc, 'mo_all_20240513_seuratv5_annotated.rds', sep = '')
saveRDS(seurat_object, mo_annotated_loc)

# get multimodal monocytes
mono_rna <- seurat_object[, !is.na(seurat_object@meta.data$celltype_imputed_lowerres) & seurat_object@meta.data$celltype_imputed_lowerres == 'monocyte']
mono_atac <- cell_type_objects[['monocyte']]
#rm(seurat_object)
#rm(cell_type_objects)
monocyte_multimodal <- merge_atac_and_rna(
  rna_object = mono_rna,
  atac_object = mono_atac
)
saveRDS(monocyte_multimodal, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_monocyte_1_80_20240513.rds')
rm(monocyte_multimodal)

# multimodal CD4T
cd4t_rna <- seurat_object[, !is.na(seurat_object@meta.data$celltype_imputed_lowerres) & seurat_object@meta.data$celltype_imputed_lowerres == 'CD4T']
cd4t_atac <- cell_type_objects[['CD4T']]
#rm(seurat_object)
#rm(cell_type_objects)
cd4t_multimodal <- merge_atac_and_rna(
  rna_object = cd4t_rna,
  atac_object = cd4t_atac
)
saveRDS(cd4t_multimodal, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_cd4t_1_80_20240513.rds')
rm(cd4t_multimodal)

# multimodal CD8T
cd8t_rna <- seurat_object[, !is.na(seurat_object@meta.data$celltype_imputed_lowerres) & seurat_object@meta.data$celltype_imputed_lowerres == 'CD8T']
cd8t_atac <- cell_type_objects[['CD8T']]
#rm(seurat_object)
#rm(cell_type_objects)
cd8t_multimodal <- merge_atac_and_rna(
  rna_object = cd8t_rna,
  atac_object = cd8t_atac
)
saveRDS(cd8t_multimodal, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_cd8t_1_80_20240513.rds')
rm(cd8t_multimodal)

# multimodal DC
dc_rna <- seurat_object[, !is.na(seurat_object@meta.data$celltype_imputed_lowerres) & seurat_object@meta.data$celltype_imputed_lowerres == 'DC']
dc_atac <- cell_type_objects[['DC']]
#rm(seurat_object)
#rm(cell_type_objects)
dc_multimodal <- merge_atac_and_rna(
  rna_object = dc_rna,
  atac_object = dc_atac
)
saveRDS(dc_multimodal, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_dc_1_80_20240513.rds')
rm(dc_multimodal)

# multimodal NK
nk_rna <- seurat_object[, !is.na(seurat_object@meta.data$celltype_imputed_lowerres) & seurat_object@meta.data$celltype_imputed_lowerres == 'NK']
nk_atac <- cell_type_objects[['NK']]
#rm(seurat_object)
#rm(cell_type_objects)
nk_multimodal <- merge_atac_and_rna(
  rna_object = nk_rna,
  atac_object = nk_atac
)
saveRDS(nk_multimodal, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_nk_1_80_20240513.rds')

# multimodal B
b_rna <- seurat_object[, !is.na(seurat_object@meta.data$celltype_imputed_lowerres) & seurat_object@meta.data$celltype_imputed_lowerres == 'B']
b_atac <- cell_type_objects[['B']]
#rm(seurat_object)
#rm(cell_type_objects)
b_multimodal <- merge_atac_and_rna(
  rna_object = b_rna,
  atac_object = b_atac
)
saveRDS(b_multimodal, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_b_1_80_20240513.rds')
rm(b_multimodal)
