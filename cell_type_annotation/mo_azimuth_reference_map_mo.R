#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_azimuth_reference_map_mo.R
# Function: merge the seurat objects
############################################################################################################################

####################
# libraries        #
####################

library(SeuratData)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)

####################
# Functions        #
####################


#' do Azimuth reference mapping of query onto a reference
#' 
#' @param reference the reference Seurat object
#' @param query the query Seurat object
#' @returns a Seurat object
#' 
do_reference_mapping <- function(reference, query){
  # find transfer anchors
  anchors <- FindTransferAnchors(
    reference = reference,
    query = query,
    normalization.method = "LogNormalize",
    reference.reduction = "spca",
    dims = 1:50
  )
  # do the reference mapping
  query <- MapQuery(
    anchorset = anchors,
    query = query,
    reference = reference,
    refdata = list(
      mo_10x_cell_type = "seurat_annotations"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  return(query)
}


process_lane <- function(fragments_dir, seurat_objects_dir, lane, seurat_prepend='mo_', seurat_append='.rds', fragment_append='outs/atac_fragments.tsv.gz', h5_append='outs/filtered_feature_bc_matrix.h5') {
  # read the object
  tryCatch({
    object_loc <- paste(seurat_objects_dir, '/', seurat_prepend, lane, seurat_append, sep = '')
    # read RNA object
    object_lane <- readRDS(object_loc)
    
    # read the fragments
    fragments_loc <- paste(fragments_dir, '/', lane, '/', fragment_append, sep = '')
    # and counts
    h5_loc <- paste(fragments_dir, '/', lane, '/', h5_append, sep = '')
    counts <- Read10X_h5(filename = h5_loc)
    # subset to peaks
    counts <- counts$Peaks
    
    # create the atac metadata
    barcodes_short <- gsub('(-\\d+)', '', colnames(counts))
    barcodes_lane <- paste(barcodes_short, rep(lane, times = length(barcodes_short)), sep = '_')
    metadata <- data.frame(barcode=barcodes_short,
                           barcode_1=colnames(counts),
                           barcode_lane=barcodes_lane)
    metadata$lane <- lane
    metadata$batch <- lane
    rownames(metadata) <- metadata$barcode_1
    # change to a barcode unique across lanes
    colnames(counts) <- metadata$barcode_1
    
    # create the chromatin assay
    chrom_assay <- CreateChromatinAssay(
      counts = counts,
      sep = c(":", "-"),
      fragments = fragments_loc,
      min.cells = 10,
      min.features = 200
    )
    # rename
    chrom_assay <- RenameCells(chrom_assay, new.names = metadata[match(colnames(chrom_assay), metadata[['barcode_1']]), 'barcode_lane'])
    rownames(metadata) <- metadata$barcode_lane
    # and temporary object
    atac <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "peaks",
      meta.data = metadata[colnames(chrom_assay), ]
    )
    # now get the barcodes present in both
    barcodes_both <- intersect(atac@meta.data[['barcode_lane']], object_lane@meta.data[['barcode_lane']])
    
    # put all in one object
    cbmc <- CreateSeuratObject(counts = CreateAssayObject(object_lane[, barcodes_both]@assays$RNA$counts))
    # add this assay to the previously created Seurat object
    cbmc[["peaks"]] <- atac[, barcodes_both][['peaks']]
    # with metadata
    cbmc <- AddMetaData(cbmc, metadata = metadata)
    
    # try to do the peaks again
    DefaultAssay(cbmc) <- "peaks"
    # We exclude the first dimension as this is typically correlated with sequencing depth
    cbmc <- RunTFIDF(cbmc)
    cbmc <- FindVariableFeatures(cbmc, assay = 'peaks')
    cbmc <- RunSVD(cbmc)
    cbmc <- RunUMAP(cbmc, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_", return.model = T)
    cbmc <- FindNeighbors(cbmc, dims = 2:30, graph.name = 'peaks_snn', reduction = 'lsi')
    cbmc <- FindClusters(cbmc, resolution = 1.2, verbose = FALSE, graph.name = 'peaks_snn')
    
    # first do RNA
    DefaultAssay(cbmc) <- "RNA"
    # perform visualization and clustering steps
    cbmc <- NormalizeData(cbmc)
    cbmc <- FindVariableFeatures(cbmc)
    cbmc <- ScaleData(cbmc)
    cbmc <- RunPCA(cbmc, verbose = FALSE)
    cbmc <- FindNeighbors(cbmc, dims = 1:30, reduction = 'pca')
    cbmc <- FindClusters(cbmc, resolution = 1.2, verbose = FALSE)
    cbmc <- RunUMAP(cbmc, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", return.model = T)
    
    # do WNN
    cbmc <- FindMultiModalNeighbors(
      cbmc, reduction.list = list("pca", "lsi"), 
      dims.list = list(1:30, 2:30), modality.weight.name = c("RNA.weight", "ATAC.weigth")
    )
    cbmc <- RunUMAP(cbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = T)
    cbmc <- FindClusters(cbmc, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)
    # calculate spca 
    cbmc <- RunSPCA(cbmc, assay = 'RNA', graph = 'wsnn')
    # and cache neighbourhood index
    cbmc <- FindNeighbors(
      object = cbmc,
      reduction = "spca",
      dims = 1:50,
      graph.name = "spca.annoy.neighbors", 
      k.param = 50,
      cache.index = TRUE,
      return.neighbor = TRUE,
      l2.norm = TRUE
    )
    DefaultAssay(cbmc) <- 'RNA'
    return(cbmc)
  }, error = function(err) {
    print(paste('error in processing', err))
    return(NULL)
  })
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

lanes <- c('230105_lane1', '230105_lane2', '230105_lane3', '230105_lane4',
           '230105_lane5', '230105_lane6', '230105_lane7', '230105_lane8',
           '230112_lane1', '230112_lane2', '230112_lane3', '230112_lane4',
           '230112_lane5', '230112_lane6', '230112_lane7', '230112_lane8',
           '230120_lane1', '230120_lane2', '230120_lane3', '230120_lane4',
           '230120_lane5', '230120_lane6', '230120_lane7', '230120_lane8',
           '230127_lane1', '230127_lane2', '230127_lane3', '230127_lane4',
           '230127_lane5', '230127_lane6', '230127_lane7', '230127_lane8',
           '230202_lane1', '230202_lane2', '230202_lane3', '230202_lane4',
           '230202_lane5', '230202_lane6', '230202_lane7', '230202_lane8',
           '230209_lane1', '230209_lane2', '230209_lane3', '230209_lane4',
           '230209_lane5', '230209_lane6', '230209_lane7', '230209_lane8',
           '230216_lane1', '230216_lane2', '230216_lane3', '230216_lane4',
           '230216_lane5', '230216_lane6', '230216_lane7', '230216_lane8',
           '230223_lane1', '230223_lane2', '230223_lane3', '230223_lane4',
           '230223_lane5', '230223_lane6', '230223_lane7', '230223_lane8',
           '230302_lane1', '230302_lane2', '230302_lane3', '230302_lane4',
           '230302_lane5', '230302_lane6', '230302_lane7', '230302_lane8',
           '230316_lane1', '230316_lane2', '230316_lane3', '230316_lane4',
           '230316_lane5', '230316_lane6', '230316_lane7'
)

# location of the reference
reference_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cell_type_assignment/azimuth/10x_mo_reference.rds'
# read the reference
reference <- readRDS(reference_loc)

# location of fragments
fragments_dir <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/'
seurat_objects_dir <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/'

# do each lane
for (lane in lanes) {
  processed_object <- process_lane(fragments_dir = fragments_dir, seurat_objects_dir = seurat_objects_dir, lane = lane)
  # try to do mapping
  if (!is.null(processed_object)) {
    processed_object <- do_reference_mapping(reference = reference, query = processed_object)
    # save result
    result_loc <- paste(seurat_objects_dir, 'mo_', lane, 'multimodal_azi_mapped.rds', sep = '')
    saveRDS(processed_object, result_loc)
  }
}

# now summarize all
all_ct_predictions_per_lane <- list()
for (lane in lanes) {
  # load the object
  result_loc <- paste(seurat_objects_dir, 'mo_', lane, 'multimodal_azi_mapped.rds', sep = '')
  try({
    processed_object <- readRDS(result_loc)
    # extract the data we want
    annotation_lane <- data.frame(barcode = processed_object@meta.data[['barcode_lane']],
                                  predicted.mo_10x_cell_type = processed_object@meta.data[['predicted.mo_10x_cell_type']],
                                  predicted.mo_10x_cell_type.score = processed_object@meta.data[['predicted.mo_10x_cell_type.score']])
    # add to the list
    all_ct_predictions_per_lane[[lane]] <- annotation_lane
  })
}
# merge all
all_ct_predictions <- do.call('rbind', all_ct_predictions_per_lane)
# write the results
write.table(all_ct_predictions, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cell_type_assignment/azimuth/10x_multiome_PBMCs/mo_azimuth_ct_10xmultiome.tsv', sep = '\t', row.names = F, col.names = T)
