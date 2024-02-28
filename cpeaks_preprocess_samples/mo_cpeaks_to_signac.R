#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cpeaks_to_signac.R
# Function: perform normalization
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Signac)
# this database is needed for annotations
library(EnsDb.Hsapiens.v86)


####################
# Functions        #
####################


read_cpeaks <- function(cpeaks_loc, cellranger_loc, lane, annotations) {
  # read the features
  features <- read.table(paste(cpeaks_loc, '/', lane, '/features.tsv.gz', sep = ''), sep = '\t', header = F)
  features[['peak']] <- paste(features[[1]], ':',  features[[2]], '-', features[[3]], sep = '')
  write.table(features, gzfile(paste(cpeaks_loc, '/', lane, '/named_features.tsv.gz', sep = '')), sep = '\t', row.names = F, col.names = F, quote = F)
  # read the actual counts
  peaks <- ReadMtx(mtx = paste(cpeaks_loc, '/', lane, '/atac_matrix.mtx.gz', sep = ''),
                    cells = paste(cpeaks_loc, '/', lane, '/barcodes.tsv.gz', sep = ''),
                    features = paste(cpeaks_loc, '/', lane, '/named_features.tsv.gz', sep = ''),
                    feature.column = 4)
  
  # create some metadata, for now, we'll first just store the lane here
  barcodes_short <- gsub('(-\\d+)', '', colnames(peaks))
  barcodes_lane <- paste(barcodes_short, rep(lane, times = length(barcodes_short)), sep = '_')
  # create metadata
  metadata <- data.frame(barcode=barcodes_short,
                         barcode_1=colnames(peaks),
                         barcode_lane=barcodes_lane)
  # save the lane the samples are from
  metadata$lane <- lane
  metadata$batch <- lane
  #rownames(metadata) <- metadata$barcode_lane
  # change to a barcode unique across lanes
  #colnames(peaks) <- metadata$barcode_lane
  # get the fragments
  fragments_loc <- paste(cellranger_loc, '/', lane, '/outs/atac_fragments.tsv.gz', sep = '')
  # create the chromatin assay
  chrom_assay <- CreateChromatinAssay(
    counts = peaks,
    sep = c(":", "-"),
    fragments = fragments_loc,
    min.cells = 10,
    min.features = 200
  )
  # create object
  seurat_object <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata,
    project = 'wijst_multiome'
  )
  # set the annotations to the object now
  Annotation(seurat_object) <- annotations
  # rename to lane and barcode
  seurat_object <- RenameCells(seurat_object, new.names = seurat_object@meta.data[['barcode_lane']])
  # return result
  return(seurat_object)
}

####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
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

# these are the lanes
lanes <- c('230105_lane1', '230105_lane2', '230105_lane3', '230105_lane4',
           '230105_lane5', '230105_lane6', '230105_lane7', '230105_lane8'
           # '230112_lane1', '230112_lane2', '230112_lane3', '230112_lane4',
           # '230112_lane5', '230112_lane6', '230112_lane7', '230112_lane8',
           # '230120_lane1', '230120_lane2', '230120_lane3', '230120_lane4',
           # '230120_lane5', '230120_lane6', '230120_lane7', '230120_lane8',
           # '230127_lane1', '230127_lane2', '230127_lane3', '230127_lane4',
           # '230127_lane5', '230127_lane6', '230127_lane7', '230127_lane8',
           # '230202_lane1', '230202_lane2', '230202_lane3', '230202_lane4',
           # '230202_lane5', '230202_lane6', '230202_lane7', '230202_lane8',
           # '230209_lane1', '230209_lane2', '230209_lane3', '230209_lane4',
           # '230209_lane5', '230209_lane6', '230209_lane7', '230209_lane8',
           # '230216_lane1', '230216_lane2', '230216_lane3', '230216_lane4',
           # '230216_lane5', '230216_lane6', '230216_lane7', '230216_lane8',
           # '230223_lane1', '230223_lane2', '230223_lane3', '230223_lane4',
           # '230223_lane5', '230223_lane6', '230223_lane7', '230223_lane8',
           # '230302_lane1', '230302_lane2', '230302_lane3', '230302_lane4',
           # '230302_lane5', '230302_lane6', '230302_lane7', '230302_lane8',
           # '230316_lane1', '230316_lane2', '230316_lane3', '230316_lane4',
           # '230316_lane5', '230316_lane6', '230316_lane7', '230316_lane8'
)

# location of the cpeaks objects
cpeaks_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/output/default/'
# get cellranger loc
cellranger_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/joint/alignment/b38/'
# location of the gene annotations
gtf_loc <- '/groups/umcg-franke-scrna/tmp03/external_datasets/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz'

# check each lane
for (lane in lanes) {
  # read the object
  object_lane <- read_cpeaks(cpeaks_loc, cellranger_loc, lane = lane, annotations = annotations)
  # save the result
  cpeaks_object_loc <- paste(cpeaks_loc, '/', lane, '/', 'mo_cpeaks_', lane, '.rds', sep = '')
  saveRDS(object_lane, cpeaks_object_loc)
  rm(object_lane)
}
