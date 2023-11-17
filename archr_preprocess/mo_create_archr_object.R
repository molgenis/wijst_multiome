#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_archr_object.R
# Function: 
############################################################################################################################


####################
# libraries        #
####################

library(ArchR)


####################
# Functions        #
####################


####################
# Settings         #
####################

# this is the reference genome
addArchRGenome('hg38')
# the number of threads to use
addArchRThreads(threads = 4) 


####################
# Main Code        #
####################

# these are the lanes
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

# location of cellranger output
cellranger_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/processed/joint/alignment/b38/'
cellranger_append <- '/outs/atac_fragments.tsv.gz'

# locations of the inputs
inputs_lanes <- paste(cellranger_loc, lanes, cellranger_append, sep = '')

# create everything
ArrowFiles <- createArrowFiles(
  inputFiles = inputs_lanes,
  sampleNames = lanes,
  filterTSS = 4, # don't set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# add doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, # how many cells near a "pseudo-doublet" to count.
  knnMethod = 'UMAP', # embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

# create the project
mo_peaks <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/',
  copyArrows = TRUE # so that if you modify the Arrow files, you have an original copy for later usage
)

# read the metadata we have based on the RNA-seq
rna_metadata_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv'
rna_metadata <- read.table(rna_metadata_loc, header = T, sep = '\t')
# add column as it is in the ArchR data
rna_metadata[['archr_barcode']] <- paste(rna_metadata[['lane']], rna_metadata[['barcode_1']], sep = '#')

# add some metadata we got from the rna metadata
for (metadata_col in c('lane', 'soup_best_match_sample', 'soup_best_match_correlation', 'cell_type', 'cell_type_score', 'cell_type_lowerres', 'cell_type_lowerres_imputed', 'condition', 'condition_imputed')) {
  mo_peaks <- addCellColData(ArchRProj = mo_peaks, data = rna_metadata[[metadata_col]],
                             cells = rna_metadata[['archr_barcode']], name = metadata_col)
}

# filter the doublets
mo_peaks <- filterDoublets(mo_peaks)

# do the LSI in lieu of PCA
mo_peaks <- addIterativeLSI(
  ArchRProj = mo_peaks,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(1.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

# then do clustering
mo_peaks <- addClusters(
  input = mo_peaks,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "archr_clusters",
  resolution = 1.2
)