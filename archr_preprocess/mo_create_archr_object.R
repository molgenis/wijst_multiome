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
#devtools::install_version("Matrix", version = "1.6-1", repos = "http://cran.us.r-project.org")
#library(Matrix)

####################
# Functions        #
####################

#' get the number eGenes per cell type from EMP output
#' 
#' @param freemuxlet_output_loc base location of the demuxlet output
#' @param lanes which lanes to analyse (optional, if omitted, will do a list.files to discover lanes)
#' @param freemux_prepend what the output file per lane is prepended with (optional, empty string by default)
#' @param freemux_append what the output file per lane is appended with (optional, '.best')
#' @returns a dataframe with the demuxlet output of all the lanes, combined with a lane and lane+barcode column added
#' 
get_freemux_assignments <- function(freemuxlet_output_loc, lanes=NULL, freemux_prepend='', freemux_append='_mo_freemux_atac_cluster_correlations.tsv') {
  # first put each result in a list of tables
  result_per_lane <- list()
  # get which lanes to check
  if (is.null(lanes)) {
    # if not supplied we will try to use regex
    lane_files <- list.files(path = freemuxlet_output_loc, pattern = paste(freemux_prepend, '*', freemux_append, sep = ''), full.names = F)
    lanes <- sub(pattern = freemux_prepend, replacement = "", lane_files)
    lanes <- sub(pattern = freemux_append, replacement = "", lane_files)
  }
  # check each lane
  for (lane in lanes) {
    # get the full path
    full_freemux_path <- paste(freemuxlet_output_loc, '/', freemux_prepend, lane, freemux_append, sep = '')
    # check if the file exists
    if (file.exists(full_freemux_path)) {
      # read the file
      freemux_lane <- read.table(full_freemux_path, header = T, sep = '\t', row.names = 1)
      # get the bare barcode and the barcode+lane
      barcodes_short <- gsub('(-\\d+)', '', freemux_lane$BARCODE)
      barcodes_lane <- paste(barcodes_short, rep(lane, times = length(barcodes_short)), sep = '_')
      # get the final assignment as well
      final_assignments <- as.character(freemux_lane$cluster)
      # but overwrite where it should be a doublet
      final_assignments[!is.na(freemux_lane$DROPLET.TYPE) & freemux_lane$DROPLET.TYPE == 'DBL'] <- 'multiplet'
      final_assignments[!is.na(freemux_lane$DROPLET.TYPE) & freemux_lane$DROPLET.TYPE == 'AMB'] <- 'ambiguous'
      # put the new variables in a dataframe
      freemux_additional <- data.frame(
        lane = rep(lane, times = length(barcodes_short)),
        barcode_1 = freemux_lane$BARCODE,
        barcode_bare = barcodes_short,
        barcode_lane = barcodes_lane,
        assignment = final_assignments
      )
      # now add this to the existing freemuxlet output
      freemux_lane <- cbind(freemux_additional, freemux_lane)
      result_per_lane[[lane]] <- freemux_lane
    }
    else{
      message(paste('lane ', lane, ' does not exist at location ', full_freemux_path, ', skipped', sep = ''))
    }
  }
  # combine the outputs
  all_freemux <- do.call('rbind', result_per_lane)
  return(all_freemux)
}


add_conditions_demuxed <- function(demuxing_output, condition_mapping, lane_column_mapping='lane', sample_column_mapping='sample', condition_column_mapping='condition', lane_column_metadata='lane', sample_column_metadata='best_match_sample', condition_column_metadata='condition') {
  # add the metadata column
  demuxing_output[[condition_column_metadata]] <- NA
  # check each row
  for (row_i in 1 : nrow(condition_mapping)) {
    # get the lane
    lane <- condition_mapping[row_i, lane_column_mapping]
    sample <- condition_mapping[row_i, sample_column_mapping]
    condition <- condition_mapping[row_i, condition_column_mapping]
    if (!is.na(lane) & !is.na(sample)) {
      # set the condition
      demuxing_output[!is.na(demuxing_output[[lane_column_metadata]]) & demuxing_output[[lane_column_metadata]] == lane &
                                !is.na(demuxing_output[[sample_column_metadata]]) & demuxing_output[[sample_column_metadata]] == sample
                              , condition_column_metadata] <- condition
    }
  }
  return(demuxing_output)
}


####################
# Settings         #
####################

# this is the reference genome
addArchRGenome('hg38')
# the number of threads to use
addArchRThreads(threads = 4) 
# we need some more memory
options(future.globals.maxSize = 190 * 1000 * 1024^2)
# set seed
set.seed(7777)
# set locking
#addArchRLocking(locking = TRUE)

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
  addGeneScoreMat = TRUE,
  threads = 1
)

write.table(ArrowFiles, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_arrowfiles.txt', row.names = F, col.names = F)
ArrowFiles <- read.table('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_arrowfiles.txt', header = F)$V1

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

saveRDS(mo_peaks, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered.rds')


#mo_peaks <- readRDS('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered.rds')
# read the metadata we have based on the RNA-seq
rna_metadata_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv'
rna_metadata <- read.table(rna_metadata_loc, header = T, sep = '\t')
# add column as it is in the ArchR data
rna_metadata[['archr_barcode']] <- paste(rna_metadata[['lane']], rna_metadata[['barcode_1']], sep = '#')
# add some metadata we got from the rna metadata
for (metadata_col in c('lane', 'soup_best_match_sample', 'soup_best_match_correlation', 'cell_type', 'cell_type_score', 'cell_type_lowerres', 'cell_type_lowerres_imputed', 'condition', 'condition_imputed', 'seurat_clusters')) {
  mo_peaks <- addCellColData(ArchRProj = mo_peaks, data = rna_metadata[[metadata_col]],
                             cells = rna_metadata[['archr_barcode']], name = metadata_col)
}

# add the ATAC metadata
freemuxlet_output_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/demultiplexing/freemuxlet/atac/assignments/'
freemux_prepend <- ''
freemux_append <- '_mo_freemux_atac_cluster_correlations.tsv'
freemux_assignments <- get_freemux_assignments(freemuxlet_output_loc, freemux_prepend = '', freemux_append = freemux_append)
# get the condition mapping
condition_mapping_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'
condition_mapping <- read.table(condition_mapping_loc, header = T, sep = '\t')
freemux_assignments <- add_conditions_demuxed(freemux_assignments, condition_mapping)
# add column as it is in the ArchR data
freemux_assignments[['archr_barcode']] <- paste(freemux_assignments[['lane']], freemux_assignments[['barcode_1']], sep = '#')
for (metadata_col in c('condition')) {
  mo_peaks <- addCellColData(ArchRProj = mo_peaks, data = freemux_assignments[[metadata_col]],
                             cells = freemux_assignments[['archr_barcode']], name = paste('freemux', metadata_col, sep = '_'))
}

# ncell 1307396
# filter the doublets
mo_peaks <- filterDoublets(mo_peaks)
# ncell 1262102

# get the cells that are not doublets
non_doublets <- rna_metadata[!is.na(rna_metadata[['soup_best_match_sample']]), 'archr_barcode']
# filter ones that we have in the peaks
non_doublets <- non_doublets[non_doublets %in% mo_peaks$cellNames]


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
  dimsToUse = 1:30,
  sampleCellsFinal = 50000,
  projectCellsPre = T
)
saveRDS(mo_peaks, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered_lsi.rds')

#mo_peaks <- readRDS('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered_lsi.rds')
# and do dimensional reduction using UMAP
mo_peaks <- addUMAP(
  ArchRProj = mo_peaks, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)
saveRDS(mo_peaks, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered_umap.rds')

# then do clustering
mo_peaks <- addClusters(
  input = mo_peaks,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "archr_clusters",
  resolution = 0.8,
  sampleCells = 25000
)
saveRDS(mo_peaks, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered_clus.rds')

# do scran clustering
mo_peaks <- addClusters(
  input = mo_peaks,
  reducedDims = "IterativeLSI",
  method = "scran",
  name = "scran_clusters",
  k = 15,
  sampleCells = 25000
)
saveRDS(mo_peaks, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered_clus_scran.rds')


# get UT
mo_peaks <-  readRDS('/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/archr_preprocess_samples/objects/merged/mo_archr_unfiltered.rds')
for (metadata_col in c('lane', 'soup_best_match_sample', 'soup_best_match_correlation', 'cell_type', 'cell_type_score', 'cell_type_lowerres', 'cell_type_lowerres_imputed', 'condition', 'condition_imputed', 'seurat_clusters')) {
  mo_peaks <- addCellColData(ArchRProj = mo_peaks, data = rna_metadata[[metadata_col]],
                             cells = rna_metadata[['archr_barcode']], name = metadata_col)
}
for (metadata_col in c('condition')) {
  mo_peaks <- addCellColData(ArchRProj = mo_peaks, data = freemux_assignments[[metadata_col]],
                             cells = freemux_assignments[['archr_barcode']], name = paste('freemux', metadata_col, sep = '_'))
}
ut_cells <- freemux_assignments[!is.na(freemux_assignments[['condition']]) & freemux_assignments[['condition']] == 'UT' & !is.na(freemux_assignments[['DROPLET.TYPE']]) & freemux_assignments[['DROPLET.TYPE']] == 'SNG', 'archr_barcode']
ut_cells <- ut_cells[!is.na(ut_cells) & ut_cells %in% mo_peaks$cellNames]
mo_peaks$cellNames[!is.na(mo_peaks$freemux_condition) & mo_peaks$freemux_condition == 'UT' & !is.na(mo_peaks$freemux_DROPLET.TYPE) & mo_peaks$freemux_DROPLET.TYPE == 'SNG']
ut_peaks <- subsetArchRProject(mo_peaks, cells = ut_cells, outputDirectory = '/local/2339855/projects/multiome/ongoing/archr_preprocess_samples/objects/UT/', force = T, threads = 1)
# # reread all the files
# ut_peaks <- ArchRProject(
#   ArrowFiles = list.files('/local/2339855/projects/multiome/ongoing/archr_preprocess_samples/objects/UT/ArrowFiles/', full.names = T, pattern = '*.arrow'), 
#   outputDirectory = '/local/2339855/projects/multiome/ongoing/archr_preprocess_samples/objects/UT/',
#   copyArrows = FALSE # so that if you modify the Arrow files, you have an original copy for later usage
# )
# for (metadata_col in c('lane', 'soup_best_match_sample', 'soup_best_match_correlation', 'cell_type', 'cell_type_score', 'cell_type_lowerres', 'cell_type_lowerres_imputed', 'condition', 'condition_imputed', 'seurat_clusters')) {
#   ut_peaks <- addCellColData(ArchRProj = ut_peaks, data = rna_metadata[[metadata_col]],
#                              cells = rna_metadata[['archr_barcode']], name = metadata_col)
# }
# for (metadata_col in c('condition')) {
#   ut_peaks <- addCellColData(ArchRProj = ut_peaks, data = freemux_assignments[[metadata_col]],
#                              cells = freemux_assignments[['archr_barcode']], name = paste('freemux', metadata_col, sep = '_'))
# }
# # and bad subset
# ut_peaks <- ut_peaks[!is.na(ut_peaks$freemux_condition) & ut_peaks$freemux_condition == 'UT' & !is.na(ut_peaks$freemux_DROPLET.TYPE) & ut_peaks$freemux_DROPLET.TYPE == 'SNG', ]
# # do the LSI in lieu of PCA
# ut_peaks <- addIterativeLSI(
#   ArchRProj = ut_peaks,
#   useMatrix = "TileMatrix", 
#   name = "IterativeLSI", 
#   iterations = 2, 
#   clusterParams = list(
#     resolution = c(1.2), 
#     sampleCells = 10000, 
#     n.start = 10
#   ), 
#   varFeatures = 25000, 
#   dimsToUse = 1:30,
#   sampleCellsFinal = 50000
# )
# # then do clustering
# ut_peaks <- addClusters(
#   input = ut_peaks,
#   reducedDims = "IterativeLSI",
#   method = "Seurat",
#   name = "archr_clusters",
#   resolution = 0.8,
#   sampleCells = 25000
# )
# # and do dimensional reduction using UMAP
# ut_peaks <- addUMAP(
#   ArchRProj = ut_peaks, 
#   reducedDims = "IterativeLSI", 
#   name = "UMAP", 
#   nNeighbors = 30, 
#   minDist = 0.5, 
#   metric = "cosine"
# )
# 
