#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Martijn van der Werf, Roy Oelen
# Name: iterative-LSI.R
# Function: do LSI dimensional reduction, then perform clustering and further dimension reduction with UMAP
#
############################################################################################################################


setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)

addArchRThreads(threads=8)
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject('Project_all_cells')

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(resolution = c(2), sampleCells = 10000, maxClusters = 6, n.start
    = 10),
  firstSelection = "top",
  depthCol = "nFrags",
  varFeatures = 25000,
  dimsToUse = 1:30,
  LSIMethod = 2,
  scaleDims = TRUE,
  corCutOff = 0.75,
  binarize = TRUE,
  outlierQuantiles = c(0.02, 0.98),
  filterBias = TRUE,
  sampleCellsPre = 100000,
  projectCellsPre = FALSE,
  sampleCellsFinal = 100000,
  selectionMethod = "var",
  scaleTo = 10000,
  totalFeatures = 5e+05,
  filterQuantile = 0.995,
  excludeChr = c(),
  saveIterations = TRUE,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose =
    FALSE, fast_sgd = TRUE),
  nPlot = 10000,
  threads = getArchRThreads(),
  seed = 1,
  verbose = TRUE,
  force = TRUE,
  logFile = createLogFile("addIterativeLSI")
)

saveArchRProject(ArchRProj=proj, outputDirectory='Project_all_cells/')

#########################################
# IMPUTE CELL TYPES BASED ON LSI CLUSTERS
#########################################

# Add clusters
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    sampleCells=10000
)

proj_common <- addClusters(
    input = proj_common,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    sampleCells=100000
)
saveArchRProject(ArchRProj=proj, outputDirectory='Project_all_cells/')

# Add UMAP embedding

proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

proj_common <- addUMAP(
    ArchRProj = proj_common, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
# Plot embedding
umap_clusters <- plotEmbedding(ArchRProj = proj_common, colorBy = "cellColData", name = "celltype", embedding = "UMAP")

saveArchRProject(ArchRProj=proj, outputDirectory='Project_all_cells/')

plotPDF(umap_clusters, name = "Plot-UMAP-Sample-Imputed_celltype.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)