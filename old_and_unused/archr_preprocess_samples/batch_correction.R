#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Martijn van der Werf
# Name: batch_correction.R
# Function: perform batch correction on the ATAC data using Harmony, correcting across lanes
#
############################################################################################################################


setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)

addArchRThreads(threads=4)
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject('Project_all_cells')

proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

saveArchRProject(ArchRProj=proj, outputDirectory='Project_all_cells/')