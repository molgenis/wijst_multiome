#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Martijn van der Werf, Roy Oelen
# Name: build-project.R
# Function: create the ArchR object
#
############################################################################################################################

setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)

addArchRThreads(threads=16)
addArchRGenome("hg38")
h5disableFileLocking()

ArrowFiles <- list.files('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/Arrowfiles_raw', full.names=T)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Project_all_cells",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)