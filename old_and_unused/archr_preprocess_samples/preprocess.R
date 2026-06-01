#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Martijn van der Werf, Roy Oelen
# Name: preprocess.R
# Function: add demultiplexing information from the RNA and freemuxlet runs, to remove doublets
#
############################################################################################################################

setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)

addArchRThreads(threads=16)
addArchRGenome("hg38")
h5disableFileLocking()

# 945,804 cells after TSS > 4 and nFrags >= 2000 
proj <- loadArchRProject('Project_all_cells')

# Get demultiplexing information 
rna_lanes <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_sample_matched.tsv', sep = '\t', header = T) 

# Load Freemuxlet output 
atac_files <- list.files('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/freemuxlet/atac/assignments/', pattern = '.tsv', full.names = T) 
atac_tables <- lapply(atac_files, read.table, sep = '\t', header = TRUE)
atac_lanes <- do.call(rbind ,atac_tables)#[,c(1,2,3)]

# Create names matching ArchR cellColData rownames
rna_lanes$barcode_lane <-paste(rna_lanes$barcode, rna_lanes$lane,sep="_")
atac_lanes$barcode_lane <- paste(gsub("-1", "", atac_lanes$BARCODE), atac_lanes$lanes, sep = "_")

# Common barcodes

# Fetch doublets
souporcell_doublets <- rna_lanes[rna_lanes$status == "doublet",]$barcode_lane
freemuxlet_doublets <- atac_lanes[atac_lanes$DROPLET.TYPE == "DBL",]$barcode_lane
freemuxlet_no_doublets <- atac_lanes[atac_lanes$DROPLET.TYPE == "SNG",]$barcode_lane

# 853,360 cells left after Souporcell removal
proj_doublet_filtered <- proj[!paste(gsub("(.*#)", "", gsub("-1", "", rownames(proj@cellColData))), proj$Sample, sep = "_") %in% souporcell_doublets]
# 840,065 cells left after Freemuxlet removal
proj_doublet_filtered <- proj_doublet_filtered[!paste(gsub("(.*#)", "", gsub("-1", "", rownames(proj_doublet_filtered@cellColData))), proj_doublet_filtered$Sample, sep = "_") %in% freemuxlet_doublets]

# Save project
saveArchRProject(proj_doublet_filtered, outputDirectory='Project_all_cells')