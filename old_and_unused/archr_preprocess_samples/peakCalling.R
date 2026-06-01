#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Martijn van der Werf
# Name: peakCalling.R
# Function: perform peak calling using macs2 per (imputed) celltype
#
############################################################################################################################


setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads=16)
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject('Project_common_atac_rna_cells/')

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "cell_type_lowerres_imputed")

macs2_path <- "/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/miniconda3/envs/macs2/bin/macs2"

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "cell_type_lowerres_imputed", 
    pathToMacs2 = macs2_path,
    maxPeaks=2000000
)

# Save peaksets as bed files
library(Repitools)
peaks_monocytes <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/monocyte-reproduciblePeaks.gr.rds'))
peaks_b <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/B-reproduciblePeaks.gr.rds'))
peaks_cd4t <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/CD4T-reproduciblePeaks.gr.rds'))
peaks_cd8t <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/CD8T-reproduciblePeaks.gr.rds'))
peaks_dc <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/DC-reproduciblePeaks.gr.rds'))
peaks_nk <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/NK-reproduciblePeaks.gr.rds'))

write.table(peaks_monocytes, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/monocytes.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_b, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/B.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_cd4t, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD4T.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_cd8t, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD8T.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_dc, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/DC.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_nk, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/NK.bed', sep='\t', row.names=F, col.names=T, quote=F)

# Compare unique and overlapping peaks in our bed files
unique_monocytes <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/monocyte_unique_peaks.bed', sep='\t', header=F)
unique_b <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/B_unique_peaks.bed', sep='\t', header=F)
unique_cd4t <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD4T_unique_peaks.bed', sep='\t', header=F)
unique_cd8t <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD8T_unique_peaks.bed', sep='\t', header=F)
unique_dc <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/DC_unique_peaks.bed', sep='\t', header=F)
unique_nk <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/NK_unique_peaks.bed', sep='\t', header=F)

peak_comparison <- data.frame(cell_type=c('Monocyte', 'B', 'CD4T', 'CD8T', 'DC','NK'),
                   total_peaks=c(nrow(peaks_monocytes), nrow(peaks_b), nrow(peaks_cd4t), nrow(peaks_cd8t), nrow(peaks_dc), nrow(peaks_nk)),
                   overlapping_peaks=c(nrow(peaks_monocytes) - nrow(unique_monocytes), nrow(peaks_b) - nrow(unique_b), nrow(peaks_cd4t) - nrow(unique_cd4t), 
                                       nrow(peaks_cd8t) - nrow(unique_cd8t), nrow(peaks_dc) - nrow(unique_dc), nrow(peaks_nk) - nrow(unique_dc)),
                   unique_peaks=c(nrow(unique_monocytes), nrow(unique_b), nrow(unique_cd4t), nrow(unique_cd8t), nrow(unique_dc), nrow(unique_nk)))

peak_comparison$percentage_unique <- round(peak_comparison$unique_peaks / peak_comparison$total_peaks * 100, 2)

write.table(peak_comparison, 'Project_common_atac_rna_cells/PeakCalls/peak_comparison_celltypes.tsv', col.names=T,  row.names=F, quote=F, sep='\t')setwd("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/archr_preprocess_samples/")
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads=16)
addArchRGenome("hg38")
h5disableFileLocking()

proj <- loadArchRProject('Project_common_atac_rna_cells/')

proj <- addGroupCoverages(ArchRProj = proj, groupBy = "cell_type_lowerres_imputed")

macs2_path <- "/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/miniconda3/envs/macs2/bin/macs2"

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "cell_type_lowerres_imputed", 
    pathToMacs2 = macs2_path,
    maxPeaks=2000000
)

# Save peaksets as bed files
library(Repitools)
peaks_monocytes <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/monocyte-reproduciblePeaks.gr.rds'))
peaks_b <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/B-reproduciblePeaks.gr.rds'))
peaks_cd4t <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/CD4T-reproduciblePeaks.gr.rds'))
peaks_cd8t <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/CD8T-reproduciblePeaks.gr.rds'))
peaks_dc <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/DC-reproduciblePeaks.gr.rds'))
peaks_nk <- Repitools::annoGR2DF(readRDS('Project_common_atac_rna_cells/PeakCalls/NK-reproduciblePeaks.gr.rds'))

write.table(peaks_monocytes, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/monocytes.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_b, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/B.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_cd4t, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD4T.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_cd8t, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD8T.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_dc, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/DC.bed', sep='\t', row.names=F, col.names=T, quote=F)
write.table(peaks_nk, 'Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/NK.bed', sep='\t', row.names=F, col.names=T, quote=F)

# Compare unique and overlapping peaks in our bed files
unique_monocytes <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/monocyte_unique_peaks.bed', sep='\t', header=F)
unique_b <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/B_unique_peaks.bed', sep='\t', header=F)
unique_cd4t <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD4T_unique_peaks.bed', sep='\t', header=F)
unique_cd8t <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/CD8T_unique_peaks.bed', sep='\t', header=F)
unique_dc <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/DC_unique_peaks.bed', sep='\t', header=F)
unique_nk <- read.table('Project_common_atac_rna_cells/PeakCalls/Beds_celltype_peaks/NK_unique_peaks.bed', sep='\t', header=F)

peak_comparison <- data.frame(cell_type=c('Monocyte', 'B', 'CD4T', 'CD8T', 'DC','NK'),
                   total_peaks=c(nrow(peaks_monocytes), nrow(peaks_b), nrow(peaks_cd4t), nrow(peaks_cd8t), nrow(peaks_dc), nrow(peaks_nk)),
                   overlapping_peaks=c(nrow(peaks_monocytes) - nrow(unique_monocytes), nrow(peaks_b) - nrow(unique_b), nrow(peaks_cd4t) - nrow(unique_cd4t), 
                                       nrow(peaks_cd8t) - nrow(unique_cd8t), nrow(peaks_dc) - nrow(unique_dc), nrow(peaks_nk) - nrow(unique_dc)),
                   unique_peaks=c(nrow(unique_monocytes), nrow(unique_b), nrow(unique_cd4t), nrow(unique_cd8t), nrow(unique_dc), nrow(unique_nk)))

peak_comparison$percentage_unique <- round(peak_comparison$unique_peaks / peak_comparison$total_peaks * 100, 2)

write.table(peak_comparison, 'Project_common_atac_rna_cells/PeakCalls/peak_comparison_celltypes.tsv', col.names=T,  row.names=F, quote=F, sep='\t')