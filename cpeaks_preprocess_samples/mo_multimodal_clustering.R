#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_multimodal_clustering.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)

####################
# Functions        #
####################



####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
mem_requested <- 1900
options(future.globals.maxSize = (mem_requested * .95) * 1000 * 1024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

mo_multimodal_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
mo_multi_monocyte_loc <- paste(mo_multimodal_loc, 'mo_multimodal_monocyte_1_80_20240521.rds', sep = '')

# read object
mo_multi_monocyte <- readRDS(mo_multi_monocyte_loc)

# ATAC clustering
DefaultAssay(mo_multi_monocyte) <- "peaks"
mo_multi_monocyte <- RunTFIDF(mo_multi_monocyte)
mo_multi_monocyte <- FindTopFeatures(mo_multi_monocyte, min.cutoff = 'q0')
mo_multi_monocyte <- RunSVD(mo_multi_monocyte)
mo_multi_monocyte <- FindNeighbors(mo_multi_monocyte, dims = 2:30, graph.name = 'peaks_snn', reduction = 'lsi')
mo_multi_monocyte <- FindClusters(mo_multi_monocyte, resolution = 1.2, verbose = T, graph.name = 'peaks_snn')
mo_multi_monocyte <- RunUMAP(mo_multi_monocyte, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_", return.model = T)

# RNA clustering
DefaultAssay(mo_multi_monocyte) <- "RNA"
mo_multi_monocyte <- NormalizeData(mo_multi_monocyte)
mo_multi_monocyte <- FindVariableFeatures(mo_multi_monocyte)
mo_multi_monocyte <- ScaleData(mo_multi_monocyte)
mo_multi_monocyte <- RunPCA(mo_multi_monocyte, verbose = T)
mo_multi_monocyte <- FindNeighbors(mo_multi_monocyte, dims = 1:30)
mo_multi_monocyte <- FindClusters(mo_multi_monocyte, resolution = 1.2, verbose = T)
mo_multi_monocyte <- RunUMAP(mo_multi_monocyte, dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", return.model = T)

# multimodal clustering
mo_multi_monocyte <- FindMultiModalNeighbors(
  mo_multi_monocyte, reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30), modality.weight.name = c("RNA.weight", "ATAC.weigth")
)
mo_multi_monocyte <- RunUMAP(mo_multi_monocyte, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = T)
mo_multi_monocyte <- FindClusters(mo_multi_monocyte, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = T)

# set NA LONG COVID to 'unknown'
mo_multi_monocyte@meta.data[is.na(mo_multi_monocyte@meta.data[['LONG_COVID']]), 'LONG_COVID'] <- 'unknown'

# make some dimplots
p_peaks_clusters <- DimPlot(mo_multi_monocyte, reduction = 'umap.atac', group.by = 'peaks_snn_res.1.2')
p_peaks_condition <- DimPlot(mo_multi_monocyte, reduction = 'umap.atac', group.by = 'condition_final')
p_peaks_covid <- DimPlot(mo_multi_monocyte, reduction = 'umap.atac', group.by = 'LONG_COVID') + scale_color_manual(values = list('case' = 'darkred', control = 'darkcyan', 'unknown' = 'gray'))
p_rna_clusters <- DimPlot(mo_multi_monocyte, reduction = 'umap.rna', group.by = 'RNA_snn_res.1.2')
p_rna_condition <- DimPlot(mo_multi_monocyte, reduction = 'umap.rna', group.by = 'condition_final')
p_rna_covid <- DimPlot(mo_multi_monocyte, reduction = 'umap.rna', group.by = 'LONG_COVID') + scale_color_manual(values = list('case' = 'darkred', control = 'darkcyan', 'unknown' = 'gray'))
p_mm_clusters <- DimPlot(mo_multi_monocyte, reduction = 'wnn.umap', group.by = 'wsnn_res.2')
p_mm_condition <- DimPlot(mo_multi_monocyte, reduction = 'wnn.umap', group.by = 'condition_final')
p_mm_covid <- DimPlot(mo_multi_monocyte, reduction = 'wnn.umap', group.by = 'LONG_COVID') + scale_color_manual(values = list('case' = 'darkred', control = 'darkcyan', 'unknown' = 'gray'))
# all together
plot_grid(p_peaks_clusters + theme(legend.position = "none"), 
        p_peaks_condition, 
        p_peaks_covid, 
        p_rna_clusters + theme(legend.position = "none"), 
        p_rna_condition, 
        p_rna_covid, 
        p_mm_clusters + theme(legend.position = "none"), 
        p_mm_condition, 
        p_mm_covid,
        nrow = 3, 
        ncol = 3
)
# save the result
saveRDS(mo_multi_monocyte, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_monocyte_1_80_20240521_clust.rds')
