#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cluster_seurat_objects.R
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

ref10xmo_predictions_to_lower_res_mapping <- function() {
  high_to_low <- list()
  high_to_low[['CD4 Naive']] <- 'CD4T'
  high_to_low[['CD4 TCM']] <- 'CD4T'
  high_to_low[['CD8 Naive']] <- 'CD8T'
  high_to_low[['CD16 Mono']] <- 'monocyte'
  high_to_low[['NK']] <- 'NK'
  high_to_low[['Treg']] <- 'T_other'
  high_to_low[['CD14 Mono']] <- 'monocyte'
  high_to_low[['cDC']] <- 'DC'
  high_to_low[['CD8 TEM_1']] <- 'CD8T'
  high_to_low[['Intermediate B']] <- 'B'
  high_to_low[['Naive B']] <- 'B'
  high_to_low[['Plasma']] <- 'plasmablast'
  high_to_low[['CD4 TEM']] <- 'CD4T'
  high_to_low[['MAIT']] <- 'T_other'
  high_to_low[['Memory B']] <- 'B'
  high_to_low[['gdT']] <- 'T_other'
  high_to_low[['pDC']] <- 'DC'
  high_to_low[['CD8 TEM_2']] <- 'CD8T'
  high_to_low[['HSPC']] <- 'hemapoietic_stem'
  return(high_to_low)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- "lightgrey"
  color_coding[["24hCA"]] <- "forestgreen"
  color_coding[["24hCa"]] <- "forestgreen"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  color_coding[["T_other"]] <- "#FF63B6"
  color_coding[["hemapoietic_stem"]] <- "#8B8000"
  color_coding[["hemapoietic stem"]] <- "#8B8000"
  return(color_coding)
}


label_dict <- function(){
  label_dict <- list()
  # condition combinations
  label_dict[["UT"]] <- "UT"
  label_dict[["24hCa"]] <- "24hCA"
  label_dict[["24hCA"]] <- "24hCA"
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["T_other"]] <- "other T"
  label_dict[["hemapoietic_stem"]] <- "hemapoietic stem"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  return(label_dict)
}


#' get a dataframe with the proportons and cell numbers per sample
#' 
#' @param metadata the metadata to calculate the cell type numbers and proportions from
#' @param sample_column the sample denoting the sample
#' @param cell_type_column the column denoting the cell type
#' @returns a dataframe with the proportion and number of each cell type per sample
#' 
get_celltype_numbers <- function(metadata, sample_column='sample', cell_type_column='cell_type') {
  # get the numbers in total
  numbers_totals <- data.frame(table(metadata[, c(sample_column, cell_type_column)]))
  # init proportion
  numbers_totals[['prop']] <- 0
  # now add proportions
  for (sample_name in unique(metadata[[sample_column]])) {
    # subset to that sample
    numbers_sample <- numbers_totals[!is.na(numbers_totals[[sample_column]]) & numbers_totals[[sample_column]] == sample_name, ]
    # get the total for that sample
    numbers_sample_total <- sum(numbers_sample[['Freq']])
    # check each cell type
    for (cell_type in unique(numbers_sample[[cell_type_column]])) {
      # get the number
      numbers_celltype_sample <- numbers_sample[!is.na(numbers_sample[[cell_type_column]]) & numbers_sample[[cell_type_column]] == cell_type, 'Freq']
      # get the proportion
      prop_celltype_sample <- numbers_celltype_sample / numbers_sample_total
      # add to the table
      numbers_totals[!is.na(numbers_totals[[sample_column]]) & 
                       numbers_totals[[sample_column]] == sample_name &
                       !is.na(numbers_totals[[cell_type_column]]) & 
                       numbers_totals[[cell_type_column]] == cell_type, 'prop'] <- prop_celltype_sample
    }
  }
  # set the column names
  colnames(numbers_totals) <- c(sample_column, cell_type_column, 'number', 'proportion')
  return(numbers_totals)
}

#' get a dataframe with the proportons and cell numbers per sample
#' 
#' @param celltype_numbers the metadata to calculate the cell type numbers and proportions from
#' @param sample_column the sample denoting the sample
#' @param cell_type_column the column denoting the cell type
#' @param number_column the column denoting number to plot
#' @param use_label_dict update labels with nicer ones, using the label dict
#' @param pointless plot the ticks and names on the x axis
#' @param legendless plot the legend
#' @param paper_style use paper style plot, with more whitespace
#' @param angle_x_labels angle the x label 90 degrees
#' @param use_colour_coding_dict use colors from the color dict
#' @returns a plot of cell numbers/proportions per sample
#' 
plot_celltype_abundance <- function(celltype_numbers, sample_column='sample', cell_type_column='cell_type', number_column='proportion', use_label_dict=T, pointless=F, legendless=F, paper_style=T, angle_x_labels=F, use_colour_coding_dict=T) {
  # use some better labels if requested
  if (use_label_dict) {
    celltype_numbers[[cell_type_column]] <- as.vector(unlist(label_dict()[as.character(celltype_numbers[[cell_type_column]])]))
  }
  
  # init plot
  p <- ggplot(data = NULL, mapping = aes(x = celltype_numbers[[sample_column]], y = celltype_numbers[[number_column]], fill = celltype_numbers[[cell_type_column]])) +
    geom_bar(stat = 'identity', position = 'stack') +
    xlab(sample_column) +
    ylab(number_column)
  
  # add colours if requested
  if (use_colour_coding_dict) {
    p <- p + scale_fill_manual(name = cell_type_column, values = get_color_coding_dict())
  }
  
  # all our options
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if(paper_style){
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if(angle_x_labels){
    p <- p + theme(axis.text.x = element_text(angle = 90))
  }
  return(p)
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

# install the dataset
#InstallData("pbmcMultiome")

# load both modalities
pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")

# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc.atac) <- annotations

# repeat QC steps performed in the WNN vignette
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
#pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")

#pbmc.rna@meta.data[['cell_type_lowerres']] <- as.vector(unlist(ref10xmo_predictions_to_lower_res_mapping()[pbmc.rna@meta.data[['seurat_annotations']]]))
#celllevel_ct_numbers_10x <- get_celltype_numbers(pbmc.rna@meta.data, sample_column = 'orig.ident', cell_type_column = 'cell_type_lowerres')
#plot_celltype_abundance(celllevel_ct_numbers_10x, sample_column = 'orig.ident', cell_type_column = 'cell_type_lowerres', angle_x_labels = T) + guides(fill = guide_legend(title = 'cell type'))

# We exclude the first dimension as this is typically correlated with sequencing depth
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")

# instead just get the chromatin assay data
counts <- GetAssayData(pbmc.atac, slot = "counts")
# create a new assay, while filtering out what we already removed in the RNA data
chromatinassay <- CreateChromatinAssay(counts = counts[, !is.na(pbmc.atac@meta.data[['seurat_annotations']]) & pbmc.atac@meta.data[['seurat_annotations']] != 'filtered'], 
                                       genome = "hg38")

# put all in one object
cbmc <- CreateSeuratObject(counts = CreateAssayObject(pbmc.rna@assays$RNA$counts))
# add this assay to the previously created Seurat object
cbmc[["peaks"]] <- chromatinassay

# add the metadata back
cbmc <- AddMetaData(cbmc, metadata = pbmc.rna@meta.data['seurat_annotations'])

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
cbmc <- FindNeighbors(cbmc, dims = 1:30)
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

# add lower resolution
cbmc@meta.data[['cell_type_lowerres']] <- as.vector(unlist(ref10xmo_predictions_to_lower_res_mapping()[cbmc@meta.data[['seurat_annotations']]]))

# save result
saveRDS(cbmc, '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/cell_type_assignment/azimuth/10x_mo_reference.rds')

plot_grid(
  DimPlot(cbmc, reduction = 'umap.atac', group.by = 'peaks_snn_res.1.2'),
  DimPlot(cbmc, reduction = 'umap.atac', group.by = 'seurat_annotations'),
  DimPlot(cbmc, reduction = 'umap.atac', group.by = 'cell_type_lowerres'),
  DimPlot(cbmc, reduction = 'umap.rna', group.by = 'RNA_snn_res.1.2'),
  DimPlot(cbmc, reduction = 'umap.rna', group.by = 'seurat_annotations'),
  DimPlot(cbmc, reduction = 'umap.rna', group.by = 'cell_type_lowerres'),
  DimPlot(cbmc, reduction = 'wnn.umap', group.by = 'wsnn_res.2'),
  DimPlot(cbmc, reduction = 'wnn.umap', group.by = 'seurat_annotations'),
  DimPlot(cbmc, reduction = 'wnn.umap', group.by = 'cell_type_lowerres'),
  nrow = 3
)

ggsave('~/mo_10x_reference_clusters_celltype.pdf', width = 12, height = 12)