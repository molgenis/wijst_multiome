#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_plot_seurat_object.R
# Function: various plot of the MO object
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(ggplot2)
library(RColorBrewer)


####################
# Functions        #
####################


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
  color_coding[["unannotated"]] <- "gray"
  return(color_coding)
}



#' get a vector of as distinct possible colours
#' 
#' @param number_of_colours how many colours to return
#' @param use_sampling whether or not to randomly extract the colours instead of grabbing the first n colours
#' @param color_indices (optional, not used by default) if specific colours are needed, supply the indices of the colours here. Use 'get_available_colours_grid' to get the colours and their indices
#' @returns a vector of colours
#' 
sample_many_colours <- function(number_of_colours, use_sampling=F, color_indices=NULL) {
  # get all colours from the 'quality' palettes
  quality_colour_palettes <- brewer.pal.info[brewer.pal.info[['category']] == 'qual', ]
  # save each palette
  colours_per_palette <- list()
  # apply over each palette
  for (i in 1:nrow(quality_colour_palettes)) {
    # get the name of the palette
    palette_name <- rownames(quality_colour_palettes)[i]
    # get the number of colours in the palette
    palette_max_colours <- quality_colour_palettes[i, 'maxcolors']
    # use brewer.pal to get all colours
    colours_palette <- brewer.pal(palette_max_colours, palette_name)
    # put result in the list
    colours_per_palette[[palette_name]] <- colours_palette
  }
  # merge all palettes
  all_colours <- do.call('c', colours_per_palette)
  # randomly get colours from the palette
  max_possible_colours <- length(all_colours)
  if (is.null(number_of_colours)) {
    message('no number of colors supplied, assuming color indices have been')
  }
  else if (number_of_colours > max_possible_colours) {
    message(paste('requesting more colours than is possible: ', as.character(number_of_colours), ' vs ', max_possible_colours, ', returning max possible', sep = ''))
    number_of_colours <- max_possible_colours
  }
  colours_to_return <- NULL
  # specific colours we like (the indices)
  if (!is.null(color_indices)) {
    colours_to_return <- all_colours[color_indices]
  }
  # or use sampling
  else if (use_sampling) {
    colours_to_return <- sample(all_colours, number_of_colours)
  }
  # or the first x colours
  else {
    colours_to_return <- all_colours[1 : number_of_colours]
  }
  return(colours_to_return)
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
  # unannotated ones
  label_dict[['unannotated']] <- 'unannotated'
  return(label_dict)
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


# location of where to place the objects
seurat_objects_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
mo_normalized_loc <- paste(seurat_objects_loc, 'mo_all_20240223_seuratv5_normalized.rds', sep = '')

# read object
mo <- readRDS(mo_normalized_loc)

# get the number of cells QC
rna_filter_table <- data.frame(filter_step = factor(c('cellbender', 'doublet', 'correlation', 'min_200_umi', 'min_200_unique_genes', 'umi_mad', 'unique_gene_mad'), levels = c('cellbender', 'doublet', 'correlation', 'min_200_umi', 'min_200_unique_genes', 'umi_mad', 'unique_gene_mad')), 
                               cells_left = c(1423526, 1037166, 992701, 914447, 893663, 874201, 874200))
# get some colours
colours_rna_filter <- sample_many_colours(n = length(unique(rna_filter_table[['filter_step']])))
names(colours_rna_filter) <- as.character(unique(rna_filter_table[['filter_step']]))
# make the plot
ggplot(data = rna_filter_table, mapping = aes(x = filter_step, y = cells_left, fill = filter_step)) + 
  # as bar
  geom_bar(stat = 'identity') + 
  # add text as well
  geom_text(mapping = aes(label = cells_left), vjust = 1.6, color = 'black', size = 5) +
  # with own colours
  scale_fill_manual(values = colours_rna_filter) + 
  # with clean background
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # no legend
  theme(legend.position = 'none') +
  # better x label
  xlab('filtering step') +
  # better y label
  ylab('cells left after filter') +
  # better title
  ggtitle('cells left after each filtering step') + 
  # change x axis text size
  theme(axis.text = element_text(size = 12), axis.title = element_text(size=14), axis.title.x = element_text(size = 16))

# get celltype numbers
celltypes_lane_numbers <- get_celltype_numbers(mo@meta.data, sample_column = 'lane', cell_type_column = 'predicted.mo_10x_cell_type.lowerres')
# per stim as well
celltypes_stim_numbers <- get_celltype_numbers(mo@meta.data[!is.na(mo@meta.data[['final_condition']]), ], sample_column = 'final_condition', cell_type_column = 'predicted.mo_10x_cell_type.lowerres')
celltypes_stim_numbers[['final_condition']] <- factor(celltypes_stim_numbers[['final_condition']], levels = c('UT', '24hCa'))
# and per sample+lane
mo@meta.data[['lane_sample']] <- paste(mo@meta.data[['lane']], mo@meta.data[['sample_final']])
celltypes_sample_numbers <- get_celltype_numbers(mo@meta.data[!is.na(mo@meta.data[['lane_sample']]), ], sample_column = 'lane_sample', cell_type_column = 'predicted.mo_10x_cell_type.lowerres')
# make plots of it
plot_celltype_abundance(celltypes_lane_numbers[!is.na(celltypes_lane_numbers[['predicted.mo_10x_cell_type.lowerres']]),], sample_column = 'lane', cell_type_column = 'predicted.mo_10x_cell_type.lowerres', number_column = 'proportion', angle_x_labels = T)
plot_celltype_abundance(celltypes_stim_numbers[!is.na(celltypes_stim_numbers[['predicted.mo_10x_cell_type.lowerres']]) & !is.na(celltypes_stim_numbers[['final_condition']]),], sample_column = 'final_condition', cell_type_column = 'predicted.mo_10x_cell_type.lowerres', number_column = 'proportion', angle_x_labels = T, legendless = T)
plot_celltype_abundance(celltypes_sample_numbers[!is.na(celltypes_sample_numbers[['predicted.mo_10x_cell_type.lowerres']]),], sample_column = 'lane_sample', cell_type_column = 'predicted.mo_10x_cell_type.lowerres', number_column = 'proportion', pointless = T, angle_x_labels = F, legendless = T)

# plot cell type assignment
DimPlot(mo, group.by = 'predicted.mo_10x_cell_type') + scale_colour_manual(values = unlist(as.vector(sample_many_colours(length(unique(mo@meta.data[['predicted.mo_10x_cell_type']]))))))
# or low level
DimPlot(mo, group.by = 'predicted.mo_10x_cell_type.lowerres') + scale_colour_manual(values = get_color_coding_dict())
FeaturePlot(mo, features = c('predicted.mo_10x_cell_type.score'))

# plot by condition
DimPlot(mo[, !is.na(mo@meta.data[['final_condition']])], group.by = 'final_condition') + scale_colour_manual(values = get_color_coding_dict())

# get ncells per donor and lane
ncell_donor_lane <- data.frame(table(mo@meta.data[, c('lane', 'unconfined_best_match_sample', 'unconfined_best_match_correlation')]))
ncell_donor_lane <- ncell_donor_lane[ncell_donor_lane[['Freq']] > 0, ]
ncell_donor_lane[['unconfined_best_match_correlation']] <- as.numeric(as.character(ncell_donor_lane[['unconfined_best_match_correlation']]))
# plot the correlation vs the number of cells
ggplot(data = ncell_donor_lane, mapping = aes(x = Freq, y = unconfined_best_match_correlation)) + 
  # do points
  geom_point() +
  # from 0 to 1
  ylim(0, 1) +
  # and from zero to max
  xlim(0, max(ncell_donor_lane[['Freq']])) + 
  # with clean background
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # no legend
  theme(legend.position = 'none') +
  # better x label
  xlab('number of cells') +
  # better y label
  ylab('best correlation') +
  # better title
  ggtitle('number of cells vs best correlation') + 
  # change x axis text size
  theme(axis.text = element_text(size = 12), axis.title = element_text(size=14), axis.title.x = element_text(size = 16))
