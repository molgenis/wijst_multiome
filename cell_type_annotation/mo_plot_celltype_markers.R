#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_plot_celltype_markers.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# for reading the Seurat object
library(Seurat)
# for plotting
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

label_dict <- function(){
  label_dict <- list()
  # condition combinations
  label_dict[['UT24hCA']] <- 'UT-24hCA'
  # conditions
  label_dict[['UT']] <- 'C'
  label_dict[['C']] <- 'C'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't6-8w'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["bulk"]] <- "bulk-like"
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
  # do the datasets
  label_dict[["mo"]] <- "multiome"
  label_dict[["1m"]] <- "NC 2022"
  label_dict[["1M"]] <- "NC 2022"
  return(label_dict)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UTBaseline"]] <- "khaki2"
  color_coding[["UTt24h"]] <- "khaki4"
  color_coding[["UTt8w"]] <- "paleturquoise1"
  color_coding[["Baselinet24h"]] <- "paleturquoise3"
  color_coding[["Baselinet8w"]] <- "rosybrown1"
  color_coding[["t24ht8w"]] <- "rosybrown3"
  color_coding[["UT\nBaseline"]] <- "khaki2"
  color_coding[["UT\nt24h"]] <- "khaki4"
  color_coding[["UT\nt8w"]] <- "paleturquoise1"
  color_coding[["Baseline\nt24h"]] <- "paleturquoise3"
  color_coding[["Baseline\nt8w"]] <- "rosybrown1"
  color_coding[["t24h\nt8w"]] <- "rosybrown3"
  color_coding[["UT-Baseline"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["Baseline-t24h"]] <- "paleturquoise3"
  color_coding[["Baseline-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  color_coding[["UT-t0"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["HC-t0"]] <- "khaki2"
  color_coding[["t0-HC"]] <- "khaki2"
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["t24h-HC"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t8w-HC"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t24h-t0"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t0-t8w"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t8w-t0"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t24h-t8w"]] <- "#C00040" #"rosybrown3"
  color_coding[["t8w-t24h"]] <- "#C00040" #"rosybrown3"
  # set condition colors
  color_coding[["HC"]] <- "grey"
  color_coding[["C"]] <- "grey"
  color_coding[["Controls"]] <- "grey"
  color_coding[["t0"]] <- "pink"
  color_coding[["t24h"]] <- "red"
  color_coding[["t8w"]] <- "purple"
  color_coding[["t6-8w"]] <- "purple"
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
  # up and down
  color_coding[["Bulk up"]] <- "black"
  color_coding[["CD4T up"]] <- "#153057"
  color_coding[["CD8T up"]] <- "#009DDB"
  color_coding[["monocyte up"]] <- "#EDBA1B"
  color_coding[["NK up"]] <- "#E64B50"
  color_coding[["B up"]] <- "#71BC4B"
  color_coding[["DC up"]] <- "#965EC8"
  color_coding[["CD4+ T up"]] <- "#153057"
  color_coding[["CD8+ T up"]] <- "#009DDB"
  # percentage of whitening 
  pct_whitening=40
  color_coding[["Bulk down"]] <- colorRampPalette(c(color_coding[["Bulk up"]], "white"))(100)[pct_whitening]
  color_coding[["CD4T down"]] <- colorRampPalette(c(color_coding[["CD4T up"]], "white"))(100)[pct_whitening]
  color_coding[["CD8T down"]] <- colorRampPalette(c(color_coding[["CD8T up"]], "white"))(100)[pct_whitening]
  color_coding[["monocyte down"]] <- colorRampPalette(c(color_coding[["monocyte up"]], "white"))(100)[pct_whitening]
  color_coding[["NK down"]] <- colorRampPalette(c(color_coding[["NK up"]], "white"))(100)[pct_whitening]
  color_coding[["B down"]] <- colorRampPalette(c(color_coding[["B up"]], "white"))(100)[pct_whitening]
  color_coding[["DC down"]] <- colorRampPalette(c(color_coding[["DC up"]], "white"))(100)[pct_whitening]
  color_coding[["CD4+ T down"]] <- colorRampPalette(c(color_coding[["CD4+ T up"]], "white"))(100)[pct_whitening]
  color_coding[["CD8+ T down"]] <- colorRampPalette(c(color_coding[["CD8+ T up"]], "white"))(100)[pct_whitening]
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}



rename_labels <- function(vector_to_rename) {
  # get the labels that are present
  label_renaming <- get_label_dict()
  # convert vector
  vector_to_rename <- as.character(vector_to_rename)
  # now check which labels we are missing
  missing_renames <- setdiff(unique(vector_to_rename), names(label_renaming))
  # add those renames as not being renames
  for (missing_rename in missing_renames) {
    label_renaming[[missing_rename]] <- missing_rename
  }
  # now replace each value with the rename
  renamed_vector <- as.vector(unlist(label_renaming[vector_to_rename]))
  # and return that
  return(renamed_vector)
}


wide_to_long_table <- function(wide_table, column_name='seurat_clusters') {
  # the genes are the rownames
  genes <- rownames(wide_table)
  # the colnames we will make into a new column
  colname_values <- colnames(wide_table)
  # create a new dataframe
  long_table <- data.frame(ident = rep(NA, times = length(genes) * length(colname_values)), gene = rep(NA, times = length(genes) * length(colname_values)), expression = rep(NA, times = length(genes) * length(colname_values)))
  # we will fill the new table row by row, so we need to keep an index
  long_table_row_i <- 1
  # now we'll check each column in the original
  for (colvalue in colname_values) {
    # check each gene
    for (gene in genes) {
      # get the actual value from the table
      value_combination <- wide_table[gene, colvalue]
      # add to the long table
      long_table[long_table_row_i, 'ident'] <- colvalue
      long_table[long_table_row_i, 'gene'] <- gene
      long_table[long_table_row_i, 'expression'] <- value_combination
      # increase the index
      long_table_row_i <- long_table_row_i + 1
    }
  }
  # reset column names
  colnames(long_table) <- c(column_name, 'gene', 'expression')
  return(long_table)
}


plot_expression_and_pct <- function(expression_pct_table, ident_column='seurat_clusters', gene_column='gene', expression_column='expression', percentage_column='percentage_expressed', colour_ticks=T, colour_mapping=NULL,ident_order=NULL, gene_order=NULL, paper_style=T, align_right=T, legendless=F) {
  # order the idents if the order was supplied
  if (!is.null(ident_order)) {
    expression_pct_table[[ident_column]] <- factor(expression_pct_table[[ident_column]], levels = rev(ident_order))
  }
  # and the genes
  if(!is.null(gene_order)) {
    expression_pct_table[[gene_column]] <- factor(expression_pct_table[[gene_column]], levels = gene_order)
  }
  # make the plot
  p <- ggplot(data=NULL, mapping = aes(x = expression_pct_table[[gene_column]], y = expression_pct_table[[ident_column]], size = expression_pct_table[[percentage_column]], color = expression_pct_table[[expression_column]])) + geom_point() + scale_colour_gradient(name = 'Average expression', low = 'lightblue', high = 'darkblue') + labs(size = percentage_column)
  # align to the right side
  if (align_right) {
    p <- p + scale_y_discrete(position = 'right')
  }
  # and labels
  p <- p + xlab(gene_column) + ylab(ident_column)
  # paper style if requested
  if (paper_style) {
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  # remove legend if requested
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  # rotate labels
  p <- p + theme(axis.text.x = element_text(angle = 90))
  # colour the ticks if requested
  if (colour_ticks) {
    # get the labels
    y_labels <- levels(expression_pct_table[[ident_column]])
    # get the appropriate colours
    y_colours <- as.vector(unlist(colour_mapping[as.character(y_labels)]))
    p <- p + theme(axis.ticks = element_line(),      # Change ticks line fo all axes
                   axis.ticks.y = element_line(linewidth = 5, colour = y_colours),    # Change x axis ticks only
    )
  }
  return(p)
}


####################
# Main code        #
####################

# location of the Seurat object
mo_rna_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds'
# read the object
mo_rna <- readRDS(mo_rna_object_loc)

# the list of features to plot
ng2018_marker_genes <- c(
  'CD3D',
  'CD8A',
  'GZMB',
  'KLRC1',
  'FCGR3A',
  'MS4A1',
  'CD79A',
  'CD14',
  'ITGAX',
  'PPBP'
)

cell_type_order <- c(
  'CD4T',
  'CD4+ T', 
  'CD8T', 
  'CD8+ T', 
  'NK',
  'B',
  'plasmablast', 
  'Plasmablast', 
  'monocyte', 
  'Monocyte', 
  'DC',
  'megakaryocyte', 
  'Megakaryocyte', 
  'T_other', 
  'T other', 
  'other T', 
  'Other T',
  'unannotated'
)

# calculate the average expression
mo_rna_exp <- AverageExpression(mo_rna, group.by = 'celltype_imputed_lowerres', features = ng2018_marker_genes)
# extract which ones we care about
mo_rna_exp_sct_genes <- data.frame(mo_rna_exp[['SCT']][ng2018_marker_genes, ])
# fix changes to spaces
colnames(mo_rna_exp_sct_genes) <- gsub('\\.', '_', colnames(mo_rna_exp_sct_genes))
# turn into a long table
mo_rna_ctlr_avg_exp_sct_genes_long <- wide_to_long_table(mo_rna_exp_sct_genes, column_name = 'celltype_imputed_lowerres')
# get pct info
mo_rna_dp <- DotPlot(mo_rna, features = ng2018_marker_genes, assay = 'SCT', group.by = 'celltype_imputed_lowerres') + theme(axis.text.x = element_text(angle = 90)) + xlab('gene') + ylab('Seurat cluster')
mo_rna_dp
# extract the pct.exp from the dotplot
feature_and_ident_ctlr_dotplot <- paste(as.character(mo_rna_dp$data$id), as.character(mo_rna_dp$data$features.plot), sep = '_')
feature_and_ident_ctlr_avgexp <- paste(as.character(mo_rna_ctlr_avg_exp_sct_genes_long[['celltype_imputed_lowerres']]), as.character(mo_rna_ctlr_avg_exp_sct_genes_long[['gene']]), sep = '_')
# now add the pct to the table
mo_rna_ctlr_avg_exp_sct_genes_long[['percentage_expressed']] <- mo_rna_dp$data[match(feature_and_ident_ctlr_avgexp, feature_and_ident_ctlr_dotplot), 'pct.exp']
# rename cell types
mo_rna_ctlr_avg_exp_sct_genes_long[['cell_type']] <- rename_labels(mo_rna_ctlr_avg_exp_sct_genes_long[['celltype_imputed_lowerres']])
# plot the percentages
p_pct_exp <- plot_expression_and_pct(mo_rna_ctlr_avg_exp_sct_genes_long[mo_rna_ctlr_avg_exp_sct_genes_long[['celltype_imputed_lowerres']] != 'unannotated', ], ident_column = 'cell_type', ident_order = cell_type_order, gene_order = ng2018_marker_genes, colour_ticks = T, colour_mapping = get_color_coding_dict(), align_right = F) + ylab('cell type') + labs(size = 'Percentage cells expressing', colour = 'Average expression') + xlab('Gene')
p_pct_exp
# save this
ggsave('~/plots/mo_celltype_markers.pdf', plot = p_pct_exp, width = 9, height = 6)
