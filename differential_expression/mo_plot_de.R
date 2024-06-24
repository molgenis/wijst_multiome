#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_plot_de.R
# Function: plot DE results
############################################################################################################################


####################
# libraries        #
####################

# plotting
library(ggplot2)
library(ggvenn)
library(cowplot)

####################
# Functions        #
####################

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


get_de_genes <- function(mast_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', lfc_cutoff=NULL, to_ens=F, symbols.to.ensg.mapping='genes.tsv', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w'), de_method='mast_seurat', condition_name='timepoint'){
  # set up per cell type
  de_per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # set up per stim combination
    de_per_condition <- list()
    # check each stim
    for(stim in stims){
      for(stim2 in stims){
        try({
          if(stim != stim2){
            print(paste(cell_type, stim, stim2, sep = ' '))
            # paste the filepath together
            filepath <- NULL
            if ("limma" == de_method) {
              #filepath <- paste(mast_output_loc, cell_type, '_', condition_name,  stim, '_', stim2, '.tsv.gz', sep = '')
              filepath <- paste(mast_output_loc, cell_type, '_', condition_name, '.tsv.gz', sep = '')
            }
            else {
              filepath <- paste(mast_output_loc, cell_type, stim,stim2, '.tsv.gz', sep = '')
            }
            if (file.exists(filepath)) {
              # read the mast output
              mast <- read.table(filepath, header=T)
              # filter to only include the significant results
              mast <- mast[mast[[pval_column]] <= 0.05, ]
              # filter for only the positive lfc if required
              if(only_positive){
                mast <- mast[mast[[lfc_column]] < 0, ]
              }
              # filter for only the positive lfc if required
              if(only_negative){
                mast <- mast[mast[[lfc_column]] > 0, ]
              }
              # filter only ones with strong enough effect if required
              if (!is.null(lfc_cutoff)) {
                mast <- mast[abs(mast[[lfc_column]]) >= lfc_cutoff, ]
              }
              # confine in some way if reporting a max number of genes
              if(!is.null(max)){
                # by p if required
                if(max_by_pval){
                  mast <- mast[order(mast[[p_val_column]]), ]
                }
                # by lfc otherwise
                else{
                  mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
                }
                # subset to the number we requested if max was set
                mast <- mast[1:max,]
              }
              # grab the genes from the column names
              genes <- rownames(mast)
              if ("limma" == de_method) {
                genes <- mast[['feature']]
              }
              # convert the symbols to ensemble IDs
              if (to_ens) {
                mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
                mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
                genes <- mapping[match(genes, mapping$V2),"V1"]
              }
              # otherwise change the Seurat replacement back
              else{
                #genes <- gsub("-", "_", genes)
              }
              de_per_condition[[paste(stim, stim2, sep = '')]] <- genes
            }
          }
        })
      }
    }
    de_per_ct[[cell_type]] <- de_per_condition
  }
  return(de_per_ct)
}



de_genes_number_to_table <- function(mast_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', lfc_cutoff=NULL, to_ens=F, symbols.to.ensg.mapping='genes.tsv', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w'), remove_na_cols=T, de_method='mast_seurat', condition_name='timepoint'){
  # get the DE genes
  de_genes_per_ct <- get_de_genes(mast_output_loc, pval_column=pval_column, sig_pval=sig_pval, max=max, max_by_pval=max_by_pval, only_positive=only_positive, only_negative=only_negative, lfc_column=lfc_column, to_ens=to_ens, symbols.to.ensg.mapping=symbols.to.ensg.mapping, cell_types=cell_types, stims=stims, de_method=de_method, lfc_cutoff=lfc_cutoff, condition_name=condition_name)
  # make all possible combinations of stims
  combs <- paste(rep(stims, each = length(stims)), stims, sep = '')
  # create matrix to store results
  number_table <- matrix(, ncol=length(combs), nrow=length(cell_types), dimnames = list(cell_types, combs))
  # check each cell type
  for(cell_type in intersect(cell_types, names(de_genes_per_ct))){
    # get for specific cell type
    de_genes_per_conditin <- de_genes_per_ct[[cell_type]]
    # check each condition combination
    for(comb in intersect(combs, names(de_genes_per_conditin))){
      # get the number of genes
      nr_of_de_genes <- length(de_genes_per_conditin[[comb]])
      # add to matrix
      number_table[cell_type, comb] <- nr_of_de_genes
    }
  }
  # remove na column
  if(remove_na_cols){
    number_table <- number_table[, colSums(is.na(number_table)) < nrow(number_table)]
  }
  # i like dataframes
  number_table <- data.frame(number_table)
  return(number_table)
}


numbers_table_to_plot <- function(numbers_table, cols_include=NULL, use_label_dict=T, use_groups_dict=T, title=NULL, pointless=F, legendless=F, grid_single=T, paper_style=F, angle_x_labels=F, colour_top_label=F, colour_celltypes=T){
  numbers_table_to_do <- numbers_table
  if(!is.null(cols_include)){
    numbers_table_to_do <- numbers_table_to_do[, cols_include, drop = F]
  }
  plot_data <- NULL
  for(cell_type in rownames(numbers_table_to_do)){
    for(condition_comb in colnames(numbers_table_to_do)){
      val <- numbers_table_to_do[cell_type, condition_comb]
      if(is.na(val)){
        val <- 0
      }
      # set labels to use in plot
      conditions_label <- condition_comb
      cell_type_label <- cell_type
      # create nicer labels if requested
      if(use_label_dict){
        conditions_label <- label_dict()[[conditions_label]]
        cell_type_label <- label_dict()[[cell_type]]
      }
      # create the row
      row_plot_data <- data.frame(de_genes=val, cell_type=cell_type_label, conditions=conditions_label)
      if(use_groups_dict){
        row_plot_data$comparison <- groups_dict()[[condition_comb]]
      }
      if(is.null(plot_data)){
        plot_data <- row_plot_data
      }
      else{
        plot_data <- rbind(plot_data, row_plot_data)
      }
    }
  }
  # create the ylim
  ylims <- c(0, max(plot_data$de_genes*1.1))
  # make the plots
  p <- NULL
  if(use_groups_dict){
    p <- ggplot(data=plot_data, aes(x=comparison, y=de_genes)) + geom_point(aes(color=conditions), size=6) + facet_grid(. ~ cell_type) + scale_color_manual(name = 'condition\ncombination', values = unlist(get_color_coding_dict()[unique(plot_data$conditions)])) + ylim(ylims)
  }
  else if(length(unique(plot_data$conditions))==1 | colour_celltypes){
    p <- NULL
    if(grid_single){
      p <- ggplot(data=plot_data, aes(x=conditions, y=de_genes)) + geom_point(aes(color=cell_type), size=6) + facet_grid(. ~ cell_type)  + scale_color_manual(values = unlist(get_color_coding_dict()[unique(plot_data$cell_type)]))
    }
    else{
      p <- ggplot(data=plot_data, aes(x=cell_type, y=de_genes)) + geom_point(aes(color=cell_type), size=6) + scale_color_manual(name = 'cell type', values = unlist(get_color_coding_dict()[unique(plot_data$cell_type)]))
    }
    #
    p <- p +
      xlab('cell type') + 
      ylab('number of DE genes') +
      theme(legend.title = element_text(size=14), 
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            axis.text.y = element_text(size=12),
            strip.text.x = element_text(size=12)) +
      ylim(ylims)
  }
  else{
    #ggplot(data=plot_data, aes(x=conditions, y=de_genes)) + geom_point(aes(color=conditions), size=6) + facet_grid(. ~ cell_type)  + scale_color_manual(name = 'condition\ncombination', values = unlist(get_color_coding_dict()[unique(plot_data$conditions)])) + theme(legend.position = 'none') +
    p <- ggplot(data=plot_data, aes(x=conditions, y=de_genes)) + geom_point(aes(color=conditions), size=6) + facet_grid(. ~ cell_type)  + scale_color_manual(name = 'condition\ncombination', values = unlist(get_color_coding_dict()[unique(plot_data$conditions)])) +
      xlab('condition combination') + 
      ylab('number of DE genes') +
      theme(#axis.text.x=element_blank(), 
        #axis.ticks = element_blank(), 
        legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(size=12)) + ylim(ylims)
  }
  if(!is.null(title)){
    p <- p + ggtitle(title)
  }
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
  if(colour_top_label){
    striprt <- which( grepl('strip-r', p$layout$name) | grepl('strip-t', p$layout$name) )
    fills <- unlist(get_color_coding_dict()[unique(plot_data$cell_type)])
    k <- 1
    for (i in striprt) {
      j <- which(grepl('rect', p$grobs[[i]]$grobs[[1]]$childrenOrder))
      p$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
  return(p)
}


wide_to_long_table <- function(wide_table) {
  # turn into plottable dataframe
  category_combinations_n <- length(unique(colnames(wide_table))) * length(unique(rownames(wide_table)))
  # make new table
  long_table <- data.frame(x = rep(NA, times = category_combinations_n), y = rep(NA, times = category_combinations_n), number = rep(NA, times = category_combinations_n))
  # keep an index
  i_new <- 1
  # check each row
  for (i_row in 1:nrow(wide_table)) {
    # check each column
    for (i_column in 1:ncol(wide_table)) {
      # get the x, the row
      x <- rownames(wide_table)[i_row]
      # get the y, the column
      y <- colnames(wide_table)[i_column]
      # get the value
      number <- wide_table[i_row, i_column]
      # add in the new table
      long_table[i_new, 'x'] <- x
      long_table[i_new, 'y'] <- y
      long_table[i_new, 'number'] <- number
      # update new index
      i_new <- i_new + 1
    }
  }
  return(long_table)
}


plot_de_number_bars <- function(numbers_table, cell_type_column='cell_type', combination_column='combination', number_column='number', direction_column='direction', cell_types_include=NULL, combinations_include=NULL, grid_single=T, use_label_dict=T, pointless=F, legendless=F, paper_style=T, angle_x_labels=F, colour_top_label=T, use_colour_coding_dict=T){
  numbers_table_to_do <- numbers_table
  # filter if requested
  if(!is.null(cell_types_include)){
    numbers_table_to_do <- numbers_table_to_do[numbers_table_to_do[[cell_type_column]] %in% cell_types_include, ]
  }
  if (!is.null(combinations_include)) {
    numbers_table_to_do <- numbers_table_to_do[numbers_table_to_do[[combination_column]] %in% combinations_include, ]
  }
  # use some better labels if requested
  if (use_label_dict) {
    numbers_table_to_do[[combination_column]] <- as.vector(unlist(label_dict()[as.character(numbers_table_to_do[[combination_column]])]))
    numbers_table_to_do[[cell_type_column]] <- as.vector(unlist(label_dict()[as.character(numbers_table_to_do[[cell_type_column]])]))
  }
  # standardize some columns
  numbers_table_to_do[['combination']] <- as.character(numbers_table_to_do[[combination_column]])
  numbers_table_to_do[['cell_type']] <- numbers_table_to_do[[cell_type_column]]
  numbers_table_to_do[['number']] <- numbers_table_to_do[[number_column]]
  numbers_table_to_do[['cell_type_direction']] <- paste(numbers_table_to_do[['cell_type']], numbers_table_to_do[[direction_column]])
  # make sure that up is always before down
  ct_direction_unique <- unique(numbers_table_to_do[, c('cell_type', direction_column, 'cell_type_direction')])
  ct_direction_unique[['direction_numeric']] <- as.numeric(as.factor(ct_direction_unique[[direction_column]]))
  numbers_table_to_do[['cell_type_direction']] <- factor(numbers_table_to_do[['cell_type_direction']], levels = ct_direction_unique[order(ct_direction_unique[['cell_type']], ct_direction_unique[['direction_numeric']]), 'cell_type_direction'])

  # get the y limits
  ct_condition <- unique(numbers_table_to_do[, c('cell_type', 'combination')])
  ct_condition[['total']] <- apply(ct_condition, 1, function(x) {
    return(sum(numbers_table_to_do[numbers_table_to_do[['cell_type']] == x['cell_type'] & numbers_table_to_do[['combination']] == x['combination'], 'number']))
  })
  # get the max number per celltype
  max_per_ct <- data.frame(cell_type = unique(numbers_table_to_do[['cell_type']]))
  ylims <- c(0, (max(ct_condition[['total']])*1.1))
  
  # initalize
  p <- NULL
  
  # split if requested
  if (!grid_single) {
    p <- ggplot(data = numbers_table_to_do, mapping = aes(x = combination, y = number, fill = cell_type_direction)) + facet_grid(. ~ cell_type)
  }
  else {
    p <- ggplot(data = numbers_table_to_do, mapping = aes(x = cell_type, y = number, fill = cell_type_direction))
  }
  # add to the plot
  p <- p + 
    geom_bar(stat = 'identity', position = 'stack') + 
    ylim(ylims)
  # add colours if requested
  if (use_colour_coding_dict) {
    p <- p + scale_fill_manual(name = paste(cell_type_column, 'direction'), values = get_color_coding_dict())
  }
  else {
    # p <- p + scale_fill_manual(name = paste(cell_type_column, 'direction'))
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
  if(colour_top_label){
    striprt <- which( grepl('strip-r', p$layout$name) | grepl('strip-t', p$layout$name) )
    fills <- unlist(get_color_coding_dict()[unique(numbers_table_to_do$cell_type)])
    k <- 1
    for (i in striprt) {
      j <- which(grepl('rect', p$grobs[[i]]$grobs[[1]]$childrenOrder))
      p$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
      k <- k+1
    }
  }
  return(p)
  
}



####################
# Settings.        #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 500 * 1000 * 1024^2)

# set seed
set.seed(7777)

####################
# Main Code        #
####################

# get the locations of the DE output
limma_output_mo_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/'
limma_output_1m_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation_1m/'

# get a table of the number of DE genes
de_numbers_table_mo_limma_up <- de_genes_number_to_table(limma_output_mo_loc, de_method = 'limma', pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_positive = T, stims = c('24hCA', 'UT'))
de_numbers_table_mo_limma_up[is.na(de_numbers_table_mo_limma_up)] <- 0
de_numbers_table_mo_limma_down <- de_genes_number_to_table(limma_output_mo_loc, de_method = 'limma', pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_negative = T, stims = c('24hCA', 'UT'))
de_numbers_table_mo_limma_down[is.na(de_numbers_table_mo_limma_down)] <- 0
# turn into long tables
de_numbers_table_mo_limma_up_wide <- wide_to_long_table(de_numbers_table_mo_limma_up)
de_numbers_table_mo_limma_down_wide <- wide_to_long_table(de_numbers_table_mo_limma_down)
# set column names
colnames(de_numbers_table_mo_limma_up_wide) <- c('cell_type', 'combination', 'number')
colnames(de_numbers_table_mo_limma_down_wide) <- c('cell_type', 'combination', 'number')
# add new column
de_numbers_table_mo_limma_up_wide[['direction']] <- 'up'
de_numbers_table_mo_limma_down_wide[['direction']] <- 'down'
# merge
de_numbers_table_mo_limma_wide <- rbind(de_numbers_table_mo_limma_down_wide, de_numbers_table_mo_limma_up_wide)
de_numbers_table_mo_limma_wide <- de_numbers_table_mo_limma_wide[de_numbers_table_mo_limma_wide$combination == 'UT24hCA', ]

# get a table of the number of DE genes
de_numbers_table_1m_limma_up <- de_genes_number_to_table(limma_output_1m_loc, de_method = 'limma', pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_positive = T, stims = c('24hCA', 'UT'), condition_name = 'timepoint')
de_numbers_table_1m_limma_up[is.na(de_numbers_table_1m_limma_up)] <- 0
de_numbers_table_1m_limma_down <- de_genes_number_to_table(limma_output_1m_loc, de_method = 'limma', pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_negative = T, stims = c('24hCA', 'UT'), condition_name = 'timepoint')
de_numbers_table_1m_limma_down[is.na(de_numbers_table_1m_limma_down)] <- 0
# turn into long tables
de_numbers_table_1m_limma_up_wide <- wide_to_long_table(de_numbers_table_1m_limma_up)
de_numbers_table_1m_limma_down_wide <- wide_to_long_table(de_numbers_table_1m_limma_down)
# set column names
colnames(de_numbers_table_1m_limma_up_wide) <- c('cell_type', 'combination', 'number')
colnames(de_numbers_table_1m_limma_down_wide) <- c('cell_type', 'combination', 'number')
# add new column
de_numbers_table_1m_limma_up_wide[['direction']] <- 'up'
de_numbers_table_1m_limma_down_wide[['direction']] <- 'down'
# merge
de_numbers_table_1m_limma_wide <- rbind(de_numbers_table_1m_limma_down_wide, de_numbers_table_1m_limma_up_wide)
de_numbers_table_1m_limma_wide <- de_numbers_table_1m_limma_wide[de_numbers_table_1m_limma_wide$combination == 'UT24hCA', ]

# add the dataset
de_numbers_table_mo_limma_wide[['dataset']] <- 'mo'
de_numbers_table_1m_limma_wide[['dataset']] <- '1M'
# merge
de_numbers_table <- rbind(de_numbers_table_mo_limma_wide, de_numbers_table_1m_limma_wide)
plot_de_number_bars(de_numbers_table, combination_column = 'dataset', grid_single = F, use_label_dict = T, use_colour_coding_dict = T, angle_x_labels = T) + ylab('Number of DE genes') + xlab('cell type') + guides(fill = guide_legend(title = 'cell type direction'))

# just get the DE genes
de_mo <- get_de_genes(limma_output_mo_loc, de_method = 'limma', pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_positive = T, stims = c('24hCA', 'UT'), condition_name = 'condition_final')
de_1m <- get_de_genes(limma_output_1m_loc, de_method = 'limma', pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_positive = T, stims = c('24hCA', 'UT'), condition_name = 'timepoint')

# do the ggvenn diagrams
plot_grid(
  ggvenn(data = list('multiome' = de_mo$B$UT24hCA, 'NC2022' = de_1m$B$UT24hCA)),
  ggvenn(data = list('multiome' = de_mo$CD4T$UT24hCA, 'NC2022' = de_1m$CD4T$UT24hCA)),
  ggvenn(data = list('multiome' = de_mo$CD8T$UT24hCA, 'NC2022' = de_1m$CD8T$UT24hCA)),
  ggvenn(data = list('multiome' = de_mo$DC$UT24hCA, 'NC2022' = de_1m$DC$UT24hCA)),
  ggvenn(data = list('multiome' = de_mo$monocyte$UT24hCA, 'NC2022' = de_1m$monocyte$UT24hCA)),
  ggvenn(data = list('multiome' = de_mo$NK$UT24hCA, 'NC2022' = de_1m$NK$UT24hCA)),
  nrow = 3,
  ncol = 2,
  labels = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
)
