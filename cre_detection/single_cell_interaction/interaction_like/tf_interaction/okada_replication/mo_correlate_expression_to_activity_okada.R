#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_correlate_expression_to_activity_okada.R
# Function: correlate TF expression to inferred TF activity in the Okada dataset
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(Seurat)
library(progress)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggpubr)


####################
# Functions        #
####################

normalize_mj <- function(seurat_object) {
  # get the count matrix where we have the correct cell type
  count_matrix <- GetAssayData(seurat_object, slot = "counts")
  # ignore genes that are never expressed
  count_matrix <-  count_matrix[which(rowSums(count_matrix) != 0), ]
  # create new object to store the counts in
  norm_count_matrix <- count_matrix
  # do mean sample-sum normalization
  sample_sum_info = colSums(norm_count_matrix)
  mean_sample_sum = mean(sample_sum_info)
  sample_scale = sample_sum_info / mean_sample_sum
  # divide each column by sample_scale
  norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
  if ('layers' %in% slotNames(seurat_object[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    seurat_object[['MJ']] <- CreateAssay5Object(data = norm_count_matrix)
    
  } else {
    print('using Seurat v3/4 style \'slot\'')
    seurat_object[['MJ']] <- CreateAssayObject(data = norm_count_matrix)
  }
  return(seurat_object)
}


correlate_ereg_to_tf <- function(seurat_object, mapping, expression_assay='MJ', activity_assay='tf_activity', expression_slot='data', activity_slot='data') {
  # get the number of eregulons
  n_tf <- nrow(mapping)
  # store result
  cor_tf_ereg <- data.frame(
    matrix(, nrow = n_tf, ncol = 3, dimnames = list(mapping$eRegulon, c('tf', 'eregulon', 'correlation')))
  )
  # make a progress bar
  pb <- progress_bar$new(total = n_tf)
  # initialize the progress bar
  pb$tick(0)
  # check each eregulon
  for (i in 1:n_tf) {
    # update the progress bar
    pb$tick()
    # extract the TF name and the eRegulon name
    tf <- mapping$TF[i]
    ereg <- mapping$eRegulon[i]
    # set these in the matrix
    cor_tf_ereg[['tf']][i] <- tf
    cor_tf_ereg[['eregulon']][i] <- ereg
    # check if there is a TF in the MJ assay that matches the TF in the AUC assay
    if (tf %in% rownames(seurat_object@assays[[expression_assay]]@features)) {
      # get the index of the TF
      tf_index <- which(rownames(seurat_object@assays[[expression_assay]]@features) == tf)
      # get the index of the eregulon
      ereg_index <- which(rownames(seurat_object@assays[[activity_assay]]@features) == ereg)
      # get the expression and activity values
      expr_values <- seurat_object@assays[[expression_assay]]@layers[[expression_slot]][tf_index, ]
      activity_values <- seurat_object@assays[[activity_assay]]@layers[[activity_slot]][ereg_index, ]
      # calculate the correlation
      cor_value <- cor(expr_values, activity_values, method = 'spearman')
      # put in the matrix
      cor_tf_ereg[['correlation']][i] <- cor_value
    } else {
      cor_tf_ereg[['correlation']][i] <- NA
    }
  }
  return(cor_tf_ereg)
}


correlate_ereg_to_tf_v4 <- function(seurat_object, mapping, expression_assay='MJ', activity_assay='tf_activity') {
  # get the number of eregulons
  n_tf <- nrow(mapping)
  # store result
  cor_tf_ereg <- data.frame(
    matrix(, nrow = n_tf, ncol = 3, dimnames = list(mapping$eRegulon, c('tf', 'eregulon', 'correlation')))
  )
  # make a progress bar
  pb <- progress_bar$new(total = n_tf)
  # initialize the progress bar
  pb$tick(0)
  # check each eregulon
  for (i in 1:n_tf) {
    # update the progress bar
    pb$tick()
    # extract the TF name and the eRegulon name
    tf <- mapping$TF[i]
    ereg <- mapping$eRegulon[i]
    # set these in the matrix
    cor_tf_ereg[['tf']][i] <- tf
    cor_tf_ereg[['eregulon']][i] <- ereg
    # check if there is a TF in the MJ assay that matches the TF in the AUC assay
    if (tf %in% rownames(seurat_object@assays[[expression_assay]]@data)) {
      # get the index of the TF
      tf_index <- which(rownames(seurat_object@assays[[expression_assay]]@data) == tf)
      # get the index of the eregulon
      ereg_index <- which(rownames(seurat_object@assays[[activity_assay]]@data) == ereg)
      # get the expression and activity values
      expr_values <- seurat_object@assays[[expression_assay]]@data[tf_index, ]
      activity_values <- seurat_object@assays[[activity_assay]]@data[ereg_index, ]
      # calculate the correlation
      cor_value <- cor(expr_values, activity_values, method = 'spearman')
      # put in the matrix
      cor_tf_ereg[['correlation']][i] <- cor_value
    } else {
      cor_tf_ereg[['correlation']][i] <- NA
    }
  }
  return(cor_tf_ereg)
}

get_color_coding_dict <- function() {
  # medhigh
  color_coding_dict <- list()
  color_coding_dict[["B"]] <- "#71BC4B"
  #color_coding_dict[['CD4_T_cells']] <- '#7FC97F'
  color_coding_dict[['CD4_T_cells']] <- '#153057'
  color_coding_dict[['CD4T']] <- '#153057'
  #color_coding_dict[['CD8_T_cells']] <- '#BEAED4'
  color_coding_dict[['CD8_T_cells']] <- '#009DDB'
  color_coding_dict[['CD8T']] <- '#009DDB'
  #color_coding_dict[['Dendritic_cells']] <- '#FDC086'
  color_coding_dict[['Dendritic_cells']] <- '#965EC8'
  color_coding_dict[['DC']] <- '#965EC8'
  color_coding_dict[['Endothelial_cells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia_cells']] <- '#F0027F'
  color_coding_dict[['Mast_cells']] <- '#BF5B17'
  color_coding_dict[['Mature_absorptive_enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature_secretory_enterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory_B']] <- '#D95F02'
  color_coding_dict[['Microfold_cell']] <- '#BEAED4'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Monocyte']] <- '#EDBA1B'
  color_coding_dict[['Naive_B_cells']] <- '#FDC086'
  color_coding_dict[['NK']] <- '#E64B50'
  #color_coding_dict[['Plasma_cells']] <- '#E7298A'
  color_coding_dict[['Plasma_cells']] <- '#DB8E00'
  color_coding_dict[['Stem_cells']] <- '#66A61E'
  color_coding_dict[['Stromal_cells']] <- '#8DD3C7'
  #color_coding_dict[['T_others']] <- '#A6761D'
  color_coding_dict[['T_others']] <- '#FF63B6'
  color_coding_dict[['Transit_amplifying_cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  #color_coding_dict[['CD4+ T cells']] <- '#7FC97F'
  color_coding_dict[['CD4+ T cells']] <- '#153057'
  color_coding_dict[['CD4+ T']] <- '#153057'
  #color_coding_dict[['CD8+ T cells']] <- '#BEAED4'
  color_coding_dict[['CD8+ T cells']] <- '#009DDB'
  color_coding_dict[['CD8+ T']] <- '#009DDB'
  #color_coding_dict[['Dendritic cells']] <- '#FDC086'
  color_coding_dict[['Dendritic cells']] <- '#965EC8'
  color_coding_dict[['Endothelial cells']] <- '#FFFFB3'
  color_coding_dict[['Endothelial\ncells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia cells']] <- '#F0027F'
  color_coding_dict[['MAST cells']] <- '#BF5B17'
  color_coding_dict[['Mature absorptive enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature\nabsorptive\nenterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature secretory enterocytes']] <- '#1B9E77'
  color_coding_dict[['Mature secretory\nenterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory B cells']] <- '#D95F02'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Microfold cells']] <- '#BEAED4'
  color_coding_dict[['Monocytes']] <- '#EDBA1B'
  color_coding_dict[['Naive B cells']] <- '#FDC086'
  #color_coding_dict[['Plasma cells']] <- '#E7298A'
  color_coding_dict[['Plasma cells']] <- '#DB8E00'
  color_coding_dict[['Stem cells']] <- '#66A61E'
  color_coding_dict[['Stromal cells']] <- '#8DD3C7'
  #color_coding_dict[['other T cells']] <- '#A6761D'
  color_coding_dict[['other T cells']] <- '#FF63B6'
  color_coding_dict[['Transit amplifying cells']] <- '#FF7F00'
  color_coding_dict[['Transit\namplifying cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  color_coding_dict[['UT']] <- 'gray'
  color_coding_dict[['24hCA']] <- 'darkgreen'
  # up and down regulation will be added to, we need a whitening percentage
  pct_whitening <- 40
  # then we will check each cell type
  for (cell_type in names(color_coding_dict)) {
    # the up color is the same as the regular one
    color_coding_dict[[paste(cell_type, 'up')]] <- color_coding_dict[[cell_type]]
    # but the down one will have a more faded colour
    color_coding_dict[[paste(cell_type, 'down')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    # we'll do something similiar when we have multiple conditions
    color_coding_dict[[paste(cell_type, 'combined')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'UT')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, '24hCA')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "black"))(100)[pct_whitening]
  }
  # general
  color_coding_dict[['AI']] <- 'darkblue'
  color_coding_dict[['NI']] <- 'darkred'
  color_coding_dict[['Actively Inflamed']] <- 'darkblue'
  color_coding_dict[['Non-Inflamed']] <- 'darkred'
  return(color_coding_dict)
}


################################
# Main code #
################################

# get AUC location
cells_AUC_notop_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/okada//tf_interaction/okada_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed.tsv.gz'
cells_AUC_notop <- read.table(cells_AUC_notop_loc, header = T, sep = '\t')

# create a mapping of the eregulon+iegene to just the eregulon
mapping <- data.frame(
  'eRegulon_iegene' = cells_AUC_notop[[1]], 
  'eRegulon' = gsub('_[^_]+$', '', cells_AUC_notop[[1]])
)
# add the TF based on the eregulon
mapping[['TF']] <- gsub("_(extended|direct).*", "", mapping[['eRegulon']])
# order
mapping <- mapping[order(mapping[['eRegulon']]), ]
# then keep only the first eregulon
mapping_tops <- mapping[!duplicated(mapping[['eRegulon']]), ]
# then filter the auc matrix to only the first eregulon
cells_AUC_notop <- cells_AUC_notop[cells_AUC_notop[[1]] %in% mapping_tops[['eRegulon_iegene']], ]

# turn into a matrix
auc_mtx <- Matrix::as.matrix(cells_AUC_notop[, 2:ncol(cells_AUC_notop)])
# and set the features as rownames instead of as the first column
rownames(auc_mtx) <- cells_AUC_notop[[1]]

# location of the Seurat object
seurat_object_loc <- paste0('/groups/umcg-franke-scrna/tmp02/external_datasets/okada/seurat_objects/okada_major_cts.rds')
# read the object
seurat_object <- readRDS(seurat_object_loc)
# add MJ normalization
seurat_object <- normalize_mj(seurat_object)

# next make into an assay
auc_assay <- CreateAssayObject(data = auc_mtx)
# add the AUC assay to the seurat object
seurat_object[['tf_activity']] <- auc_assay

# make a dataframe to store results
cor_tf_ereg <- correlate_ereg_to_tf_v4(seurat_object, mapping = data.frame('TF' = mapping_tops[['TF']], 'eRegulon' = gsub('_', '-', mapping_tops[['eRegulon_iegene']])))
# add dataset label
cor_tf_ereg[['dataset']] <- 'okada'
# rename to eregulon to eregulon_iegene
colnames(cor_tf_ereg) <- c('tf', 'eregulon_iegene', 'correlation', 'dataset')
# then add actual ereg name
cor_tf_ereg[['eregulon']] <- mapping[match(cor_tf_ereg[['eregulon_iegene']], gsub('_', '-', mapping[['eRegulon_iegene']])), ][['eRegulon']]
# then a nicer name
cor_tf_ereg[['eregulon_nice']] <- gsub('_', ' ', cor_tf_ereg[['eregulon']])
# then plot
p_scenic_tf_expression_activity_correlations_okada <- plot_grid(
  ggplot(
  data = cor_tf_ereg,
  mapping = aes(
    x=dataset, 
    y=eregulon_nice, 
    fill=correlation)
  ) + geom_tile() + scale_fill_gradient2(low='darkblue', mid = 'white', high='darkred', midpoint = 0)  + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
    xlab('Dataset') + ylab('eRegulon') + 
    theme(legend.position = 'none'), 
  get_legend(
    ggplot(data = data.frame(x = c('a', 'a'), y = c('w', 'z'), z = c(-1,1)), mapping = aes(x = x, y = y, fill = z)) + 
      geom_tile() +
      scale_fill_gradient2(low='darkblue', mid = 'white', high='darkred', midpoint = 0) + 
      labs(fill = 'rho')
  ), 
  rel_widths = c(3,1)
)
# show the plot
p_scenic_tf_expression_activity_correlations_okada
# and save the result
ggsave(filename = '~/multiome/plots/mo_scenic_tf_expression_activity_correlations_okada.pdf', plot = p_scenic_tf_expression_activity_correlations_okada, width = 5, height = 18)

# also check the correlation difference between the conditions
p_scenic_tf_expression_activity_correlation_density_okada <- ggplot(
  data = cor_tf_ereg, 
  aes(x = correlation, fill = dataset)
) + geom_density(alpha = .5) + scale_fill_manual(values = list('okada' = 'lightgreen')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('Rho') + ylab('Density') + labs(fill = 'Dataset')
# show the plot
p_scenic_tf_expression_activity_correlation_density_okada
# and save the result
ggsave(filename = '~/multiome/plots/mo_scenic_tf_expression_activity_correlation_density_okada.pdf', plot = p_scenic_tf_expression_activity_correlation_density_okada, width = 6, height = 4)

# save the object temporarily
seurat_object_tfa_loc <- paste0('/groups/umcg-franke-scrna/tmp02/external_datasets/okada/seurat_objects/okada_major_cts_tfa.rds')
saveRDS(seurat_object, seurat_object_tfa_loc)
mdfiver::create_sha256_for_file(seurat_object_tfa_loc)

# save result as well
okada_tfe_tfa_correlations_loc <- '~/multiome/tables/mo_tfe_vs_tfa_correlation_okada.tsv.gz'
write.table(cor_tf_ereg, gzfile(okada_tfe_tfa_correlations_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_sha256_for_file(okada_tfe_tfa_correlations_loc)
