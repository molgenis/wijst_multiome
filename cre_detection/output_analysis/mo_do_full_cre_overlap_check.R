#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_do_full_cre_overlap_check.R
# Function: compare SCENIC+, eQTL+caQTL, pseudobulk and binomial CRE detection methods
#
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(mdfiver)
library(ggplot2)
library(cowplot)
library(UpSetR)
library(roycols)


####################
# Functions        #
####################


read_pseudobulk_cre_output_per_celltype <- function(pseudobulk_output_folder, cell_types=NULL, filename_output='qtl_results_all.txt.gz', significance_column='empirical_feature_p_value', significance_cutoff=0.05) {
  # list all the files in the directory
  cell_type_folders <- list.dirs(pseudobulk_output_folder, full.names = F, recursive = F)
  # intersect the cell type folders with the cell types we are interested in
  if (!is.null(cell_types)) {
    cell_type_folders <- intersect(cell_type_folders, cell_types)
  }
  # save the results in a list
  output_per_celltype <- list()
  # now check each cell type
  for (cell_type in cell_type_folders) {
    # paste together the full file path
    cell_type_output_loc <- paste(pseudobulk_output_folder, cell_type, filename_output, sep = '/')
    # read this file
    cell_type_output <- fread(cell_type_output_loc, header = T, sep = '\t')
    # filter the file if requested
    if (!is.null(significance_column)) {
      cell_type_output <- cell_type_output[
        cell_type_output[[significance_column]] < significance_cutoff, 
      ]
    }
    # add the cell type
    cell_type_output[['cell_type']] <- cell_type
    # put in the list
    output_per_celltype[[cell_type]] <- cell_type_output
  }
  return(output_per_celltype)
}


read_binomial_output_per_celltype <- function(binomial_output_folder, cell_types=NULL, filename_output='meta_result.tsv.gz', significance_column='meta_q', significance_cutoff=0.05) {
  # just use the pseudobulk function
  output_per_celltype <- read_pseudobulk_cre_output_per_celltype(binomial_output_folder, cell_types = cell_types, filename_output = filename_output, significance_column = significance_column, significance_cutoff = significance_cutoff)
  return(output_per_celltype)
}


plot_concondance <- function(dataset_to_compare, d1_effect_column='d1_zscore', d2_effect_column='d2_zscore') {
  # get the minimal significant z for k1
  min_sig_z_d1 <- min(abs(dataset_to_compare[[d1_effect_column]]))
  min_sig_z_d2 <- min(abs(dataset_to_compare[[d2_effect_column]]))
  # and the max z
  max_sig_z_d1 <- max(abs(dataset_to_compare[[d1_effect_column]]))
  max_sig_z_d2 <- max(abs(dataset_to_compare[[d2_effect_column]]))
  
  # calculate concordance
  n_significant_both <- nrow(dataset_to_compare)
  n_significant_directed_both <- sum(sign(dataset_to_compare[[d1_effect_column]]) == sign(dataset_to_compare[[d2_effect_column]]))
  concordance <- round(n_significant_directed_both / n_significant_both, digits = 3)
  
  # plot the concordance of the two
  p <- ggplot(data = dataset_to_compare, mapping = aes(x = !! rlang::sym(d1_effect_column), y = !! rlang::sym(d2_effect_column))) + 
    geom_point(size = 0.2) + 
    xlim(c(max_sig_z_d1 * -1.1, max_sig_z_d1 * 1.1)) +
    ylim(c(max_sig_z_d2 * -1.1, max_sig_z_d2 * 1.1)) + 
    # left to right block of non-significant effects
    geom_rect(aes(xmin = -1 * max_sig_z_d1, xmax = max_sig_z_d1, ymin = -1 * min_sig_z_d2, ymax = min_sig_z_d2), 
              fill = "white", alpha = 0.005) +
    # bottom to top block of non-significant effects
    geom_rect(aes(xmin = -1 * min_sig_z_d1, xmax = min_sig_z_d1, ymin = -1 * max_sig_z_d2, ymax = max_sig_z_d2), 
              fill = "white", alpha = 0.005) +
    # bottom left block
    geom_rect(aes(xmin = -1 * max_sig_z_d1, xmax = -1 *min_sig_z_d1, ymin = -1 * max_sig_z_d2, ymax = -1 * min_sig_z_d2), 
              fill = "#0072B2", alpha = 0.01) +
    # bottom right block
    geom_rect(aes(xmin = min_sig_z_d1, xmax = max_sig_z_d1, ymin = -1 * max_sig_z_d2, ymax = -1 * min_sig_z_d2), 
              fill = "#D55E00", alpha = 0.01) +
    # top left block
    geom_rect(aes(xmin = -1 * max_sig_z_d1, xmax = -1 *min_sig_z_d1, ymin = min_sig_z_d2, ymax = max_sig_z_d2), 
              fill = "#D55E00", alpha = 0.01) + 
    # top right block
    geom_rect(aes(xmin = min_sig_z_d1, xmax = max_sig_z_d1, ymin = max_sig_z_d2, ymax = min_sig_z_d2), 
              fill = "#0072B2", alpha = 0.01) + 
    # labels for x and y axis
    xlab('effect in dataset1') +
    ylab('effect in dataset2') +
    # vertical negative z ccombinedoff line
    geom_vline(xintercept=c(-1 *min_sig_z_d1), color="black", size=0.5, linetype = "dashed") +
    # horizonal negative z ccombinedoff line
    geom_hline(yintercept=c(-1 *min_sig_z_d2), color="black", size=0.5, linetype="dashed") +
    # vertical positive z ccombinedoff line
    geom_vline(xintercept=c(min_sig_z_d1), color="black", size=0.5, linetype = "dashed") +
    # horizonal positive z ccombinedoff line
    geom_hline(yintercept=c(min_sig_z_d2), color="black", size=0.5, linetype="dashed") +
    # make backgrounds white
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
    # add the concordance
    annotate("label", x = max_sig_z_d1 * 0.75 , y = max_sig_z_d2 * -0.75, label = paste('concordance', concordance, sep = ':\n')) +
    # add the names of the concordant and non-concordant blocks
    annotate("text", x = max_sig_z_d1 * -0.70 , y = max_sig_z_d2 * 0.75, label = 'disconcordant', colour = '#D55E00', fontface = 'bold') +
    # add the names of the concordant and non-concordant blocks
    annotate("text", x = max_sig_z_d1 * 0.70 , y = max_sig_z_d2 * 0.75, label = 'concordant', colour = '#0072B2', fontface = 'bold') +
    # add the title
    ggtitle('Concordance of effects between datasets')
  return(p)
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


#' Plot Sharing of  Genes per Cell Type
#'
#' This function plots the sharing of differentially genes across cell types using an UpSet plot. It allows the use of custom label and color dictionaries.
#'
#' @param genes_per_ct list with the genes for each cell type
#' @param use_label_dict A logical value indicating whether to use a custom label dictionary for renaming cell types. Default is TRUE.
#' @param use_color_dict A logical value indicating whether to use a custom color dictionary for cell types. Default is TRUE.
#' @param n_intersects value describing how many intersections to plot, default is all
#' @return An UpSet plot showing the sharing of DE genes across cell types.
#'
plot_sharing_per_celltype <- function(genes_per_ct, use_label_dict=T, use_color_dict=T, n_intersects=NA){
  # get the total overlap
  total_overlap <- fromList(genes_per_ct)
  # count how many in each overlap
  combination_to_overlap <- list()
  apply(total_overlap, 1, function(x){
    # get which ones have a number
    positions_1 <- x > 0
    # get which these are, and paste together
    cols_1 <- paste(colnames(total_overlap)[positions_1], collapse = ',')
    # add to list if does not exist
    if (cols_1 %in% c(names(combination_to_overlap))) {
      combination_to_overlap[[cols_1]] <<- combination_to_overlap[[cols_1]] + 1
    }
    else{
      combination_to_overlap[[cols_1]] <<- 1
    }
  })
  # make into an ordered df
  combination_numbers <- data.frame('combination' = names(combination_to_overlap), 'nr' = as.vector(unlist(combination_to_overlap)))
  combination_numbers <- combination_numbers[order(combination_numbers[['nr']], combination_numbers[['combination']], decreasing = T), ]
  # get lowest number we have for the sets
  intersect_size_smallest <- 0
  if (is.na(n_intersects)) {
    intersect_size_smallest <- min(combination_numbers[['nr']])
  }
  else {
    intersect_size_smallest <- combination_numbers[n_intersects, 'nr']
  }
  # set queries and colours
  queries <- list()
  sets.bar.color <- 'black'
  if(use_color_dict){
    # create df to store the number of each set, so we know how to order
    nrs_df <- NULL
    # get the cell types we have
    cell_types <- names(genes_per_ct)
    # get colour codes for the cell types
    cell_type_colours <- get_color_coding_dict()
    # also get colours for the cell types we don't have colours for
    cell_type_colours_missing <- roycols::get_color_list(setdiff(cell_types, names(cell_type_colours)))
    # and merge them
    cell_type_colours <- c(cell_type_colours, cell_type_colours_missing)
    # add the colors for the cell types
    i <- 1
    for(cell_type in cell_types){
      # check if there is a singleton for this cell type, so this celltype is 1, but the total of the row is also 1
      n_singletons <- nrow(total_overlap[total_overlap[[cell_type]] == 1 & rowSums(total_overlap) == 1, ])
      # if we have singletons, colour it
      if (n_singletons > 0 & n_singletons >= intersect_size_smallest) {
        # add for the singles in the intersection sizes
        ct_list <- list(
          query = intersects,
          params = list(cell_type),
          color = cell_type_colours[[cell_type]],
          active = T)
        queries[[i]] <- ct_list
        i <- i + 1
      }
      # add for the DF to order the set sizes
      numbers_row <- data.frame(ct=c(cell_type), nr=c(length(genes_per_ct[[cell_type]])), stringsAsFactors = F)
      if(is.null(nrs_df)){
        nrs_df <- numbers_row
      }
      else{
        nrs_df <- rbind(nrs_df, numbers_row)
      }
    }
    # get the order of the sets
    ordered_cts <- nrs_df[order(nrs_df$nr, decreasing = T), 'ct']
    # add the colors for the sets
    sets.bar.color <- unlist(cell_type_colours[ordered_cts])
    # make the plot
    upset(total_overlap, order.by = 'freq', nsets = length(genes_per_ct), queries = queries, sets.bar.color=sets.bar.color, nintersects = n_intersects)
  }
  else {
    upset(total_overlap, order.by = 'freq', nsets = length(genes_per_ct), nintersects = n_intersects)
  }
  
  #return(DE_genes_per_ct)
}

create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL, angle_labels=T, premade_table=F, freq_column='freq'){
  confusion_table <- NULL
  if (!premade_table) {
    confusion_table <- create_confusion_table(assignment_table, truth_column, prediction_column)
    # round the frequency off to a sensible cutoff
    confusion_table$freq <- round(confusion_table$freq, digits=2)
  }
  else {
    confusion_table <- data.frame(
      'truth' = assignment_table[[truth_column]], 
      'prediction' = assignment_table[[prediction_column]], 
      'freq' = assignment_table[[freq_column]]
    )
  }

  # turn into plot
  p <- ggplot(data=confusion_table, aes(x=truth, y=prediction, fill=freq)) + geom_tile() + scale_fill_gradient(low='blue', high='red') + geom_text(aes(label=freq))
  # some options
  if(!is.null(truth_column_label)){
    p <- p + xlab(truth_column_label)
  }
  if(!is.null(prediction_column_label)){
    p <- p + ylab(prediction_column_label)
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p)
}


####################
# Settings         #
####################

# luck seed
set.seed(7777)


####################
# Main code        #
####################

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# location of the overlapping caQTLs and eQTLs
qtl_overlap_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl_caqtl_overlap/combined/L1/all/eqtl_caqtl_overlapping_variants.tsv.gz'
# location of the pseudobulk CRE mapping
pseudobulk_output_folder <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/eQTA/L1/'
# location of the binomial method
binomial_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/combined/betas_ps/'

# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# read the qtl overlap
qtl_overlap <- fread(qtl_overlap_loc, header = T, sep = '\t')

# read the pseudobulk outputs
pseudobulk_output_ut_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, 'UT', sep = '/'))
pseudobulk_output_24hca_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, '24hCA', sep = '/'))
# merge the cell types
pseudobulk_output_ut <- do.call('rbind', pseudobulk_output_ut_list)
pseudobulk_output_24hca <- do.call('rbind', pseudobulk_output_24hca_list)
# add the condition
pseudobulk_output_ut[['condition']] <- 'UT'
pseudobulk_output_24hca[['condition']] <- '24hCA'
# merge them
pseudobulk_output <- do.call('rbind', list(pseudobulk_output_ut, pseudobulk_output_24hca))

# read the binomial results
binomial_output_list <- read_binomial_output_per_celltype(binomial_output_loc)
# merge them
binomial_output <- do.call('rbind', binomial_output_list)

# add the Z to the pseudobulk analysis
pseudobulk_output[['zscore']] <- pseudobulk_output[['beta']] / pseudobulk_output[['beta_se']]
# and the qtl overlap
qtl_overlap[['z_caqtl']] <- qtl_overlap[['beta_caqtl']] / qtl_overlap[['se_caqtl']]
qtl_overlap[['z_eqtl']] <- qtl_overlap[['beta_eqtl']] / qtl_overlap[['se_eqtl']]
# get the sign overlap
qtl_overlap[['sign']] <- sign(qtl_overlap[['z_caqtl']]) * sign(qtl_overlap[['z_eqtl']])


# sort all of them by the Z
pseudobulk_output <- pseudobulk_output[order(abs(pseudobulk_output[['zscore']])), ]
binomial_output <- binomial_output[order(abs(binomial_output[['meta_z']])), ]
scenic_output <- scenic_output[order(scenic_output[['rho_R2G']]), ]
qtl_overlap <- qtl_overlap[order(qtl_overlap[['z_eqtl']], qtl_overlap[['z_caqtl']]), ]

# remove the entries that are more likely to be false positives
scenic_output_unique <- scenic_output[scenic_output[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]

# get unique ones
pseudobulk_output_unique <- pseudobulk_output[!duplicated(paste(pseudobulk_output[['snp_id']], pseudobulk_output[['feature_id']])), ]
binomial_output_unique <- binomial_output[!duplicated(paste(binomial_output[['region']], binomial_output[['gene']])), ]
scenic_output_unique <- scenic_output_unique[!duplicated(paste(scenic_output_unique[['Region']], scenic_output_unique[['Gene']])), ]
scenic_output_unique_unfiltered <- scenic_output[!duplicated(paste(scenic_output[['Region']], scenic_output[['Gene']])), ]
qtl_overlap_unique <- qtl_overlap[!duplicated(paste(qtl_overlap[['feature_caqtl']], qtl_overlap[['feature_eqtl']])), ]

# show how many we have in each set
nrow(binomial_output_unique)
# [1] 12710
nrow(scenic_output_unique_unfiltered)
# [1] 80440
nrow(scenic_output_unique)
# [1] 53796
nrow(pseudobulk_output_unique)
# [1] 121935
nrow(qtl_overlap_unique)
# [1] 7677

# add region to gene column
pseudobulk_output_unique[['r2g']] <- paste(pseudobulk_output_unique[['snp_id']], pseudobulk_output_unique[['feature_id']])
binomial_output_unique[['r2g']] <- paste(binomial_output_unique[['region']], binomial_output_unique[['gene']])
scenic_output_unique[['r2g']] <- paste(gsub(':', '-', scenic_output_unique[['Region']]), scenic_output_unique[['Gene']])
scenic_output_unique_unfiltered[['r2g']] <- paste(gsub(':', '-', scenic_output_unique_unfiltered[['Region']]), scenic_output_unique_unfiltered[['Gene']])
qtl_overlap_unique[['r2g']] <- paste(qtl_overlap_unique[['feature_caqtl']], qtl_overlap_unique[['feature_eqtl']])

# plot these numbers
plot_sharing_per_celltype(
  list('pseudobulk' = pseudobulk_output_unique[['r2g']], 
       'binomial' = binomial_output_unique[['r2g']], 
       'scenic' = scenic_output_unique[['r2g']], 
       'scenic unfiltered' = scenic_output_unique_unfiltered[['r2g']],
       'qtl' = qtl_overlap_unique[['r2g']]), 
  use_label_dict=F, use_color_dict=T
)


# get the percentage positive
nrow(pseudobulk_output_unique[sign(pseudobulk_output_unique[['zscore']]) == 1, ]) / nrow(pseudobulk_output_unique)
# [1] 0.8557264
nrow(binomial_output_unique[sign(binomial_output_unique[['meta_z']]) == 1, ]) / nrow(binomial_output_unique)
# [1] 0.9983478
nrow(scenic_output_unique[sign(scenic_output_unique[['rho_R2G']]) == 1, ]) / nrow(scenic_output_unique)
# [1] 0.6687718 / 1
nrow(qtl_overlap_unique[sign(qtl_overlap_unique[['sign']]) == 1, ]) / nrow(qtl_overlap_unique)
# [1] 0.7630585

# plot these numbers
n_pos_tbl <- data.frame(
  'method' = c('pseudobulk', 'binomial', 'scenic', 'QTL', 'pseudobulk', 'binomial', 'scenic', 'QTL'), 
  'direction' = c('positive', 'positive', 'positive', 'positive', 'negative', 'negative', 'negative', 'negative'), 
  'n' = c(0.8557264, 0.9983478, 0.6687718, 0.7630585, 1-0.8557264, 1-0.9983478, 1-0.6687718, 1-0.7630585)
)
ggplot(data = n_pos_tbl, mapping = aes(x = method, y = n, fill = direction)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('CRE detection method') +
  ylab('fraction') +
  ggtitle('CRE method result directions') +
  scale_fill_manual(values = list('positive' = 'darkblue', 'negative' = 'darkred')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# merge pseudobulk and binominal
pseudobulk_vs_binomial <- merge(x = pseudobulk_output_unique, y = binomial_output_unique, by = 'r2g')
# check pseudobulk and scenic
pseudobulk_vs_scenic <- merge(x = pseudobulk_output_unique, y = scenic_output_unique, by = 'r2g')
# check pseudobulk and qtls
pseudobulk_vs_qtl <- merge(x = pseudobulk_output_unique, y = qtl_overlap_unique, by = 'r2g')
# and binomial vs scenic
scenic_vs_binomial <- merge(x = scenic_output_unique, y = binomial_output_unique, by = 'r2g')
# scenic vs qtl
scenic_vs_qtl <- merge(x = scenic_output_unique, y = qtl_overlap_unique, by = 'r2g')
# binomial vs qtl
binomial_vs_qtl <- merge(x = binomial_output_unique, y = qtl_overlap_unique, by = 'r2g')


# get the concordances
nrow(pseudobulk_vs_binomial[sign(pseudobulk_vs_binomial[['zscore']]) == sign(pseudobulk_vs_binomial[['meta_z']]), ]) / nrow(pseudobulk_vs_binomial)
# [1] 0.9458333
nrow(pseudobulk_vs_scenic[sign(pseudobulk_vs_scenic[['zscore']]) == sign(pseudobulk_vs_scenic[['rho_R2G']]), ]) / nrow(pseudobulk_vs_scenic)
# [1] 0.8541311
nrow(pseudobulk_vs_qtl[sign(pseudobulk_vs_qtl[['zscore']]) == sign(pseudobulk_vs_qtl[['sign']]), ]) / nrow(pseudobulk_vs_qtl)
# [1] 0.9856528
nrow(scenic_vs_binomial[sign(scenic_vs_binomial[['meta_z']]) == sign(scenic_vs_binomial[['rho_R2G']]), ]) / nrow(scenic_vs_binomial)
# [1] 0.8156997 / 1
nrow(scenic_vs_qtl[sign(scenic_vs_qtl[['sign']]) == sign(scenic_vs_qtl[['rho_R2G']]), ]) / nrow(scenic_vs_qtl)
# [1] 0.7747748 / 0.8636364
nrow(binomial_vs_qtl[sign(binomial_vs_qtl[['meta_z']]) == sign(binomial_vs_qtl[['sign']]), ]) / nrow(binomial_vs_qtl)
# [1] 0.7746358

# in a pairwise table
pairwise_correlation_table <- data.frame(
  'method1' = c('pseudobulk', 'pseudobulk', 'pseudobulk', 'pseudobulk', 
                'scenic', 'scenic', 'scenic', 'scenic', 
                'binomial', 'binomial', 'binomial', 'binomial',
                'QTL', 'QTL', 'QTL', 'QTL'), 
  'method2' = c('pseudobulk', 'scenic', 'binomial', 'QTL', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 
                'pseudobulk', 'scenic', 'binomial', 'QTL'), 
  'correlation' = c(1, 0.85, 0.82, 0.99, 
                    0.85, 1, 0.81, 0.77, 
                    0.95, 0.81, 1, 0.77, 
                    0.99, 0.77, 0.77, 1)
)
create_confusion_matrix(pairwise_correlation_table, truth_column = 'method1', prediction_column = 'method2', freq_column = 'correlation', premade_table = T, truth_column_label = 'method 1', prediction_column_label = 'method 2') +
  ggtitle('correlations of effect sizes\nin CRE detection methods')

nrow(pseudobulk_vs_binomial)
# [1] 960
nrow(pseudobulk_vs_scenic)
# [1] 1755 / 1509
nrow(pseudobulk_vs_qtl)
# [1] 1394
nrow(scenic_vs_binomial)
# [1] 586 / 478
nrow(scenic_vs_qtl)
# [1] 666 / 572
nrow(binomial_vs_qtl)
# [1] 3501

# plot them as well
plot_grid(
  plot_concondance(pseudobulk_vs_binomial, 'zscore', 'meta_z') + ggtitle('Effects of pseudobulk vs binomial\nCRE detection') + xlab('Pseudobulk Z-score') + ylab('Binomial model Z-score'), 
  plot_concondance(pseudobulk_vs_scenic, 'zscore', 'rho_R2G') + ggtitle('Effects of pseudobulk vs SCENIC+ CRE\ndetection') + xlab('Pseudobulk Z-score') + ylab('SCENIC+ R2G Rho'), 
  plot_concondance(scenic_vs_binomial, 'rho_R2G', 'meta_z') + ggtitle('Effects of SCENIC+ vs binomial\nCRE detection') + xlab('SCENIC+ R2G Rho') + ylab('Binomial model Z-score')
)

