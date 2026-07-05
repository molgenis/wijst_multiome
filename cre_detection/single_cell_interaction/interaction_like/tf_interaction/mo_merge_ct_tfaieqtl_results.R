#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_merge_ct_tfaieqtl_results.R
# Function: merge the per-ct TFa-i-eQTL results
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(ggplot2)
library(UpSetR)
library(roycols)
library(cowplot)

####################
# functions        #
####################

read_all_tfaieqtl_files <- function(results_dir, result_append='/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz', cell_types=NULL) {
  # get all directories in source directory
  ct_dirs <- list.dirs(results_dir, recursive=F, full.names = F)
  # overlap with CTs if supplied
  if (!is.null(cell_types)) {
    ct_dirs <- intersect(ct_dirs, cell_types)
  }
  # create a list with CTs
  res_per_ct <- list()
  # check each ct directory for the results file and read it in
  for(ct in ct_dirs) {
    # construct the results file path
    results_file <- paste(results_dir, ct, result_append, sep = '/')
    # check if the file exists
    if (file.exists(results_file)) {
      # read in the results file
      ct_results <- fread(results_file, header = T, sep = '\t')
      # add a column for the cell type
      ct_results[, cell_type := ct]
      # add to the list
      res_per_ct[[ct]] <- ct_results
    } else {
      warning(paste("Results file not found for cell type:", ct))
    }
  }
  # init variable
  merged_results <- NULL
  # merge all cts
  if (length(res_per_ct) > 0) {
    merged_results <- rbindlist(res_per_ct, use.names = TRUE, fill = TRUE)
    
  } else {
    warning("No results files found for any cell type.")
  }
  return(merged_results)
}


get_qtls_per_ct <- function(qtl_output_table, category_columns=c('cell_type'), qtl_columns=c('variant', 'tfa', 'gene')) {
  # add columns
  for (category_column in category_columns) {
    # append if column exists
    if ('category' %in% colnames(qtl_output_table)) {
      qtl_output_table[['category']] <- paste(qtl_output_table[['category']], qtl_output_table[[category_column]], sep = ' ')
    }
    # add if does not
    else {
      qtl_output_table[['category']] <- qtl_output_table[[category_column]]
    }
  }
  # now do the same for the QTL columns
  for (qtl_column in qtl_columns) {
    # append if column exists
    if ('qtl' %in% colnames(qtl_output_table)) {
      qtl_output_table[['qtl']] <- paste(qtl_output_table[['qtl']], qtl_output_table[[qtl_column]], sep = ' ')
    }
    # add if does not
    else {
      qtl_output_table[['qtl']] <- qtl_output_table[[qtl_column]]
    }
  }
  # initialize a list of each category
  category_l <- list()
  # check each category
  for (tbl_category in unique(qtl_output_table[['category']])) {
    # get the subset of the table for this category
    category_subset <- qtl_output_table[qtl_output_table[['category']] == tbl_category, ]
    # get the unique QTLs for this category
    unique_qtls <- unique(category_subset[['qtl']])
    # add to the list
    category_l[[tbl_category]] <- unique_qtls
  }
  return(category_l)
}

#' get a label dict that replaces the posix safe names into printable versions
#' 
#' @returns a label dict that replaces the posix safe names into printable versions
#' label_dict_names <- get_label_dict()
get_label_dict <- function() {
  label_dict <- list()
  label_dict[['CD4_T_cells']] <- 'CD4+ T cells'
  label_dict[['CD8_T_cells']] <- 'CD8+ T cells'
  label_dict[['CD4T']] <- 'CD4+ T'
  label_dict[['CD8T']] <- 'CD8+ T'
  label_dict[['CD4_T']] <- 'CD4+ T'
  label_dict[['CD8_T']] <- 'CD8+ T'
  label_dict[['DC']] <- 'DC'
  label_dict[['Dendritic_cells']] <- 'Dendritic cells'
  label_dict[['Endothelial_cells']] <- 'Endothelial cells'
  label_dict[['Fibroblasts']] <- 'Fibroblasts'
  label_dict[['Glia_cells']] <- 'Glia cells'
  label_dict[['Mast_cells']] <- 'MAST cells'
  label_dict[['Mature_absorptive_enterocytes']] <- 'Mature absorptive enterocytes'
  label_dict[['Mature_secretory_enterocytes']] <- 'Mature secretory enterocytes'
  label_dict[['Memory_B']] <- 'Memory B cells'
  label_dict[['Monocytes']] <- 'Monocytes'
  label_dict[['monocyte']] <- 'Monocyte'
  label_dict[['Mono']] <- 'Monocyte'
  label_dict[['NK']] <- 'NK'
  label_dict[['Plasma_cells']] <- 'Plasma cells'
  label_dict[['Stem_cells']] <- 'Stem cells'
  label_dict[['Stromal_cells']] <- 'Stromal cells'
  label_dict[['T_others']] <- 'other T cells'
  label_dict[['Transit_amplifying_cells']] <- 'Transit amplifying cells'
  label_dict[['AI']] <- 'Actively Inflamed'
  label_dict[['NI']] <- 'Non-Inflamed'
  return(label_dict)
}


rename_labels <- function(vector_to_rename) {
  # get the labels that are present
  label_renaming <- get_label_dict()
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


remap_with_label_dict <- function(vector_of_names) {
  # get the label dict
  relabels <- get_label_dict()
  # get the labels available for renaming
  labels_available <- names(relabels)
  # get the ones we cant remap
  unmappable <- setdiff(unique(vector_of_names), labels_available)
  # report on those
  if (length(unmappable) > 0) {
    print(paste('cannot remap the following names, they will be returned unchanged:', paste(unmappable, collapse = ',')))
    # and put those in our remapping list as their originals
    relabels[unmappable] <- unmappable
  }
  # now actually do the remapping
  remapped <- as.vector(unlist(relabels[vector_of_names]))
  return(remapped)
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
  color_coding_dict[['monocyte']] <- '#EDBA1B'
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


split_string_with_max_length <- function(string, max_length=15) {
  # split the string on the space first
  words <- strsplit(string, ' ')[[1]]
  # and initial pasted string
  final_string <- ''
  # and the temporary string
  tmp_string <- ''
  # check each word
  for (word in words) {
    # get the lenght of the temporary string
    tmp_len <- nchar(tmp_string)
    # get the length of the word
    word_len <- nchar(word)
    # length of the final string
    final_len <- nchar(final_string)
    # check if adding the next word would make us go over the limit
    if ((tmp_len + word_len) > max_length) {
      if (final_len > 0) {
        # add the temporary word to the final string if it was not empty
        final_string <- paste(final_string, '\n', tmp_string, sep = '')
      }
      else{
        # if there was not final string yet, it will just become that temporary string
        final_string <- tmp_string
      }
      # reset the temporary string
      tmp_string <- word
    }
    else{
      # otherwise add it to the temporary string
      if (tmp_len > 0) {
        tmp_string <- paste(tmp_string, word)
      }
      # or set if the temporary string was empty
      else{
        tmp_string <- word
      }
    }
  }
  # add last word
  if (nchar(paste(final_string, tmp_string)) > max_length) {
    if (nchar(final_string) == 0) {
      final_string <- tmp_string
    }
    else if (nchar(tmp_string > 0)) {
      final_string <- paste(final_string, '\n', tmp_string, sep = '')
    }
    else{
      # nothing the final string is as it is
    }
  }
  else {
    if (nchar(final_string) == 0) {
      final_string <- tmp_string
    }
    else if(nchar(tmp_string) == 0){
      # nothing, the final string is as is
    }
    else{
      final_string <- paste(final_string, ' ', tmp_string, sep = '')
    }
  }
  # add the last word
  return(final_string)
}


split_strings_with_max_length <- function(strings, max_length=15) {
  strings_split <- sapply(strings, FUN = function(x){
    split_string_with_max_length(x, max_length)
  })
  strings_split <- as.vector(strings_split)
  return(strings_split)
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
    cell_type_colours_missing <- get_color_list(setdiff(cell_types, names(cell_type_colours)))
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
    if (length(queries) == 0) {
      upset(total_overlap, order.by = 'freq', nsets = length(genes_per_ct), sets.bar.color=sets.bar.color, nintersects = n_intersects)
    }
    else {
      upset(total_overlap, order.by = 'freq', nsets = length(genes_per_ct), queries = queries, sets.bar.color=sets.bar.color, nintersects = n_intersects)
    }
  }
  else {
    upset(total_overlap, order.by = 'freq', nsets = length(genes_per_ct), nintersects = n_intersects)
  }
  
  #return(DE_genes_per_ct)
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
    annotate("text", x = max_sig_z_d1 * -0.70 , y = max_sig_z_d2 * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold') +
    # add the names of the concordant and non-concordant blocks
    annotate("text", x = max_sig_z_d1 * 0.70 , y = max_sig_z_d2 * 0.75, label = 'concordant', colour = '#0072B2', fontface = 'bold') +
    # add the title
    ggtitle('Concordance of effects between datasets')
  return(p)
}


####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- F


####################
# Main code        #
####################

# location of output per cell type
tfaieqtl_results_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/'
tfaieqtl_combined_append <- '/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz'
tfaieqtl_ut_append <- 'ut/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz'
tfaieqtl_24hca_append <- '24hca/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz'
# location of the tf-i-ieqtls
tfaieqtl_all_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr_with_replication.tsv.gz'

# get combined results
tfaieqtl_combined <- read_all_tfaieqtl_files(tfaieqtl_results_dir, tfaieqtl_combined_append, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# UT results
tfaieqtl_ut <- read_all_tfaieqtl_files(tfaieqtl_results_dir, tfaieqtl_ut_append, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# 24hCA results
tfaieqtl_24hca <- read_all_tfaieqtl_files(tfaieqtl_results_dir, tfaieqtl_24hca_append, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# add condition to UT and 24hCA
tfaieqtl_ut[['condition']] <- 'UT'
tfaieqtl_24hca[['condition']] <- '24hCA'
# and merge them
tfaieqtl_separate <- rbindlist(list(tfaieqtl_ut, tfaieqtl_24hca), fill = TRUE)

# rename columns
colnames(tfaieqtl_separate) <- gsub('region', 'tfa', colnames(tfaieqtl_separate))
colnames(tfaieqtl_combined) <- gsub('region', 'tfa', colnames(tfaieqtl_combined))
# add explicit tf column
tfaieqtl_separate[['tf']] <- gsub("_(extended|direct).*", "", tfaieqtl_separate[['tfa']])
tfaieqtl_combined[['tf']] <- gsub("_(extended|direct).*", "", tfaieqtl_combined[['tfa']])

# filter
tfaieqtl_separate_sig <- tfaieqtl_separate[
    !is.na(tfaieqtl_separate[['tfa:genotype_bh']]) &
      tfaieqtl_separate[['tfa:genotype_bh']] < 0.05 &
      tfaieqtl_separate[['anova_bh']] < 0.05 &
      tfaieqtl_separate[['genotype_bh']] < 0.05 &
      tfaieqtl_separate[['tfa_bh']] < 0.05,
  ]
tfaieqtl_combined_sig <- tfaieqtl_combined[
  !is.na(tfaieqtl_combined[['tfa:genotype_bh']]) &
    tfaieqtl_combined[['tfa:genotype_bh']] < 0.05 &
    tfaieqtl_combined[['anova_bh']] < 0.05 &
    tfaieqtl_combined[['genotype_bh']] < 0.05 &
    tfaieqtl_combined[['tfa_bh']] < 0.05,
]
# get combined at triplet level
combined_qtls_triplet <- get_qtls_per_ct(tfaieqtl_combined_sig)
# and plot
plot_sharing_per_celltype(combined_qtls_triplet, use_color_dict = T, use_label_dict = T)

# get combined at triplet level
separate_qtls_triplet <- get_qtls_per_ct(tfaieqtl_separate_sig, category_columns=c('cell_type', 'condition'))
# and plot
plot_sharing_per_celltype(separate_qtls_triplet, use_color_dict = T, use_label_dict = T)


# save the result as well
tfaieqtl_combined_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/merged/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr_allct.tsv.gz'
write.table(tfaieqtl_combined, gzfile(tfaieqtl_combined_out_loc), sep = '\t', row.names = F, col.names = T)
mdfiver::create_sha256_for_file(tfaieqtl_combined_out_loc)

# read the old inputs
tfaieqtl_all <- fread(tfaieqtl_all_loc, header = T, sep = '\t')
# make tfa
colnames(tfaieqtl_all) <- gsub('tf', 'tfa', colnames(tfaieqtl_all))
# # get what is significant
tfaieqtl_all_sig <- tfaieqtl_all[
  !is.na(tfaieqtl_all[['tfa:genotype_bh']]) &
    tfaieqtl_all[['tfa:genotype_bh']] < 0.05 &
    tfaieqtl_all[['anova_bh']] < 0.05 &
    # tf_ieqtls[['region_bh']] < 0.05 &
    # tf_ieqtls[['genotype_bh']] < 0.05,
    tfaieqtl_all[['tfa_bh']] < 0.05,
]
# add the old run
combined_qtls_triplet[['all']] <- paste(tfaieqtl_all_sig[['variant']], tfaieqtl_all_sig[['tfa']], tfaieqtl_all_sig[['gene']])
# and plot
plot_sharing_per_celltype(combined_qtls_triplet, use_color_dict = T, use_label_dict = T)

# add the old run
separate_qtls_triplet[['all']] <- paste(tfaieqtl_all_sig[['variant']], tfaieqtl_all_sig[['tfa']], tfaieqtl_all_sig[['gene']])
# and plot
plot_sharing_per_celltype(separate_qtls_triplet, use_color_dict = T, use_label_dict = T)

# look at TFa-gene level
combined_qtls_tfagene <- get_qtls_per_ct(tfaieqtl_combined_sig, qtl_columns = c('tfa', 'gene'))
combined_qtls_tfagene[['all']] <- paste(tfaieqtl_all_sig[['tfa']], tfaieqtl_all_sig[['gene']])
plot_sharing_per_celltype(combined_qtls_tfagene, use_color_dict = T, use_label_dict = T)

# look at TFa-gene level
separate_qtls_tfagene <- get_qtls_per_ct(tfaieqtl_separate_sig, qtl_columns = c('tfa', 'gene'), category_columns = c('cell_type', 'condition'))
separate_qtls_tfagene[['all']] <- paste(tfaieqtl_all_sig[['tfa']], tfaieqtl_all_sig[['gene']])
plot_sharing_per_celltype(separate_qtls_tfagene, use_color_dict = T, use_label_dict = T)


# add the sign
tfaieqtl_combined_sig[['sign']] <- ifelse(tfaieqtl_combined_sig[['tfa:genotype_beta']] > 0, 'up', 'down')
tfaieqtl_separate_sig[['sign']] <- ifelse(tfaieqtl_separate_sig[['tfa:genotype_beta']] > 0, 'up', 'down')
tfaieqtl_all_sig[['sign']] <- ifelse(tfaieqtl_all_sig[['tfa:genotype_beta']] > 0, 'up', 'down')
# get combined at triplet level
combined_qtls_triplet_signed <- get_qtls_per_ct(tfaieqtl_combined_sig, qtl_columns = c('variant', 'tfa', 'gene'),category_columns = c('cell_type', 'sign'))
combined_qtls_triplet_signed[['all_down']] <- paste(tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'down', ][['variant']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'down', ][['tf']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'down', ][['gene']])
combined_qtls_triplet_signed[['all_up']] <- paste(tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'up', ][['variant']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'up', ][['tf']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'up', ][['gene']])
# and plot
plot_sharing_per_celltype(combined_qtls_triplet_signed, use_color_dict = F, use_label_dict = F)
# get separate at triplet level
separate_qtls_triplet_signed <- get_qtls_per_ct(tfaieqtl_separate_sig, qtl_columns = c('variant', 'tfa', 'gene'),category_columns = c('cell_type', 'condition', 'sign'))
separate_qtls_triplet_signed[['all_down']] <- paste(tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'down', ][['variant']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'down', ][['tf']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'down', ][['gene']])
separate_qtls_triplet_signed[['all_up']] <- paste(tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'up', ][['variant']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'up', ][['tf']], tfaieqtl_all_sig[tfaieqtl_all_sig[['sign']] == 'up', ][['gene']])
# and plot
plot_sharing_per_celltype(separate_qtls_triplet_signed, use_color_dict = F, use_label_dict = F)

# merge the old and new results
tfaieqtl_combined_to_all_sig <- merge(tfaieqtl_all_sig, tfaieqtl_combined_sig, by.x = c('variant', 'tfa', 'gene'), by.y = c('variant', 'tfa', 'gene'), suffixes = c('_all', '_ct'))
# show all effects
plot_concondance(tfaieqtl_combined_to_all_sig, d1_effect_column = 'tfa:genotype_beta_all', d2_effect_column = 'tfa:genotype_beta_ct') +
  xlab('all cells TFa-i-eQTL beta') +
  ylab('per celltype TFa-i-eQTL beta')
# show top
tfaieqtl_combined_to_all_sig <- tfaieqtl_combined_to_all_sig[order(abs(tfaieqtl_combined_to_all_sig[['tfa:genotype_beta_ct']]), decreasing = T), ]
tfaieqtl_combined_to_all_sig_top <- tfaieqtl_combined_to_all_sig[!duplicated(paste(tfaieqtl_combined_to_all_sig[['variant']], tfaieqtl_combined_to_all_sig[['tfa']], tfaieqtl_combined_to_all_sig[['gene']])), ]
plot_concondance(tfaieqtl_combined_to_all_sig_top, d1_effect_column = 'tfa:genotype_beta_all', d2_effect_column = 'tfa:genotype_beta_ct') +
  xlab('all cells TFa-i-eQTL beta') +
  ylab('top per celltype TFa-i-eQTL beta')

# repeate for the per-ct
tfaieqtl_separate_to_all_sig <- merge(tfaieqtl_all_sig, tfaieqtl_separate_sig, by.x = c('variant', 'tfa', 'gene'), by.y = c('variant', 'tfa', 'gene'), suffixes = c('_all', '_ct'))
# show all effects
plot_concondance(tfaieqtl_separate_to_all_sig, d1_effect_column = 'tfa:genotype_beta_all', d2_effect_column = 'tfa:genotype_beta_ct') +
  xlab('all cells TFa-i-eQTL beta') +
  ylab('per celltype TFa-i-eQTL beta')
# show top
tfaieqtl_separate_to_all_sig <- tfaieqtl_separate_to_all_sig[order(abs(tfaieqtl_separate_to_all_sig[['tfa:genotype_beta_ct']]), decreasing = T), ]
tfaieqtl_separate_to_all_sig_top <- tfaieqtl_separate_to_all_sig[!duplicated(paste(tfaieqtl_separate_to_all_sig[['variant']], tfaieqtl_separate_to_all_sig[['tfa']], tfaieqtl_separate_to_all_sig[['gene']])), ]
plot_concondance(tfaieqtl_separate_to_all_sig_top, d1_effect_column = 'tfa:genotype_beta_all', d2_effect_column = 'tfa:genotype_beta_ct') +
  xlab('all cells TFa-i-eQTL beta') +
  ylab('top per celltype TFa-i-eQTL beta')
# finally, do UT and 24hCA separately
tfaieqtl_separate_to_all_sig_top_ut <- tfaieqtl_separate_to_all_sig[tfaieqtl_separate_to_all_sig[['condition']] == 'UT', ]
tfaieqtl_separate_to_all_sig_top_ut <- tfaieqtl_separate_to_all_sig_top_ut[!duplicated(paste(tfaieqtl_separate_to_all_sig_top_ut[['variant']], tfaieqtl_separate_to_all_sig_top_ut[['tfa']], tfaieqtl_separate_to_all_sig_top_ut[['gene']])), ]
tfaieqtl_separate_to_all_sig_top_24hca <- tfaieqtl_separate_to_all_sig[tfaieqtl_separate_to_all_sig[['condition']] == '24hCA', ]
tfaieqtl_separate_to_all_sig_top_24hca <- tfaieqtl_separate_to_all_sig_top_24hca[!duplicated(paste(tfaieqtl_separate_to_all_sig_top_24hca[['variant']], tfaieqtl_separate_to_all_sig_top_24hca[['tfa']], tfaieqtl_separate_to_all_sig_top_24hca[['gene']])), ]
# and show them
plot_grid(
  plot_concondance(tfaieqtl_separate_to_all_sig_top_ut, d1_effect_column = 'tfa:genotype_beta_all', d2_effect_column = 'tfa:genotype_beta_ct') +
    xlab('all cells TFa-i-eQTL beta') +
    ylab('top per celltype TFa-i-eQTL beta') + 
    ggtitle('UT'), 
  plot_concondance(tfaieqtl_separate_to_all_sig_top_24hca, d1_effect_column = 'tfa:genotype_beta_all', d2_effect_column = 'tfa:genotype_beta_ct') +
    xlab('all cells TFa-i-eQTL beta') +
    ylab('top per celltype TFa-i-eQTL beta') +
    ggtitle('24hCA')
)
# combined condition and cell type
tfaieqtl_separate_sig[['cell_type_condition']] <- paste(tfaieqtl_separate_sig[['cell_type']], tfaieqtl_separate_sig[['condition']], sep = ' ')
# take the 'all'
tfaieqtl_all_sig[['cell_type_condition']] <- 'all'
tfaieqtl_all_sig[['cell_type']] <- 'all'
tfaieqtl_all_sig[['condition']] <- 'all'
tfaieqtl_all_sig[['tf']] <- gsub("_(extended|direct).*", "", tfaieqtl_all_sig[['tfa']])
# merge that onto the separate sig
tfaieqtl_separate_sig_wall <- rbind(
  tfaieqtl_separate_sig,
  tfaieqtl_all_sig[, colnames(tfaieqtl_separate_sig), with = F]
)
# save plots per combination
ct_comb_plots <- list()
# also keep one with only the ones that have concordances that are not 1
ct_comb_plots_conc <- list()
# check each
for (i in 1 : (length(unique(tfaieqtl_separate_sig_wall[['cell_type_condition']]))-1)) {
  # get the condition and ct
  ct_condition <- unique(tfaieqtl_separate_sig_wall[['cell_type_condition']])[i]
  # get the subset
  ct_condition_subset <- tfaieqtl_separate_sig_wall[tfaieqtl_separate_sig_wall[['cell_type_condition']] == ct_condition, ]
  # check each other condition subset
  for (i2 in (i+1):length(unique(tfaieqtl_separate_sig_wall[['cell_type_condition']]))) {
    # get the condition and ct
    ct_condition2 <- unique(tfaieqtl_separate_sig_wall[['cell_type_condition']])[i2]
    # get the subset
    ct_condition_subset2 <- tfaieqtl_separate_sig_wall[tfaieqtl_separate_sig_wall[['cell_type_condition']] == ct_condition2, ]
    # merge them
    ct_condition_subset_merged <- merge(ct_condition_subset, ct_condition_subset2, by.x = c('variant', 'tfa', 'gene'), by.y = c('variant', 'tfa', 'gene'), suffixes = c('_ct1', '_ct2'))
    # plot if anything overlaps
    if (nrow(ct_condition_subset_merged) > 0) {
      p <- plot_concondance(ct_condition_subset_merged, d1_effect_column = 'tfa:genotype_beta_ct1', d2_effect_column = 'tfa:genotype_beta_ct2') +
        xlab(paste('TFa-i-eQTL beta in', ct_condition)) +
        ylab(paste('TFa-i-eQTL beta in', ct_condition2)) +
        ggtitle(paste(ct_condition, 'and', ct_condition2))
      # add to list
      ct_comb_plots[[paste(ct_condition, ct_condition2, sep = 'to')]] <- p
      # and specifically to the non-concordant list if the concordance is not 1
      n_significant_both <- nrow(ct_condition_subset_merged)
      n_significant_directed_both <- sum(sign(ct_condition_subset_merged[['tfa:genotype_beta_ct1']]) == sign(ct_condition_subset_merged[['tfa:genotype_beta_ct2']]))
      concordance <- round(n_significant_directed_both / n_significant_both, digits = 3)
      if (concordance < 1) {
        ct_comb_plots_conc[[paste(ct_condition, ct_condition2, sep = 'to')]] <- p
      }
    }
  }
}
# show them
p_all <- plot_grid(
  plotlist = ct_comb_plots
)
ggsave('~/multiome/plots/mo_ct_tfaieqtl_overlap_concordances.pdf', p_all, width = 20, height = 20)
# show them
p_all_noconcordant <- plot_grid(
  plotlist = ct_comb_plots_conc
)
ggsave('~/multiome/plots/mo_ct_tfaieqtl_overlap_concordances_not1.pdf', p_all_noconcordant, width = 10, height = 10)

# add condition of "both" to the combined sig
tfaieqtl_combined_sig[['condition']] <- 'both'
# now just compare cell types
tfaieqtl_separate_and_combined_sig <- rbind(
  tfaieqtl_separate_sig,
  tfaieqtl_combined_sig, 
  fill = T
)
# order by significance
tfaieqtl_combined_sig <- tfaieqtl_combined_sig[order(tfaieqtl_combined_sig[['tfa:genotype_bh']]), ]
# save plot per combination
ct_comparisons_per_condition <- list()
# save the concordances
ct_concordances_per_condition <- list()
# check each cell type
for (cell_type_value in unique(tfaieqtl_separate_and_combined_sig[['cell_type']])) {
  # get the subset
  ct_subset <- tfaieqtl_separate_and_combined_sig[tfaieqtl_separate_and_combined_sig[['cell_type']] == cell_type_value, ]
  # do the comparisons per condition
  for (idx1 in 1: (length(unique(ct_subset[['condition']]))-1)) {
    # extract condition
    condition_value <- unique(ct_subset[['condition']])[idx1]
    # get the subset
    condition_subset <- ct_subset[ct_subset[['condition']] == condition_value, ]
    # be sure to have only one row per triplet, as the chunking of MJ does some weird stuff
    condition_subset <- condition_subset[!duplicated(paste(condition_subset[['variant']], condition_subset[['tfa']], condition_subset[['gene']])), ]
    # check each other condition subset
    for (idx2 in (idx1+1):length(unique(ct_subset[['condition']]))) {
      # extract condition
      condition2_value <- unique(ct_subset[['condition']])[idx2]
      # get the subset
      condition_subset2 <- ct_subset[ct_subset[['condition']] == condition2_value, ]
      # be sure to have only one row per triplet, as the chunking of MJ does some weird stuff
      condition_subset2 <- condition_subset2[!duplicated(paste(condition_subset2[['variant']], condition_subset2[['tfa']], condition_subset2[['gene']])), ]
      # merge them
      condition_subset_merged <- merge(condition_subset, condition_subset2, by.x = c('variant', 'tfa', 'gene'), by.y = c('variant', 'tfa', 'gene'), suffixes = c(paste0('_', condition_value), paste0('_', condition2_value)))
      # plot if anything overlaps
      if (nrow(condition_subset_merged) > 0) {
        p <- plot_concondance(condition_subset_merged, d1_effect_column = paste('tfa:genotype_beta', condition_value, sep = '_'), d2_effect_column = paste('tfa:genotype_beta', condition2_value, sep = '_')) +
          xlab(paste('TFa-i-eQTL beta in', cell_type_value, condition_value)) +
          ylab(paste('TFa-i-eQTL beta in', cell_type_value, condition2_value)) +
          ggtitle(paste(cell_type_value, condition_value, 'and', condition2_value))
        # add to list
        ct_comparisons_per_condition[[paste(cell_type_value, condition_value, condition2_value, sep = '_')]] <- p
        # show plot
        print(p)
        # and specifically to the non-concordant list if the concordance is not 1
        n_significant_both <- nrow(condition_subset_merged)
        n_significant_directed_both <- sum(sign(condition_subset_merged[[paste('tfa:genotype_beta', condition_value, sep = '_')]]) == sign(condition_subset_merged[[paste('tfa:genotype_beta', condition2_value, sep = '_')]]))
        concordance <- round(n_significant_directed_both / n_significant_both, digits = 3)
        # put in list
        ct_concordances_per_condition[[paste(cell_type_value, condition_value, condition2_value, sep = '_')]] <- data.frame(
          'cell_type' = cell_type_value,
          'condition1' = condition_value,
          'condition2' = condition2_value,
          'n_significant_both' = n_significant_both,
          'n_significant_directed_both' = n_significant_directed_both,
          'concordance' = concordance
        )
      }
    }
  }
}
# show them
p_ct_condition_comparisons <- plot_grid(
  plotlist = ct_comparisons_per_condition
)
p_ct_condition_comparisons
ggsave('~/multiome/plots/mo_ct_tfaieqtl_ct_condition_concordances.pdf', p_ct_condition_comparisons, width = 20, height = 20)
# show the concordances
ct_concordances_per_condition_all <- do.call(rbind, ct_concordances_per_condition)
