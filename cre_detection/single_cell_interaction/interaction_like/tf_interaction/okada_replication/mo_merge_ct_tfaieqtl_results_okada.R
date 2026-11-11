#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_merge_ct_tfaieqtl_results_okada.R
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
library(stringr)

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
      warning(paste("Results file not found for cell type:", ct, 'at', results_file))
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
    color_coding_dict[[paste(cell_type, 'shared')]] <- color_coding_dict[[cell_type]]
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


plot_concondance_dual <- function(dataset_to_compare, d1_effect_column='d1_zscore', d2_effect_column='d2_zscore', subset_column='genotype_concordant', second_concordance_value='concordant') {
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
  # and the second concordance
  concordant_subset <- dataset_to_compare[dataset_to_compare[[subset_column]] == second_concordance_value, ]
  subset_n_significant_both <- nrow(concordant_subset)
  subset_n_significant_directed_both <- sum(sign(concordant_subset[[d1_effect_column]]) == sign(concordant_subset[[d2_effect_column]]))
  subset_concordance <- round(subset_n_significant_directed_both / subset_n_significant_both, digits = 3)
  
  # plot the concordance of the two
  p <- ggplot(data = dataset_to_compare, mapping = aes(x = !! rlang::sym(d1_effect_column), y = !! rlang::sym(d2_effect_column), shape = !! rlang::sym(subset_column))) + 
    geom_point(size = 1) + 
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
    annotate("label", x = max_sig_z_d1 * 0.75 , y = max_sig_z_d2 * -0.75, label = paste('concordance:\n', concordance, '/', subset_concordance)) +
    # add the names of the concordant and non-concordant blocks
    annotate("text", x = max_sig_z_d1 * -0.70 , y = max_sig_z_d2 * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold') +
    # add the names of the concordant and non-concordant blocks
    annotate("text", x = max_sig_z_d1 * 0.70 , y = max_sig_z_d2 * 0.75, label = 'concordant', colour = '#0072B2', fontface = 'bold') +
    # add the title
    ggtitle('Concordance of effects between datasets')
  return(p)
}


gene_lists_to_overlaps <- function(features_per_celltype_condition, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', '24hCA')) {
  # initialize a dataframe with results
  overlaps_df <- data.frame(matrix(, nrow = length(cell_types)*(length(conditions) + 1), ncol = 3))
  # set the columns names
  colnames(overlaps_df) <- c('cell_type', 'cell_type_condition', 'number')
  # keep an index
  i <- 1
  # check each of the cell types
  for (cell_type in cell_types) {
    # get the features per condition
    features_per_condition <- list()
    # check each of the conditions
    for (condition in conditions) {
      # get those features
      features <- features_per_celltype_condition[[paste(cell_type, condition, sep = ' ')]]
      # put in the list
      features_per_condition[[condition]] <- features
    }
    # get the shared and unique features
    features_shared_unique <- list()
    # now check each condition again
    for (condition in conditions) {
      # get the features for this condition
      features <- features_per_condition[[condition]]
      # get the other condition
      other_condition <- setdiff(conditions, condition)
      # get the features for the other condition
      other_features <- do.call('c', features_per_condition[other_condition])
      # get the unique features for this condition
      unique_features <- setdiff(features, other_features)
      # add to the list
      features_shared_unique[[paste(cell_type, condition, sep = ' ')]] <- unique_features
    }
    # what we don't have yet, is shared
    features_shared_unique[[paste(cell_type, 'shared', sep = ' ')]] <- setdiff(do.call('c', features_per_condition), do.call('c', features_shared_unique))
    # now add the numbers to the dataframe
    for (feature_set in names(features_shared_unique)) {
      overlaps_df[i, 'cell_type'] <- cell_type
      overlaps_df[i, 'cell_type_condition'] <- feature_set
      overlaps_df[i, 'number'] <- length(features_shared_unique[[feature_set]])
      i <- i + 1
    }
  }
  # sometimes some conditions are not present, so then the table will be smaller
  overlaps_df <- overlaps_df[!is.na(overlaps_df$cell_type), ]
  return(overlaps_df)
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
tfaieqtl_results_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/okada/'
tfaieqtl_ut_append <- '/merged/results_fdr.tsv.gz'
# location of the frequence table
tfaieqtl_freq_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/okada/genotype/'
# read freqeuences
tfaieqtl_freq_l <- list() 
for (i in 1:22){
  tfaieqtl_freq_l[[i]] <- fread(paste0(tfaieqtl_freq_loc, '/EAS_imputed_hg38_varFiltered_chr', i, '.afreq'), header = T, sep = '\t', comment.char = '')
}
tfaieqtl_freq <- do.call('rbind', tfaieqtl_freq_l)
# update column names
colnames(tfaieqtl_freq) <- c('CHROM','POS','ID','REF','ALT','ALT_FREQS','OBS_CT')
# UT results
tfaieqtl_ut <- read_all_tfaieqtl_files(tfaieqtl_results_dir, result_append = tfaieqtl_ut_append, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# add condition to UT and 24hCA
tfaieqtl_ut[['condition']] <- 'UT'
# and merge them
tfaieqtl_separate <- rbindlist(list(tfaieqtl_ut), fill = TRUE)

# rename columns
colnames(tfaieqtl_separate) <- gsub('^region', 'tfa', colnames(tfaieqtl_separate))
# remove the i-egene from the eregulon (added in replication step)
tfaieqtl_separate[['tfa']] <- gsub('\\)_.*', ')', tfaieqtl_separate[['tfa']])
# add the signatures that come from the direction of the region to tf and gene
tfaieqtl_separate[['tfa_gene_direction']] <- str_extract(tfaieqtl_separate[['tfa']], '\\+\\/\\+|\\-\\/\\-|\\+\\/\\-|\\-\\/\\+')

# add explicit tf column
tfaieqtl_separate[['tf']] <- gsub("_(extended|direct).*", "", tfaieqtl_separate[['tfa']])

# add variant info
tfaieqtl_separate <- merge(tfaieqtl_separate, tfaieqtl_freq, by.x = 'variant', by.y = 'ID', all.x = T)

# filter
tfaieqtl_separate_sig <- tfaieqtl_separate[
  !is.na(tfaieqtl_separate[['tfa:genotype_bh']]) &
    tfaieqtl_separate[['tfa:genotype_bh']] < 0.05 &
    tfaieqtl_separate[['anova_bh']] < 0.05 &
    tfaieqtl_separate[['genotype_p']] < 0.05 &
    tfaieqtl_separate[['tfa_p']] < 0.05,
]
# tfaieqtl_separate_sig <- tfaieqtl_separate[
#   !is.na(tfaieqtl_separate[['tfa:genotype_p']]) &
#     tfaieqtl_separate[['tfa:genotype_p']] < 0.05 &
#     tfaieqtl_separate[['anova_p']] < 0.05 &
#     tfaieqtl_separate[['genotype_p']] < 0.05 &
#     tfaieqtl_separate[['tfa_p']] < 0.05,
# ]

# save the result as well
tfaieqtl_either_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/okada/results_fdr_allct.tsv.gz'
write.table(tfaieqtl_separate, gzfile(tfaieqtl_either_out_loc), sep = '\t', row.names = F, col.names = T)
mdfiver::create_sha256_for_file(tfaieqtl_either_out_loc)

# read the rep for UT and 24hCA mo as well
mo_tfaieqtl_results_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/'
mo_tfaieqtl_ut_append <- 'UT/lane_donor_countrna_inclcaqtls_varinregion_significant//merged/results_fdr.tsv.gz'
#mo_tfaieqtl_ut_append <- 'ut/lane_donor_countrna_inclcaqtls_varinregion//merged/results_fdr.tsv.gz'
mo_tfaieqtl_24hca_append <- '24hCA/lane_donor_countrna_inclcaqtls_varinregion_significant//merged/results_fdr.tsv.gz'
# location of genotype frequences
mo_tfaieqtl_freq_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/'
# read freqeuences
mo_tfaieqtl_freq_l <- list() 
for (i in 1:22){
  mo_tfaieqtl_freq_l[[i]] <- fread(paste0(mo_tfaieqtl_freq_loc, '/EUR_imputed_hg38_varFiltered_chr', i, '.afreq'), header = T, sep = '\t', comment.char = '')
}
mo_tfaieqtl_freq <- do.call('rbind', mo_tfaieqtl_freq_l)
# update columns
colnames(mo_tfaieqtl_freq) <- c('CHROM','POS','ID','REF','ALT','ALT_FREQS','OBS_CT')
# UT results
mo_tfaieqtl_ut <- read_all_tfaieqtl_files(mo_tfaieqtl_results_dir, result_append = mo_tfaieqtl_ut_append, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# add condition to UT and 24hCA
mo_tfaieqtl_ut[['condition']] <- 'UT'
# 24hCA results
mo_tfaieqtl_24hca <- read_all_tfaieqtl_files(mo_tfaieqtl_results_dir, result_append = mo_tfaieqtl_24hca_append, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# add condition to UT and 24hCA
mo_tfaieqtl_24hca[['condition']] <- '24hCA'
# and merge them
mo_tfaieqtl_separate <- rbindlist(list(mo_tfaieqtl_ut, mo_tfaieqtl_24hca), fill = TRUE)
# make region into tfa
colnames(mo_tfaieqtl_separate) <- gsub('^region', 'tfa', colnames(mo_tfaieqtl_separate))
# remove the i-egene from the eregulon (added in replication step)
mo_tfaieqtl_separate[['tfa']] <- gsub('\\)_.*', ')', mo_tfaieqtl_separate[['tfa']])
# add the signatures that come from the direction of the region to tf and gene
mo_tfaieqtl_separate[['tfa_gene_direction']] <- str_extract(mo_tfaieqtl_separate[['tfa']], '\\+\\/\\+|\\-\\/\\-|\\+\\/\\-|\\-\\/\\+')
# add explicit tf column
mo_tfaieqtl_separate[['tf']] <- gsub("_(extended|direct).*", "", mo_tfaieqtl_separate[['tfa']])
# add variant info
mo_tfaieqtl_separate <- merge(mo_tfaieqtl_separate, mo_tfaieqtl_freq, by.x = 'variant', by.y = 'ID', all.x = T)

# also get significant variant-tfa-gene from the original
mo_tfaieqtl_separate_sig <- mo_tfaieqtl_separate[
  !is.na(mo_tfaieqtl_separate[['tfa:genotype_bh']]) &
    mo_tfaieqtl_separate[['genotype_p']] < 0.05 &
    mo_tfaieqtl_separate[['tfa_p']] < 0.05 &
    mo_tfaieqtl_separate[['tfa:genotype_bh']] < 0.05 &
    mo_tfaieqtl_separate[['anova_bh']] < 0.05,
]
# then filter the replication on what was significant in the original
tfaieqtl_separate_sig_mo <- tfaieqtl_separate[
  paste(tfaieqtl_separate[['variant']], tfaieqtl_separate[['tfa']], tfaieqtl_separate[['gene']]) %in% paste(mo_tfaieqtl_separate_sig[['variant']], mo_tfaieqtl_separate_sig[['tfa']], mo_tfaieqtl_separate_sig[['gene']]), 
]
# make into list
tfaieqtl_separate_sig_mo_l <- list()
# check each cell type
for (cell_type in unique(tfaieqtl_separate_sig_mo[['cell_type']])) {
  # subset to that cell type
  tfaieqtl_separate_sig_mo_celltype <- tfaieqtl_separate_sig_mo[tfaieqtl_separate_sig_mo[['cell_type']] == cell_type, ]
  # extract each column that has a p value
  p_columns <- colnames(tfaieqtl_separate_sig_mo_celltype)[grep('_p$', colnames(tfaieqtl_separate_sig_mo_celltype))]
  # do MTC for each of these columns
  for (p_column in p_columns) {
    # check which rows have a value for this p
    p_valid_i <- !is.na(tfaieqtl_separate_sig_mo_celltype[[p_column]])
    # add MTC column
    bh_column <- gsub('_p$', '_bh', p_column)
    tfaieqtl_separate_sig_mo_celltype[[bh_column]] <- NA
    # extract valid p values, and calculate B&H
    bhs <- p.adjust(tfaieqtl_separate_sig_mo_celltype[[p_column]][p_valid_i], method = 'BH')
    # then place those BHs
    tfaieqtl_separate_sig_mo_celltype[[bh_column]][p_valid_i] <- bhs
    # add bf column as well
    bf_column <- gsub('_p$', '_bf', p_column)
    # adjust the valid p values
    bfs <- p.adjust(tfaieqtl_separate_sig_mo_celltype[[p_column]][p_valid_i], method = 'bonferroni')
    # add the bonferroni column
    tfaieqtl_separate_sig_mo_celltype[[bf_column]] <- NA
    # and place the bonferroni adjusted p values
    tfaieqtl_separate_sig_mo_celltype[[bf_column]][p_valid_i] <- bfs
  }
  # put in the list
  tfaieqtl_separate_sig_mo_l[[cell_type]] <- tfaieqtl_separate_sig_mo_celltype
}
# then merge all
tfaieqtl_separate_sig_mo <- rbindlist(tfaieqtl_separate_sig_mo_l, fill = TRUE)

# merge again, but with the new BH and BF values
#tfaieqtl_separate_both <- merge(tfaieqtl_separate_sig_mo, mo_tfaieqtl_separate_sig, by = c('cell_type', 'tfa', 'gene', 'variant'), suffixes = c('_okada', '_mo'), all = F)
tfaieqtl_separate_both <- merge(tfaieqtl_separate_sig_mo, mo_tfaieqtl_separate_sig, by = c('tfa', 'gene', 'variant'), suffixes = c('_okada', '_mo'), all = F, allow.cartesian=TRUE)
# make tfa-genotype beta directed
tfaieqtl_separate_both[['tfa:genotype_beta_okada']] <- tfaieqtl_separate_both[['tfa:genotype_beta_okada']] * sign(tfaieqtl_separate_both[['genotype_beta_okada']])
tfaieqtl_separate_both[['tfa:genotype_beta_mo']] <- tfaieqtl_separate_both[['tfa:genotype_beta_mo']] * sign(tfaieqtl_separate_both[['genotype_beta_mo']])
# add z
tfaieqtl_separate_both[['tfa:genotype_z_okada']] <- tfaieqtl_separate_both[['tfa:genotype_beta_okada']] / tfaieqtl_separate_both[['tfa:genotype_se_okada']]
tfaieqtl_separate_both[['tfa:genotype_z_mo']] <- tfaieqtl_separate_both[['tfa:genotype_beta_mo']] / tfaieqtl_separate_both[['tfa:genotype_se_mo']]
tfaieqtl_separate_both[['genotype_z_okada']] <- tfaieqtl_separate_both[['genotype_beta_okada']] / tfaieqtl_separate_both[['genotype_se_okada']]
tfaieqtl_separate_both[['genotype_z_mo']] <- tfaieqtl_separate_both[['genotype_beta_mo']] / tfaieqtl_separate_both[['genotype_se_mo']]
# sort on strongest effect
tfaieqtl_separate_both <- tfaieqtl_separate_both[order(abs(tfaieqtl_separate_both[['tfa:genotype_z_okada']])), ]
# get significant effects
tfaieqtl_separate_both_sig <-  tfaieqtl_separate_both[
  tfaieqtl_separate_both$anova_bh_mo < 0.05 & 
    tfaieqtl_separate_both$anova_bh_okada < 0.05 &
    #tfaieqtl_separate_both$genotype_bh_mo < 0.05 & 
    tfaieqtl_separate_both$genotype_p_mo < 0.05 &
    #tfaieqtl_separate_both$genotype_bh_okada < 0.05 & 
    tfaieqtl_separate_both$genotype_p_okada < 0.05 & 
    sign(tfaieqtl_separate_both$genotype_beta_mo) == sign(tfaieqtl_separate_both$genotype_beta_okada) &
    #tfaieqtl_separate_both$tfa_bh_okada < 0.05 & 
    tfaieqtl_separate_both$tfa_p_okada < 0.05 &
    #tfaieqtl_separate_both$tfa_bh_mo < 0.05 & 
    tfaieqtl_separate_both$tfa_p_mo < 0.05 & 
    tfaieqtl_separate_both[['tfa:genotype_bh_mo']] < 0.05 &
    tfaieqtl_separate_both[['tfa:genotype_bh_okada']] < 0.05, 
]
# add concordances
tfaieqtl_separate_both_sig[['interaction_concordant']] <- sign(tfaieqtl_separate_both_sig[['tfa:genotype_beta_okada']]) == sign(tfaieqtl_separate_both_sig[['tfa:genotype_beta_mo']])
tfaieqtl_separate_both_sig[['genotype_concordant']] <- sign(tfaieqtl_separate_both_sig[['genotype_beta_okada']]) == sign(tfaieqtl_separate_both_sig[['genotype_beta_mo']])
tfaieqtl_separate_both_sig[['celltype_concordant']] <- tfaieqtl_separate_both_sig[['cell_type_okada']] == tfaieqtl_separate_both_sig[['cell_type_mo']]
# keep the top
tfaieqtl_separate_both_sig_top <- tfaieqtl_separate_both_sig[
  !duplicated(paste(tfaieqtl_separate_both_sig[['variant']], tfaieqtl_separate_both_sig[['tfa']], tfaieqtl_separate_both_sig[['gene']])), 
]
# then plot top only
p_ut24hac_conc_top <- plot_concondance(
  tfaieqtl_separate_both_sig_top, 
  d1_effect_column='tfa:genotype_beta_mo', d2_effect_column='tfa:genotype_beta_okada') + 
  xlab('TFa-genotype beta scMO') + ylab('TFa-genotype beta Okada')
p_ut24hac_conc_top
# 91.2 % concordant
# get the replication
ut24hac_conc_top_cond <- nrow(tfaieqtl_separate_both_sig_top) / nrow(tfaieqtl_separate_both[!duplicated(paste(tfaieqtl_separate_both[['variant']], tfaieqtl_separate_both[['tfa']], tfaieqtl_separate_both[['gene']])), ])
ut24hac_conc_top_cond
# [1] 0.588

# get nominal effects
tfaieqtl_separate_both_nom <-  tfaieqtl_separate_both[
  tfaieqtl_separate_both$anova_p_mo < 0.05 & 
    tfaieqtl_separate_both$anova_p_okada < 0.05 &
    tfaieqtl_separate_both$genotype_p_mo < 0.05 &
    tfaieqtl_separate_both$genotype_p_okada < 0.05 &
    sign(tfaieqtl_separate_both$genotype_beta_mo) == sign(tfaieqtl_separate_both$genotype_beta_okada) &
    tfaieqtl_separate_both$tfa_p_okada < 0.05 &
    tfaieqtl_separate_both$tfa_p_mo < 0.05 & 
    tfaieqtl_separate_both[['tfa:genotype_p_mo']] < 0.05 &
    tfaieqtl_separate_both[['tfa:genotype_p_okada']] < 0.05, 
]
# keep the top
tfaieqtl_separate_both_nom_top <- tfaieqtl_separate_both_nom[
  !duplicated(paste(tfaieqtl_separate_both_nom[['variant']], tfaieqtl_separate_both_nom[['tfa']], tfaieqtl_separate_both_nom[['gene']])), 
]
# plot concordance
p_ut24hac_conc_nom_top <- plot_concondance(
  tfaieqtl_separate_both_nom_top, 
  d1_effect_column='tfa:genotype_beta_mo', d2_effect_column='tfa:genotype_beta_okada') + 
  xlab('TFa-genotype beta scMO') + ylab('TFa-genotype beta Okada')
p_ut24hac_conc_nom_top
# show the ad for what is concordant
ggplot(data = tfaieqtl_separate_both_sig, mapping = aes(x = gene_ad_a, fill = interaction_concordant)) + 
  geom_density(alpha = .25) + 
  scale_fill_manual(values = c('orange', 'darkgreen')) + 
  xlab('Anderson-Darling gene score') +
  ylab('Density') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
ggplot(data = tfaieqtl_separate_both_sig, mapping = aes(x = tfa_ad_a, fill = interaction_concordant)) + 
  geom_density(alpha = .25) + 
  scale_fill_manual(values = c('orange', 'darkgreen')) + 
  xlab('Anderson-Darling TFa score') +
  ylab('Density') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

