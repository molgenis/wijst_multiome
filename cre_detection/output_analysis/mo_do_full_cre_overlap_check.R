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
library(qvalue)


####################
# Functions        #
####################

#' Calculate Nominal Thresholds
#'
#' This function calculates nominal thresholds for p-values based on a given false discovery rate (FDR).
#'
#' @param res_df A data frame containing the results with p-values and other relevant columns.
#' @param fdr A numeric value specifying the false discovery rate threshold. Default is 0.05.
#' @param pval_col A character string specifying the name of the column with p-values. Default is 'p_value'.
#' @param nominal_threshold_column A character string specifying the name of the column to store the nominal thresholds. Default is 'pval_nominal_threshold'.
#' @param cutoff_column A character string specifying the name of the column with feature q-values. Default is 'feature_q_value'.
#' @param alpha_column A character string specifying the name of the column with alpha parameters for the beta distribution. Default is 'alpha_param'.
#' @param beta_column A character string specifying the name of the column with beta parameters for the beta distribution. Default is 'beta_param'.
#'
#' @return A data frame with an additional column for nominal thresholds.
#' @export
#'
#' @examples
#' \dontrun{
#'   res_df <- data.frame(
#'     p_value = runif(100),
#'     feature_q_value = runif(100),
#'     alpha_param = rep(1, 100),
#'     beta_param = rep(1, 100)
#'   )
#'   calculate_nominal_thresholds(res_df)
#' }
calculate_nominal_thresholds <- function(res_df, fdr=0.05, pval_col='p_value', nominal_threshold_column='pval_nominal_threshold', cutoff_column='feature_q_value', alpha_column='alpha_param', beta_column='beta_param') {
  # get the lowerbound p values, so the ones that are smaller than the FDR
  indices_lb <- res_df[[cutoff_column]] < fdr
  lb <- as.vector(res_df[indices_lb, ][[pval_col]])
  # put then in ascending order
  lb <- lb[order(lb)]
  # get the upperbound p values, so the ones that are bigger than the FDR
  indices_ub <- res_df[[cutoff_column]] > fdr
  ub <- as.vector(res_df[indices_ub, ][[pval_col]])
  # and order them
  ub <- ub[order(ub)]
  
  # if we have any significant effects, we can get a cutoff
  if (length(lb) > 0) {
    # get the highest (p) significant value
    highest_in_lb <- tail(lb, 1)
    # if there are any non significant effects
    if (length(ub) > 0) {
      # get the lowest (p) non-significant value
      lowest_in_ub <- head(ub, 1)
      # and calculate the threshold
      pthreshold <- (highest_in_lb + lowest_in_ub) / 2
    } else {
      # otherwise the highest effect will just be the cutoff
      pthreshold <- highest_in_lb
    }
    # ge the threshold, based on the shapes of the beta distribution and the significance threshold
    res_df[[nominal_threshold_column]] <- stats::qbeta(pthreshold, as.vector(res_df[[alpha_column]]), as.vector(res_df[[beta_column]]))
  }
  else {
    # otherwise it would have to be zero
    res_df[[nominal_threshold_column]] <- 0
  }
  return(res_df)
}



read_pseudobulk_cre_output_per_celltype <- function(pseudobulk_output_folder, cell_types=NULL, filename_output='qtl_results_all.txt.gz', significance_column='empirical_feature_p_value', significance_cutoff=0.05, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', add_global_nominal_threshold=F, add_local_nominal_threshold=F, global_nominal_threshold_column_to_add='pval_nominal_threshold_global', local_nominal_threshold_column_to_add='pval_nominal_threshold_local', alpha_column='alpha_param', beta_column='beta_param', nominal_p_column='p_value', filter_alpha=T, alpha_min=.2, alpha_max=5, filter_significance=T) {
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
    # make sure there are no duplicates
    cell_type_output <- unique(cell_type_output)
    # filter on alpha if requested
    if (filter_alpha) {
      cell_type_output <- cell_type_output[!(cell_type_output[[alpha_column]] > alpha_max | cell_type_output[[alpha_column]] < alpha_min), ]
    }
    
    # get the features and the emperical p value
    if (add_mtc) {
      # subset to what we need
      cell_type_output_features <- NULL
      # which is a bit if we care about the nominal threshold
      if (add_global_nominal_threshold) {
        cell_type_output_features <- cell_type_output[, c(..feature_mtc_column, ..mtc_column, ..nominal_p_column, ..alpha_column, ..beta_column), with = F]
      }
      # even less if we don't try to get the nominal threshold as well
      else {
        cell_type_output_features <- cell_type_output[, c(..feature_mtc_column, ..mtc_column), with = F]
      }
      # remove the wherever we dont have our significance
      cell_type_output_features <- cell_type_output_features[!is.na(cell_type_output_features[[significance_column]]) & cell_type_output_features[[significance_column]] >= 0, ]
      # order by significance
      cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
      # keep only the first entry
      cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
      # set the values that are larger than 1, to be 1, problem with precision
      cell_type_output_features[cell_type_output_features[[mtc_column]] > 1, mtc_column] <- 1
      # add multiple testing correction
      cell_type_output_features[['qvalue']] <- qvalue(cell_type_output_features[[mtc_column]])$qvalues
      # now add back to the original table
      cell_type_output[[mtc_column_to_add]] <- cell_type_output_features[match(cell_type_output[[feature_mtc_column]], cell_type_output_features[[feature_mtc_column]]), 'qvalue'][['qvalue']]
      # based on this MTC column, we can now also add a cuttoff
      if (add_local_nominal_threshold) {
        cell_type_output_local_threshold <- calculate_nominal_thresholds(cell_type_output_features, fdr=significance_cutoff, pval_col=nominal_p_column, nominal_threshold_column='nomthres', cutoff_column = 'qvalue', alpha_column = alpha_column, beta_column = beta_column)
        # now add the nominal threshold to the full table
        cell_type_output[[local_nominal_threshold_column_to_add]] <- cell_type_output_local_threshold[match(cell_type_output[[feature_mtc_column]], cell_type_output_local_threshold[[feature_mtc_column]]), 'nomthres'][['nomthres']]
      }
      if(add_global_nominal_threshold) {
        # filter the output to significant MTC hits
        cell_type_output_features_significant <- cell_type_output_features[cell_type_output_features[['qvalue']] < significance_cutoff, ]
        # and get the maximum significant nominal value
        global_p_cutoff <- max(cell_type_output_features_significant[[nominal_p_column]])
        # add that to the table
        cell_type_output[[global_nominal_threshold_column_to_add]] <- global_p_cutoff
      }
    }
    # filter the file if requested
    if (filter_significance) {
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
  output_per_celltype <- read_pseudobulk_cre_output_per_celltype(binomial_output_folder, cell_types = cell_types, filename_output = filename_output, significance_column = significance_column, significance_cutoff = significance_cutoff, add_mtc = F, add_global_nominal_threshold = F, add_local_nominal_threshold = F, filter_alpha = F)
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


get_closest_flanks <- function(position_table, left_flank_column1, right_flank_column1, left_flank_column2, right_flank_column2) {
  # get the distance between left flanks
  dist_left_flank1_to_left_flank2 <- position_table[[left_flank_column1]] - position_table[[left_flank_column2]]
  # distance between the right flanks
  dist_right_flank1_to_right_flank2 <- position_table[[right_flank_column1]] - position_table[[right_flank_column2]]
  # distance between left flank 1 and right flank 2
  dist_left_flank1_to_right_flank2 <- position_table[[left_flank_column1]] - position_table[[right_flank_column2]]
  # distance between right flank 1 and left flank 2
  dist_right_flank1_to_left_flank2 <- position_table[[right_flank_column1]] - position_table[[left_flank_column2]]
  # put in a table for convenience sake
  distances_tbl <- data.table(
    'lf1_to_lf2' = dist_left_flank1_to_left_flank2, 
    'rf1_to_rf2' = dist_right_flank1_to_right_flank2, 
    'lf1_to_rf2' = dist_left_flank1_to_right_flank2, 
    'rf1_to_lf2' = dist_right_flank1_to_left_flank2
  )
  # add the minimum absolute distance
  distances_tbl[['min_dist']] <- apply(distances_tbl, 1, function(x) {
    return(min(abs(x)))
  })
  # but set this to zero if any of the flanks end in the bodies
  #              -----
  #                 ++++
  distances_tbl[(distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['lf1_to_rf2']] > 0) |
                  #                   ----
                #                 ++++
                (distances_tbl[['rf1_to_lf2']] > 0 & distances_tbl[['lf1_to_rf2']] < 0) |
                  #                   ----
                #                 +++++++++
                (distances_tbl[['lf1_to_lf2']] > 0 & distances_tbl[['rf1_to_rf2']] < 0) |
                  #                 ---------
                #                   ++++
                (distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['rf1_to_rf2']] > 0)
                , 'min_dist'] <- 0
  return(distances_tbl)
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
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# location of the overlapping caQTLs and eQTLs
qtl_overlap_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl_caqtl_overlap/combined/L1/all/eqtl_caqtl_overlapping_variants.tsv.gz'
# location of the pseudobulk CRE mapping
pseudobulk_output_folder <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eQTA/eQTA_v2/L1/'
# location of the binomial method
binomial_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/combined/betas_ps/'
# location of the hybrid method
hybrid_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/'

# read the pseudobulk outputs
pseudobulk_output_ut_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, 'UT', sep = '/'), add_mtc = T, filter_alpha = T, add_global_nominal_threshold = T, add_local_nominal_threshold = T)
pseudobulk_output_24hca_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, '24hCA', sep = '/'), add_mtc = T, filter_alpha = T, add_global_nominal_threshold = T, add_local_nominal_threshold = T)
# merge the cell types
pseudobulk_output_ut <- do.call('rbind', pseudobulk_output_ut_list)
pseudobulk_output_24hca <- do.call('rbind', pseudobulk_output_24hca_list)
# add the condition
pseudobulk_output_ut[['condition']] <- 'UT'
pseudobulk_output_24hca[['condition']] <- '24hCA'
# add z score
pseudobulk_output_ut[['zscore']] <- pseudobulk_output_ut[['beta']] / pseudobulk_output_ut[['beta_se']]
pseudobulk_output_24hca[['zscore']] <- pseudobulk_output_24hca[['beta']] / pseudobulk_output_24hca[['beta_se']]
# add a correlation based on the Z score, taking sample size and removing 2 + 10 PCs to get the degrees of freedom
pseudobulk_output_ut[['r']] <- pseudobulk_output_ut[['zscore']] / sqrt(pseudobulk_output_ut[['zscore']]^2 + (pseudobulk_output_ut[['n_samples']][1] - 12))
pseudobulk_output_24hca[['r']] <- pseudobulk_output_24hca[['zscore']] / sqrt(pseudobulk_output_24hca[['zscore']]^2 + (pseudobulk_output_24hca[['n_samples']][1] - 12))
# add a p based z
pseudobulk_output_ut[['z_from_p']] <- qnorm(1 - pseudobulk_output_ut[['p_value']] / 2) * sign(pseudobulk_output_ut[['beta']])
pseudobulk_output_24hca[['z_from_p']] <- qnorm(1 - pseudobulk_output_24hca[['p_value']] / 2) * sign(pseudobulk_output_24hca[['beta']])
pseudobulk_output_ut[['z_from_p']] <- qnorm(pseudobulk_output_ut[['p_value']] / 2) * -1 * sign(pseudobulk_output_ut[['beta']])
pseudobulk_output_24hca[['z_from_p']] <- qnorm(pseudobulk_output_24hca[['p_value']] / 2) * -1 * sign(pseudobulk_output_24hca[['beta']])
# and a clipped z from p
# pseudobulk_output_ut_p_for_z <- pseudobulk_output_ut[['p_value']] / 2
# pseudobulk_output_ut_p_for_z[pseudobulk_output_ut_p_for_z < 1e-16] <- 1e-16
# pseudobulk_output_ut[['z_from_p_clipped']] <- qnorm(1 - pseudobulk_output_ut_p_for_z) * sign(pseudobulk_output_ut[['beta']])
# pseudobulk_output_24hca_p_for_z <- pseudobulk_output_24hca[['p_value']] / 2
# pseudobulk_output_24hca_p_for_z[pseudobulk_output_24hca_p_for_z < 1e-16] <- 1e-16
# pseudobulk_output_24hca[['z_from_p_clipped']] <- qnorm(1 - pseudobulk_output_24hca_p_for_z) * sign(pseudobulk_output_24hca[['beta']])

# merge them
pseudobulk_output <- do.call('rbind', list(pseudobulk_output_ut, pseudobulk_output_24hca))

# read hybrid method
hybrid_output_list <- read_pseudobulk_cre_output_per_celltype(hybrid_output_loc, add_mtc = T, filter_alpha = T, add_global_nominal_threshold = T, add_local_nominal_threshold = T, filename_output = 'qtl_results_all.txt', alpha_min = .8, alpha_max = 1.2)
#hybrid_output_list <- read_pseudobulk_cre_output_per_celltype(hybrid_output_loc, add_mtc = T, filter_alpha = T, add_global_nominal_threshold = T, add_local_nominal_threshold = T, filename_output = 'test_qtl_results_all.txt', alpha_min = .8, alpha_max = 1.2, filter_significance = F)
# merge them
hybrid_output <- do.call('rbind', hybrid_output_list)
# add z score
hybrid_output[['zscore']] <- hybrid_output[['beta']] / hybrid_output[['beta_se']]
# add a correlation based on the Z score, taking sample size and removing 2 + 10 PCs to get the degrees of freedom
hybrid_output[['r']] <- hybrid_output[['zscore']] / sqrt(hybrid_output[['zscore']]^2 + (hybrid_output[['n_samples']][1] - 12))
# add a p based z
hybrid_output[['z_from_p']] <- qnorm(hybrid_output[['p_value']] / 2) * -1 * sign(hybrid_output[['beta']])


# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# read the qtl overlap
qtl_overlap <- fread(qtl_overlap_loc, header = T, sep = '\t')

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/cPeaks/cPeaks_info.tsv'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = ' ')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])

# read the location of the genes
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno <- fread(gene_anno_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')

# add the location to the pseudobulk info
pseudobulk_output <- cbind(pseudobulk_output, cpeaks_anno[match(pseudobulk_output[['snp_id']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# get the distances
pseudobulk_distances <- get_closest_flanks(pseudobulk_output, 'start_hg38', 'end_hg38', 'feature_start', 'feature_end')
# add those distances
pseudobulk_output[['distance']] <- pseudobulk_distances[['min_dist']]
# and category
pseudobulk_output[['category']] <- 'pseudobulk'

# and to the hybrid method
hybrid_output <- cbind(hybrid_output, cpeaks_anno[match(hybrid_output[['snp_id']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
hybrid_distances <- get_closest_flanks(hybrid_output, 'start_hg38', 'end_hg38', 'feature_start', 'feature_end')
hybrid_output[['distance']] <- hybrid_distances[['min_dist']]
hybrid_output[['category']] <- 'hybrid'


# get extra annotations for the pseudobulk output
strand_information_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eQTA/LimixExpAnnotationFile.incStrand.txt'
strand_information <- fread(strand_information_loc, header = T, sep = '\t')
# add to the pseudobulk
pseudobulk_output[['strand']] <- strand_information[match(pseudobulk_output[['feature_id']], strand_information[['feature_id']]), ][['strand']]

# add strand info to the hybrid output as well
hybrid_output[['strand']] <- strand_information[match(hybrid_output[['feature_id']], strand_information[['feature_id']]), ][['strand']]


# read the binomial results
binomial_output_list <- read_binomial_output_per_celltype(binomial_output_loc, cell_types = c('CD4T', 'CD8T', 'NK', 'monocyte'))
# merge them
binomial_output <- do.call('rbind', binomial_output_list)
# add a correlation based on the Z score, taking sample size and removing 2 + 10 PCs to get the degrees of freedom
binomial_output[['r']] <- binomial_output[['meta_z']] / sqrt(binomial_output[['meta_z']]^2 + (binomial_output[['n_sample']][1] - 2))

# and the qtl overlap
qtl_overlap[['z_caqtl']] <- qtl_overlap[['beta_caqtl']] / qtl_overlap[['se_caqtl']]
qtl_overlap[['z_eqtl']] <- qtl_overlap[['beta_eqtl']] / qtl_overlap[['se_eqtl']]
# get the sign overlap
qtl_overlap[['sign']] <- sign(qtl_overlap[['z_caqtl']]) * sign(qtl_overlap[['z_eqtl']])

# add the location of the caQTL here as well
qtl_overlap <- cbind(qtl_overlap, cpeaks_anno[match(qtl_overlap[['feature_caqtl']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# and the locations of the genes
qtl_overlap <- cbind(qtl_overlap, gene_anno[match(qtl_overlap[['feature_eqtl']], gene_anno[['gs']]), c('chrom', 'start', 'end')])
# get the distances
qtl_distances <- get_closest_flanks(qtl_overlap, 'start_hg38', 'end_hg38', 'start', 'end')
# add that to the original table
qtl_overlap[['distance']] <- qtl_distances[['min_dist']]
# and category
qtl_overlap[['category']] <- 'qtl_overlap'

# add location for the binomial table
binomial_output <- cbind(binomial_output, cpeaks_anno[match(binomial_output[['region']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# and the locations of the genes
binomial_output <- cbind(binomial_output, gene_anno[match(binomial_output[['gene']], gene_anno[['gs']]), c('chrom', 'start', 'end')])
# get the distances
binomial_distances <- get_closest_flanks(binomial_output, 'start_hg38', 'end_hg38', 'start', 'end')
# add that to the original table
binomial_output[['distance']] <- binomial_distances[['min_dist']]
# and category
binomial_output[['category']] <- 'binomial'

# add location for the binomial table
scenic_output <- cbind(scenic_output, cpeaks_anno[match(scenic_output[['Region']], cpeaks_anno[['scenic_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# and the locations of the genes
scenic_output <- cbind(scenic_output, gene_anno[match(scenic_output[['Gene']], gene_anno[['gs']]), c('chrom', 'start', 'end')])
# get the distances again
scenic_distances <- get_closest_flanks(scenic_output, 'start_hg38', 'end_hg38', 'start', 'end')
# add that to the original table
scenic_output[['distance']] <- scenic_distances[['min_dist']]
# and category
scenic_output[['category']] <- 'scenic'


# add region to gene column
pseudobulk_output[['r2g']] <- paste(pseudobulk_output[['snp_id']], pseudobulk_output[['feature_id']])
binomial_output[['r2g']] <- paste(binomial_output[['region']], binomial_output[['gene']])
scenic_output[['r2g']] <- paste(gsub(':', '-', scenic_output[['Region']]), scenic_output[['Gene']])
qtl_overlap[['r2g']] <- paste(qtl_overlap[['feature_caqtl']], qtl_overlap[['feature_eqtl']])
hybrid_output[['r2g']] <- paste(hybrid_output[['snp_id']], hybrid_output[['feature_id']])

# sort all of them by the Z
pseudobulk_output <- pseudobulk_output[order(abs(pseudobulk_output[['zscore']]), decreasing = T), ]
binomial_output <- binomial_output[order(abs(binomial_output[['meta_z']]), decreasing = T), ]
scenic_output <- scenic_output[order(abs(scenic_output[['rho_R2G']]), decreasing = T), ]
qtl_overlap <- qtl_overlap[order(abs(qtl_overlap[['z_eqtl']]), abs(qtl_overlap[['z_caqtl']]), decreasing = T), ]
hybrid_output <- hybrid_output[order(abs(hybrid_output[['zscore']]), decreasing = T), ]

# filter on theshold
pseudobulk_output_unfiltered <- pseudobulk_output
# check significance threshold
pseudobulk_output <- pseudobulk_output[
  !is.na(pseudobulk_output[['p_value']]) & !is.na(pseudobulk_output[['pval_nominal_threshold_global']]) & 
    pseudobulk_output[['p_value']] < pseudobulk_output[['pval_nominal_threshold_global']] &
    !is.na(pseudobulk_output[['feature_q_value']]) & pseudobulk_output[['feature_q_value']] < 0.05 , ]

# export
pseudo_ext <- pseudobulk_output[pseudobulk_output[['p_value']] >= 0 & pseudobulk_output[['empirical_feature_p_value']] < 0.05, ]
pseudo_ext <- pseudo_ext[, c('snp_id', 'feature_id', 'p_value', 'zscore', 'z_from_p', 'condition', 'cell_type', 'chr_hg38', 'start_hg38', 'end_hg38', 'feature_chromosome', 'feature_start', 'feature_end', 'distance', 'strand')]
colnames(pseudo_ext) <- c('region', 'gene', 'p_value', 'zscore', 'z_from_p', 'condition', 'cell_type', 'chr_region', 'start_region', 'end_region', 'chr_gene', 'gene_start', 'gene_end', 'distance', 'strand')
pseudo_ext[['chr_gene']] <- paste0('chr', pseudo_ext[['chr_gene']])
pseudo_ext <- pseudo_ext[order(abs(pseudo_ext[['z_from_p']]), decreasing = T), ]
write.table(pseudo_ext, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eQTA/export/mo_pseudobulk_export.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
# check minimal correlation
#pseudobulk_output <- pseudobulk_output[abs(pseudobulk_output[['r']]) >= .25 , ]
# check that the region and gene do not overlap
#pseudobulk_output <- pseudobulk_output[pseudobulk_output[['distance']] > 0 , ]

# filter the hybrid one in the same way
hybrid_output_unfiltered <- hybrid_output
# check significance threshold
hybrid_output <- hybrid_output[
  !is.na(hybrid_output[['p_value']]) & !is.na(hybrid_output[['pval_nominal_threshold_global']]) & 
    hybrid_output[['p_value']] < hybrid_output[['pval_nominal_threshold_global']] &
    !is.na(hybrid_output[['feature_q_value']]) & hybrid_output[['feature_q_value']] < 0.05 , ]

# remove the entries that are more likely to be false positives
scenic_output_unfiltered <- scenic_output
scenic_output <- scenic_output_unfiltered[scenic_output_unfiltered[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# and region-gene overlaps
scenic_output <- scenic_output[scenic_output[['distance']] > 0, ]

# filter also on correlation for the binomial output
binomial_output_unfiltered <- binomial_output
#binomial_output <- binomial_output[abs(binomial_output[['r']]) >= .25, ]
# and remove region-gene overlaps
binomial_output <- binomial_output[abs(binomial_output[['distance']]) > 0, ]

# get unique ones
pseudobulk_output_unique <- pseudobulk_output[!duplicated(paste(pseudobulk_output[['snp_id']], pseudobulk_output[['feature_id']])), ]
binomial_output_unique <- binomial_output[!duplicated(paste(binomial_output[['region']], binomial_output[['gene']])), ]
scenic_output_unique <- scenic_output[!duplicated(paste(scenic_output[['Region']], scenic_output[['Gene']])), ]
scenic_output_unique_unfiltered <- scenic_output_unfiltered[!duplicated(paste(scenic_output_unfiltered[['Region']], scenic_output_unfiltered[['Gene']])), ]
qtl_overlap_unique <- qtl_overlap[!duplicated(paste(qtl_overlap[['feature_caqtl']], qtl_overlap[['feature_eqtl']])), ]
hybrid_output_unique <- hybrid_output[!duplicated(paste(hybrid_output[['snp_id']], hybrid_output[['feature_id']])), , ]

# show how many we have in each set
nrow(binomial_output_unique)
# [1] 12052
nrow(scenic_output_unique_unfiltered)
# [1] 80440
nrow(scenic_output_unique)
# [1] 53796
nrow(pseudobulk_output_unique)
# [1] 1458
nrow(qtl_overlap_unique)
# [1] 7677
nrow(hybrid_output_unique)
# [1] 7515

# make the unique region-gene numbers into a table
n_effects_region_gene <- data.frame(
  'category' = c('binomial', 'SCENIC+ filtered', 'SCENIC+ unfiltered', 'pseudobulk', 'QTL overlap', 'hybrid (chr17 only)'), 
  'neffects' = c(nrow(binomial_output_unique), nrow(scenic_output_unique), nrow(scenic_output_unique_unfiltered), nrow(pseudobulk_output_unique), nrow(qtl_overlap_unique), nrow(hybrid_output))
)
# and make into a plot
ggplot(data = n_effects_region_gene, mapping = aes(x = category, y = neffects, fill = category)) + 
  geom_bar(stat = 'identity') +
  xlab('CRE detection method') +
  ylab('Number of region-gene pairs') +
  ggtitle('Number of detected CRE-gene pairs across methods') +
  scale_fill_manual(values = list('binomial' = '#BEAED4', 'pseudobulk' = '#7FC97F', 'QTL overlap' = '#386CB0', 'SCENIC+ unfiltered' = '#FFFF99', 'SCENIC+ filtered' = '#FDC086', 'hybrid (chr17 only)' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  theme(legend.position="none")


# add region to gene column
pseudobulk_output_unique[['r2g']] <- paste(pseudobulk_output_unique[['snp_id']], pseudobulk_output_unique[['feature_id']])
binomial_output_unique[['r2g']] <- paste(binomial_output_unique[['region']], binomial_output_unique[['gene']])
scenic_output_unique[['r2g']] <- paste(gsub(':', '-', scenic_output_unique[['Region']]), scenic_output_unique[['Gene']])
scenic_output_unique_unfiltered[['r2g']] <- paste(gsub(':', '-', scenic_output_unique_unfiltered[['Region']]), scenic_output_unique_unfiltered[['Gene']])
qtl_overlap_unique[['r2g']] <- paste(qtl_overlap_unique[['feature_caqtl']], qtl_overlap_unique[['feature_eqtl']])
hybrid_output_unique[['r2g']] <- paste(hybrid_output_unique[['snp_id']], hybrid_output_unique[['feature_id']])


# plot these numbers
plot_sharing_per_celltype(
  list('pseudobulk' = pseudobulk_output_unique[['r2g']], 
       'binomial' = binomial_output_unique[['r2g']], 
       'scenic' = scenic_output_unique[['r2g']], 
       'scenic unfiltered' = scenic_output_unique_unfiltered[['r2g']],
       'qtl' = qtl_overlap_unique[['r2g']], 
       'hybrid' = hybrid_output_unique[['r2g']]), 
  use_label_dict=F, use_color_dict=T
)


# get the percentage positive
frac_pos_pseudobulk_output_unique <- nrow(pseudobulk_output_unique[sign(pseudobulk_output_unique[['zscore']]) == 1, ]) / nrow(pseudobulk_output_unique)
# [1] 0.8462757
frac_pos_binomial_output_unique <- nrow(binomial_output_unique[sign(binomial_output_unique[['meta_z']]) == 1, ]) / nrow(binomial_output_unique)
# [1] 0.9983478
frac_pos_scenic_output_unique <- nrow(scenic_output_unique[sign(scenic_output_unique[['rho_R2G']]) == 1, ]) / nrow(scenic_output_unique)
# [1] 1
frac_pos_scenic_output_unfiltered_unique <- nrow(scenic_output_unique_unfiltered[sign(scenic_output_unique_unfiltered[['rho_R2G']]) == 1, ]) / nrow(scenic_output_unique_unfiltered)
# [1] 0.6687718
frac_pos_qtl_overlap_unique <- nrow(qtl_overlap_unique[sign(qtl_overlap_unique[['sign']]) == 1, ]) / nrow(qtl_overlap_unique)
# [1] 0.767357
frac_pos_hybrid_output_unique <- nrow(hybrid_output_unique[sign(hybrid_output_unique[['zscore']]) == 1, ]) / nrow(hybrid_output_unique)
# [1] 0.9753086


# plot these numbers
n_pos_tbl <- data.frame(
  'method' = c('pseudobulk', 'binomial', 'scenic', 'scenic uf', 'QTL', 'hybrid (chr17 only)', 'pseudobulk', 'binomial', 'scenic', 'scenic uf', 'QTL', 'hybrid (chr17 only)'), 
  'direction' = c('positive', 'positive', 'positive', 'positive', 'positive', 'positive', 'negative', 'negative', 'negative', 'negative', 'negative', 'negative'), 
  'n' = c(frac_pos_pseudobulk_output_unique, frac_pos_binomial_output_unique, frac_pos_scenic_output_unique, frac_pos_scenic_output_unfiltered_unique, frac_pos_qtl_overlap_unique, frac_pos_hybrid_output_unique, 1-frac_pos_pseudobulk_output_unique, 1-frac_pos_binomial_output_unique, 1-frac_pos_scenic_output_unique, 1-frac_pos_scenic_output_unfiltered_unique, 1-frac_pos_qtl_overlap_unique, 1-frac_pos_hybrid_output_unique)
)
p_directions <- ggplot(data = n_pos_tbl, mapping = aes(x = method, y = n, fill = direction)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('CRE detection method') +
  ylab('fraction') +
  ggtitle('CRE method result directions') +
  scale_fill_manual(values = list('positive' = 'darkblue', 'negative' = 'darkred')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_directions

# show the distances in a density plot
p_region_direction_distances <- ggplot(
  data = rbind(
    pseudobulk_output_unique[pseudobulk_output_unique[['distance']] < 150000 & pseudobulk_output_unique[['distance']] > 0, c('distance', 'category'), ], 
    binomial_output_unique[binomial_output_unique[['distance']] < 150000 & binomial_output_unique[['distance']] > 0, c('distance', 'category'), ], 
    qtl_overlap_unique[qtl_overlap_unique[['distance']] < 150000 & qtl_overlap_unique[['distance']] > 0, c('distance', 'category'), ], 
    hybrid_output_unique[hybrid_output_unique[['distance']] < 150000 & hybrid_output_unique[['distance']] > 0, c('distance', 'category'), ]), 
  mapping = aes(
    x = distance, 
    fill = category
  )
) + 
  geom_density(alpha = 0.5) +
  xlab('Distance between region and gene') + 
  ylab('Density') + 
  ggtitle('Distance between region and gene\nacross different methods') + 
  scale_fill_manual(values = list('binomial' = '#BEAED4', 'pseudobulk' = '#7FC97F', 'qtl_overlap' = '#386CB0', 'hybrid' = '#F0027F' )) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# show plot
p_region_direction_distances
# do same for positive direction
p_region_direction_distances_pos <- ggplot(
  data = rbind(
    pseudobulk_output_unique[pseudobulk_output_unique[['distance']] > 0 & pseudobulk_output_unique[['r']] > 0, c('distance', 'category'), ], 
    binomial_output_unique[binomial_output_unique[['distance']] > 0 & binomial_output_unique[['r']] > 0, c('distance', 'category'), ], 
    qtl_overlap_unique[qtl_overlap_unique[['distance']] > 0 & qtl_overlap_unique[['sign']] > 0, c('distance', 'category'), ], 
    hybrid_output_unique[hybrid_output_unique[['distance']] > 0 & hybrid_output_unique[['distance']] > 0, c('distance', 'category'), ]), 
  mapping = aes(
    x = distance, 
    fill = category
  )
) + 
  geom_density(alpha = 0.5) +
  xlab('Distance between region and gene') + 
  ylab('Density') + 
  ggtitle('Distance between region and gene\nacross different methods\nfor positive associations') + 
  scale_fill_manual(values = list('binomial' = '#BEAED4', 'pseudobulk' = '#7FC97F', 'qtl_overlap' = '#386CB0', 'hybrid' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_region_direction_distances_pos
# finally for negative direction
p_region_direction_distances_neg <- ggplot(
  data = rbind(
    pseudobulk_output_unique[pseudobulk_output_unique[['distance']] < 150000 & pseudobulk_output_unique[['distance']] > 0 & pseudobulk_output_unique[['r']] < 0, c('distance', 'category'), ], 
    binomial_output_unique[binomial_output_unique[['distance']] < 150000 & binomial_output_unique[['distance']] > 0 & binomial_output_unique[['r']] < 0, c('distance', 'category'), ], 
    qtl_overlap_unique[qtl_overlap_unique[['distance']] < 150000 & qtl_overlap_unique[['distance']] > 0 & qtl_overlap_unique[['sign']] < 0, c('distance', 'category'), ], 
    hybrid_output_unique[hybrid_output_unique[['distance']] < 150000 & hybrid_output_unique[['distance']] > 0 & hybrid_output_unique[['r']] < 0, c('distance', 'category'), ]), 
  mapping = aes(
    x = distance, 
    fill = category
  )
) + 
  geom_density(alpha = 0.5) +
  xlab('Distance between region and gene') + 
  ylab('Density') + 
  ggtitle('Distance between region and gene\nacross different methods\nfor negative associations') + 
  scale_fill_manual(values = list('binomial' = '#BEAED4', 'pseudobulk' = '#7FC97F', 'qtl_overlap' = '#386CB0', 'hybrid' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_region_direction_distances_neg


# merge pseudobulk and binominal
pseudobulk_vs_binomial <- merge(x = pseudobulk_output_unique, y = binomial_output_unique, by = 'r2g')
# check pseudobulk and scenic
pseudobulk_vs_scenic_unfiltered <- merge(x = pseudobulk_output_unique, y = scenic_output_unique_unfiltered, by = 'r2g')
pseudobulk_vs_scenic <- merge(x = pseudobulk_output_unique, y = scenic_output_unique, by = 'r2g')
# check pseudobulk and qtls
pseudobulk_vs_qtl <- merge(x = pseudobulk_output_unique, y = qtl_overlap_unique, by = 'r2g')
# and binomial vs scenic
scenic_unfiltered_vs_binomial <- merge(x = scenic_output_unique_unfiltered, y = binomial_output_unique, by = 'r2g')
scenic_vs_binomial <- merge(x = scenic_output_unique, y = binomial_output_unique, by = 'r2g')
# scenic vs qtl
scenic_unfiltered_vs_qtl <- merge(x = scenic_output_unique_unfiltered, y = qtl_overlap_unique, by = 'r2g')
scenic_vs_qtl <- merge(x = scenic_output_unique, y = qtl_overlap_unique, by = 'r2g')
# binomial vs qtl
binomial_vs_qtl <- merge(x = binomial_output_unique, y = qtl_overlap_unique, by = 'r2g')
# filtered vs unfiltered scenic
scenic_unfiltered_vs_scenic_filtered <- merge(x = scenic_output_unique, y = scenic_output_unique_unfiltered, by = 'r2g')
# and vs the hybrid approach
pseudobulk_vs_hybrid <- merge(x = pseudobulk_output_unique, y = hybrid_output_unique, by = 'r2g')
hybrid_vs_binomial <- merge(x = hybrid_output_unique, y = binomial_output_unique, by = 'r2g')
hybrid_vs_scenic_unfiltered <- merge(x = hybrid_output_unique, y = scenic_output_unique_unfiltered, by = 'r2g')
hybrid_vs_scenic <- merge(x = hybrid_output_unique, y = scenic_output_unique, by = 'r2g')
hybrid_vs_qtl <- merge(x = hybrid_output_unique, y = qtl_overlap_unique, by = 'r2g')

# get the concordances
con_pseudo_bino <- nrow(pseudobulk_vs_binomial[sign(pseudobulk_vs_binomial[['zscore']]) == sign(pseudobulk_vs_binomial[['meta_z']]), ]) / nrow(pseudobulk_vs_binomial)
# [1] 0.9469027
con_pseudo_sce <- nrow(pseudobulk_vs_scenic[sign(pseudobulk_vs_scenic[['zscore']]) == sign(pseudobulk_vs_scenic[['rho_R2G']]), ]) / nrow(pseudobulk_vs_scenic)
# [1] 0.9518248
con_pseudo_sceun <- nrow(pseudobulk_vs_scenic_unfiltered[sign(pseudobulk_vs_scenic_unfiltered[['zscore']]) == sign(pseudobulk_vs_scenic_unfiltered[['rho_R2G']]), ]) / nrow(pseudobulk_vs_scenic_unfiltered)
# [1] 0.8525721
con_pseudo_qtl <- nrow(pseudobulk_vs_qtl[sign(pseudobulk_vs_qtl[['zscore']]) == sign(pseudobulk_vs_qtl[['sign']]), ]) / nrow(pseudobulk_vs_qtl)
# [1] 0.9861751
con_sceun_bino <- nrow(scenic_unfiltered_vs_binomial[sign(scenic_unfiltered_vs_binomial[['meta_z']]) == sign(scenic_unfiltered_vs_binomial[['rho_R2G']]), ]) / nrow(scenic_unfiltered_vs_binomial)
# [1] 0.8156997
con_sce_bino <- nrow(scenic_vs_binomial[sign(scenic_vs_binomial[['meta_z']]) == sign(scenic_vs_binomial[['rho_R2G']]), ]) / nrow(scenic_vs_binomial)
# [1] 1
con_sce_qtl <- nrow(scenic_vs_qtl[sign(scenic_vs_qtl[['sign']]) == sign(scenic_vs_qtl[['rho_R2G']]), ]) / nrow(scenic_vs_qtl)
# [1] 0.8723776
con_sceun_qtl <- nrow(scenic_unfiltered_vs_qtl[sign(scenic_unfiltered_vs_qtl[['sign']]) == sign(scenic_unfiltered_vs_qtl[['rho_R2G']]), ]) / nrow(scenic_unfiltered_vs_qtl)
# [1] 0.7822823
con_bino_qtl <- nrow(binomial_vs_qtl[sign(binomial_vs_qtl[['meta_z']]) == sign(binomial_vs_qtl[['sign']]), ]) / nrow(binomial_vs_qtl)
# [1] 0.7777778
con_sce_sceun <- nrow(scenic_unfiltered_vs_scenic_filtered[sign(scenic_unfiltered_vs_scenic_filtered[['rho_R2G.x']]) == sign(scenic_unfiltered_vs_scenic_filtered[['rho_R2G.y']]), ]) / nrow(scenic_unfiltered_vs_scenic_filtered)
# [1] 1

# get the concordances
con_pseudo_vs_hyb <- nrow(pseudobulk_vs_hybrid[sign(pseudobulk_vs_hybrid[['zscore.x']]) == sign(pseudobulk_vs_hybrid[['zscore.y']]), ]) / nrow(pseudobulk_vs_hybrid)
# [1] 1
con_hyb_bino <- nrow(hybrid_vs_binomial[sign(hybrid_vs_binomial[['zscore']]) == sign(hybrid_vs_binomial[['meta_z']]), ]) / nrow(hybrid_vs_binomial)
# [1] 1
con_hyb_sce <- nrow(hybrid_vs_scenic[sign(hybrid_vs_scenic[['zscore']]) == sign(hybrid_vs_scenic[['rho_R2G']]), ]) / nrow(hybrid_vs_scenic)
# [1] 0.9615385
con_hyb_sceun <- nrow(hybrid_vs_scenic_unfiltered[sign(hybrid_vs_scenic_unfiltered[['zscore']]) == sign(hybrid_vs_scenic_unfiltered[['rho_R2G']]), ]) / nrow(hybrid_vs_scenic_unfiltered)
# [1] 0.7994269
con_hyb_qtl <- nrow(hybrid_vs_qtl[sign(hybrid_vs_qtl[['zscore']]) == sign(hybrid_vs_qtl[['sign']]), ]) / nrow(hybrid_vs_qtl)
# [1] 0.9230769


# in a pairwise table
pairwise_correlation_table <- data.frame(
  'method1' = c('pseudobulk', 'pseudobulk', 'pseudobulk', 'pseudobulk', 'pseudobulk', 'pseudobulk',
                'scenic', 'scenic', 'scenic', 'scenic', 'scenic', 'scenic', 
                'binomial', 'binomial', 'binomial', 'binomial', 'binomial', 'binomial',
                'QTL', 'QTL', 'QTL', 'QTL', 'QTL', 'QTL', 
                'scenic uf', 'scenic uf', 'scenic uf', 'scenic uf', 'scenic uf', 'scenic uf', 
                'hybrid', 'hybrid', 'hybrid', 'hybrid', 'hybrid', 'hybrid'), 
  'method2' = c('pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid'), 
  'concordance' = c(1, con_pseudo_sce, con_pseudo_bino, con_pseudo_qtl, con_pseudo_sceun, con_pseudo_vs_hyb, 
                    con_pseudo_sce, 1, con_sce_bino, con_sce_qtl, con_sce_sceun, con_hyb_sceun, 
                    con_pseudo_bino, con_sce_bino, 1, con_bino_qtl, con_sceun_bino, con_hyb_bino, 
                    con_pseudo_qtl, con_sce_qtl, con_bino_qtl, 1, con_sceun_qtl, con_hyb_qtl, 
                    con_pseudo_sceun, con_sce_sceun, con_sceun_bino, con_sceun_qtl, 1, con_hyb_sceun, 
                    con_pseudo_vs_hyb, con_hyb_sce, con_hyb_bino, con_hyb_qtl, con_hyb_sceun, 1)
)
# round so it stays readable
pairwise_correlation_table[['concordance']] <- round(pairwise_correlation_table[['concordance']], digits = 2)
# make plot
p_con <- create_confusion_matrix(pairwise_correlation_table, truth_column = 'method1', prediction_column = 'method2', freq_column = 'concordance', premade_table = T, truth_column_label = 'method 1', prediction_column_label = 'method 2') +
  ggtitle('concordance of effect sizes\nin CRE detection methods')
p_con + theme(legend.position="none")

nrow(pseudobulk_vs_binomial)
# [1] 904
nrow(pseudobulk_vs_scenic)
# [1] 1370
nrow(pseudobulk_vs_qtl)
# [1] 1302
nrow(scenic_vs_binomial)
# [1] 478
nrow(scenic_vs_qtl)
# [1] 572
nrow(binomial_vs_qtl)
# [1] 3501
nrow(pseudobulk_vs_scenic_unfiltered)
# [1] 1594
nrow(scenic_unfiltered_vs_scenic_filtered)
# [1] 53796
nrow(scenic_unfiltered_vs_binomial)
# [1] 586
nrow(scenic_unfiltered_vs_qtl)
# [1] 666

# make overlap into a table
pairwise_overlap_table <- data.frame(
  'method1' = c('pseudobulk', 'pseudobulk', 'pseudobulk', 'pseudobulk', 'pseudobulk', 'pseudobulk',
                'scenic', 'scenic', 'scenic', 'scenic', 'scenic', 'scenic', 
                'binomial', 'binomial', 'binomial', 'binomial', 'binomial', 'binomial',
                'QTL', 'QTL', 'QTL', 'QTL', 'QTL', 'QTL', 
                'scenic uf', 'scenic uf', 'scenic uf', 'scenic uf', 'scenic uf', 'scenic uf', 
                'hybrid', 'hybrid', 'hybrid', 'hybrid', 'hybrid', 'hybrid'), 
  'method2' = c('pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid', 
                'pseudobulk', 'scenic', 'binomial', 'QTL', 'scenic uf', 'hybrid'), 
  'overlapping' = c(nrow(pseudobulk_output_unique), nrow(pseudobulk_vs_scenic), nrow(pseudobulk_vs_binomial), nrow(pseudobulk_vs_qtl), nrow(pseudobulk_vs_scenic_unfiltered), nrow(pseudobulk_vs_hybrid), 
                    nrow(pseudobulk_vs_scenic), nrow(scenic_output_unique), nrow(scenic_vs_binomial), nrow(scenic_vs_qtl), nrow(scenic_unfiltered_vs_scenic_filtered), nrow(hybrid_output_unique),  
                    nrow(pseudobulk_vs_binomial), nrow(scenic_vs_binomial), nrow(binomial_output_unique), nrow(binomial_vs_qtl), nrow(scenic_unfiltered_vs_binomial), nrow(hybrid_vs_binomial), 
                    nrow(pseudobulk_vs_qtl), nrow(scenic_vs_qtl), nrow(binomial_vs_qtl), nrow(qtl_overlap_unique), nrow(scenic_unfiltered_vs_qtl), nrow(hybrid_vs_qtl),  
                    nrow(pseudobulk_vs_scenic_unfiltered), nrow(scenic_unfiltered_vs_scenic_filtered), nrow(scenic_unfiltered_vs_binomial), nrow(scenic_unfiltered_vs_qtl), nrow(scenic_output_unique_unfiltered), nrow(hybrid_vs_scenic_unfiltered), 
                    nrow(pseudobulk_vs_hybrid), nrow(hybrid_output_unique), nrow(hybrid_vs_binomial), nrow(hybrid_vs_qtl), nrow(hybrid_vs_scenic_unfiltered), nrow(hybrid_output_unique))
)

# write these tables
write.table(pseudobulk_vs_scenic, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_scenic.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(pseudobulk_vs_binomial, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_binomial.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(pseudobulk_vs_qtl, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_qtl.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(pseudobulk_vs_scenic_unfiltered, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_scenic_unfiltered.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(scenic_vs_binomial, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_vs_binomial.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(scenic_vs_qtl, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_vs_qtl.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(binomial_vs_qtl, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/binomial_vs_qtl.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(scenic_unfiltered_vs_binomial, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_unfiltered_vs_binomial.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(scenic_unfiltered_vs_qtl, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_unfiltered_vs_qtl.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
# make checksums
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_scenic.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_binomial.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_qtl.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/pseudobulk_vs_scenic_unfiltered.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_vs_binomial.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_vs_qtl.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/binomial_vs_qtl.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_unfiltered_vs_binomial.tsv.gz')
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/overlaps/scenic_unfiltered_vs_qtl.tsv.gz')

# make into confusion matrix
p_n_overlaps <- create_confusion_matrix(pairwise_overlap_table, truth_column = 'method1', prediction_column = 'method2', freq_column = 'overlapping', premade_table = T, truth_column_label = 'method 1', prediction_column_label = 'method 2') +
  ggtitle('Overlapping region-gene pairs\nin CRE detection methods')
p_n_overlaps + theme(legend.position="none")

# plot them as well
plot_grid(
  plot_concondance(pseudobulk_vs_binomial, 'zscore', 'meta_z') + ggtitle('Effects of pseudobulk vs binomial\nCRE detection') + xlab('Pseudobulk Z-score') + ylab('Binomial model Z-score'), 
  plot_concondance(pseudobulk_vs_scenic, 'zscore', 'rho_R2G') + ggtitle('Effects of pseudobulk vs SCENIC+ CRE\ndetection') + xlab('Pseudobulk Z-score') + ylab('SCENIC+ R2G Rho'), 
  plot_concondance(scenic_vs_binomial, 'rho_R2G', 'meta_z') + ggtitle('Effects of SCENIC+ vs binomial\nCRE detection') + xlab('SCENIC+ R2G Rho') + ylab('Binomial model Z-score'), 
  plot_concondance(pseudobulk_vs_hybrid, 'z_from_p.x', 'z_from_p.y') + ggtitle('Effects of pseudobulk vs hybrid\nCRE detection') + xlab('Pseudobulk Z-score') + ylab('Hybrid method Z-score')
)

plot_grid(
  plot_concondance(pseudobulk_vs_binomial, 'zscore', 'meta_z') + ggtitle('Effects of pseudobulk vs binomial\nCRE detection') + xlab('Pseudobulk Z-score') + ylab('Binomial model Z-score'), 
  plot_concondance(pseudobulk_vs_scenic_unfiltered, 'zscore', 'rho_R2G') + ggtitle('Effects of pseudobulk vs SCENIC+ CRE\ndetection') + xlab('Pseudobulk Z-score') + ylab('SCENIC+ R2G Rho'), 
  plot_concondance(scenic_vs_binomial, 'rho_R2G', 'meta_z') + ggtitle('Effects of SCENIC+ vs binomial\nCRE detection') + xlab('SCENIC+ R2G Rho') + ylab('Binomial model Z-score')
)
