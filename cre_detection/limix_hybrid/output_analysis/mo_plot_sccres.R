#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_plot_sccres.R
# Function: 
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
library(grid)


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



read_pseudobulk_cre_output_per_celltype <- function(pseudobulk_output_folder, cell_types=NULL, filename_output='qtl_results_all.txt.gz', significance_column='empirical_feature_p_value', significance_cutoff=0.05, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', add_global_nominal_threshold=F, add_local_nominal_threshold=F, global_nominal_threshold_column_to_add='pval_nominal_threshold_global', local_nominal_threshold_column_to_add='pval_nominal_threshold_local', alpha_column='alpha_param', beta_column='beta_param', nominal_p_column='p_value', filter_alpha=T, alpha_min=.2, alpha_max=5, filter_significance=T, pad_columns=T) {
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
    # check if the file exists
    if (file.exists(cell_type_output_loc)) {
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
    else {
      warning(paste('folder exists at', cell_type_output_loc, 'but no file is there'))
    }
  }
  if (pad_columns) {
    # get all the columns we have
    columns_unique <- unique(as.vector(unlist(lapply(output_per_celltype, colnames))))
    # check each output
    for (ct in names(output_per_celltype)) {
      # extract that table
      ct_output <- output_per_celltype[[ct]]
      # check if we are missing any columns
      missing_columns <- setdiff(columns_unique, colnames(ct_output))
      # add those columns
      for (missing_column in missing_columns) {
        ct_output[[missing_column]] <- NA
      }
      # now make sure they are in the same order always
      ct_output <- ct_output[, ..columns_unique]
      # and put back in the list
      output_per_celltype[[ct]] <- ct_output
    }
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
    annotate("text", x = max_sig_z_d1 * -0.70 , y = max_sig_z_d2 * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold') +
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
  color_coding_dict[['all']] <- 'gray'
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
    color_coding_dict[[paste(cell_type, 'shared')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'ps')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, 'sc')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "black"))(100)[pct_whitening]
  }
  # general
  color_coding_dict[['AI']] <- 'darkblue'
  color_coding_dict[['NI']] <- 'darkred'
  color_coding_dict[['Actively Inflamed']] <- 'darkblue'
  color_coding_dict[['Non-Inflamed']] <- 'darkred'
  return(color_coding_dict)
}


get_label_dict <- function(){
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




#' Plot Sharing of  Genes per Cell Type
#'
#' This function plots the sharing of differentially genes across cell types using an UpSet plot. It allows the use of custom label and color dictionaries.
#'
#' @param genes_per_ct list with the genes for each cell type
#' @param use_label_dict A logical value indicating whether to use a custom label dictionary for renaming cell types. Default is TRUE.
#' @param use_color_dict A logical value indicating whether to use a custom color dictionary for cell types. Default is TRUE.
#' @param n_intersects value describing how many intersections to plot, default is all
#' @param use_this_color_dict list with colours for each category, if not supplied, environment default is used
#' @return An UpSet plot showing the sharing of DE genes across cell types.
#'
plot_sharing_per_celltype <- function(genes_per_ct, use_label_dict=T, use_color_dict=T, n_intersects=NA, use_this_color_dict=NULL){
  # rename labels if requested
  if (use_label_dict) {
    names(genes_per_ct) <- remap_with_label_dict(names(genes_per_ct))
  }
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
    cell_type_colours <- NULL
    # use the supplied colour coding dict
    if (!is.null(use_this_color_dict)) {
      cell_type_colours <- use_this_color_dict
    }
    # otherwise the default one
    else {
      cell_type_colours <- get_color_coding_dict()
    }
    
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


get_group_proportions <- function(named_list_of_dfs, column_to_get_proportions_from='screen', specific_trait_expression=NULL, specific_trait_name=NULL) {
  # we'll save per df
  results_per_group <- list()
  # we'll check each df
  for (group_name in names(named_list_of_dfs)) {
    # we'll need the numbers for the traits
    group_totals <- NULL
    # check a specific trait if requested
    if (!is.null(specific_trait_name) & !is.null(specific_trait_expression)) {
      # get true/false vector for that expression
      is_of_trait <- grepl(specific_trait_expression, named_list_of_dfs[[group_name]][[column_to_get_proportions_from]])
      # turn into a string
      is_trait_string <- ifelse (is_of_trait, paste('is', specific_trait_name), paste('is not', specific_trait_name))
      # make factor
      is_trait_factor <- factor(is_trait_string, levels = c(paste('is', specific_trait_name), paste('is not', specific_trait_name)))
      # get the totals from that
      group_totals <- data.frame(table(is_trait_string))
    }
    else {
      # otherwise get the totals of each group
      group_totals <- data.frame(table(named_list_of_dfs[[group_name]][[column_to_get_proportions_from]]))
    }
    # rename the columns
    colnames(group_totals) <- c('category', 'n')
    # add the fraction now as well
    group_totals[['frac']] <- group_totals[['n']] / sum(group_totals[['n']])
    # finally add the group as well
    group_totals[['group']] <- group_name
    # place in the list
    results_per_group[[group_name]] <- group_totals
  }
  # merge all the tables
  results_all <- do.call('rbind', results_per_group)
  # make sure the 'nots' are always last
  categories <- as.character(unique(results_all[['category']]))
  # get the ones that have 'not' in their name
  categories_not <- categories[grep('^is not ', categories)]
  # and the ones that are not not 'not'
  categories_other <- setdiff(categories, categories_not)
  # and make that the order
  results_all[['category']] <- factor(results_all[['category']], levels=(c(categories_not, categories_other)))
  results_all <- results_all[order(results_all[['category']]), ]
  return(results_all)
}


add_nominal_p_value_cutoff <- function(cre_sumstats, cell_type_column='cell_type', feature_column='feature_id', emperical_p_column='empirical_feature_p_value', emperical_p_cutoff=0.05, nominal_p_column='p_value', nominal_p_cutoff=0.05) {
  # convert
  cre_sumstats <- data.frame(cre_sumstats)
  # check each cell type
  cres_cutoff_l <- list()
  # check each cell type
  for (cell_type in unique(cre_sumstats[[cell_type_column]])) {
    # subset to cell type
    cre_sumstats_ct <- cre_sumstats[!is.na(cre_sumstats[[cell_type_column]]) & cre_sumstats[[cell_type_column]] == cell_type, ]
    # order by the nominal p value
    cre_sumstats_ct <- cre_sumstats_ct[order(cre_sumstats_ct[[nominal_p_column]], decreasing = F), ]
    # keep only top effects
    cre_sumstats_ct_top <- cre_sumstats_ct[!duplicated(cre_sumstats_ct[[feature_column]]), ]
    # then filter on the significant ones
    if (is.numeric(emperical_p_cutoff)) {
      # if numeric
      cre_sumstats_ct_top_significant <- cre_sumstats_ct_top[cre_sumstats_ct_top[[emperical_p_column]] < emperical_p_cutoff, ]
    }
    else if(is.logical(emperical_p_cutoff)) {
      # if logical
      cre_sumstats_ct_top_significant <- cre_sumstats_ct_top[cre_sumstats_ct_top[[emperical_p_column]] == emperical_p_cutoff, ]
    }
    else {
      stop('emperical cutoff must be logical or numeric')
    }
    # get the max p that is still significant
    cre_sumstats_ct_top_significant_top_p <- max(cre_sumstats_ct_top_significant[[nominal_p_column]])
    # add that information
    cre_sumstats_ct[['nominal_p_cutoff']] <- cre_sumstats_ct_top_significant_top_p
    # and easy cutoff info
    cre_sumstats_ct[['significant_nominal_cutoff']] <- cre_sumstats_ct[[nominal_p_column]] <= nominal_p_cutoff
    # add that to the list
    cres_cutoff_l[[cell_type]] <- cre_sumstats_ct
  }
  # merge all
  cres_cutoff <- do.call('rbind', cres_cutoff_l)
  # convert
  cres_cutoff <- data.table(cres_cutoff)
  return(cres_cutoff)
}


split_strings_with_max_length <- function(strings, max_length=15) {
  strings_split <- sapply(strings, FUN = function(x){
    split_string_with_max_length(x, max_length)
  })
  strings_split <- as.vector(strings_split)
  return(strings_split)
}

#' plot in a barplot the number of DE genes per cell type
#' 
#' @param egene_numbers the table containing the number of eenes per cell type
#' @param celltype_column which column of the supplied table has the celltype name
#' @param number_column which column of the supplied table has the number of DE genes
#' @param pointless draw the ticks on the bottom of the plot?
#' @param legendless draw the legend?
#' @param ylim the range for the y-axis
#' @param paper_style use the paper style with more whitespace
#' @param angle_labels angle the x labels 90 degrees?
#' @param use_distinct_colours use the distinct colour palette instead of ggplot defaults (optional, default T)
#' @param use_sampling randomly sample colours from the distinct colour palette (optional, default F)
#' @param color_indices supply specific indices relating to colours to use from the distinct colour palette. Use 'get_available_colours_grid' to get the colours and their indices (optional, unused by default)
#' @param use_label_dict change names to a more pleasing format
#' @param make_labels_multiline make labels multiline if they are long
#' @param x_order to order the values on the x axis
#' @returns a barplot with the number of eGenes per cell type
#' 
egene_numbers_to_plot <- function(egene_numbers, celltype_column='cell_type', number_column='nr', pointless=F, legendless=F, ylim=NULL, paper_style=T, angle_labels=T, use_distinct_colours=T, use_sampling=F, color_indices=NULL, use_color_dict=F, use_label_dict=T, split_long_labels_to_lines=T, x_order=NULL) {
  # remap names if requested
  if (use_label_dict) {
    egene_numbers[[celltype_column]] <- remap_with_label_dict(egene_numbers[[celltype_column]])
  }
  # make them multiline
  if (split_long_labels_to_lines) {
    egene_numbers[[celltype_column]] <- split_strings_with_max_length(egene_numbers[[celltype_column]])
  }
  # order the x axis if requested
  if (!is.null(x_order)) {
    egene_numbers[[celltype_column]] <- factor(egene_numbers[[celltype_column]], levels = x_order)
  }
  # start plot
  p <- ggplot(mapping = aes(x = egene_numbers[[celltype_column]], y = egene_numbers[[number_column]], fill = egene_numbers[[celltype_column]])) +
    geom_bar(stat='identity') + xlab(celltype_column) + ylab(number_column)
  # options
  if (use_distinct_colours) {
    # get the unique possible assigments
    possible_assignments <- unique(egene_numbers[[celltype_column]])
    # put into a list
    colour_mapping <- roycols::get_color_list(possible_assignments)
    # add to plot
    p <- p + scale_fill_manual(values = colour_mapping)
  }
  if (use_color_dict) {
    # get the color coding dict
    colour_mapping <- get_color_coding_dict()
    # also rename here then
    if (split_long_labels_to_lines) {
      names(colour_mapping) <- split_strings_with_max_length(names(colour_mapping))
    }
    # add to plot
    p <- p + scale_fill_manual(values = colour_mapping)
  }
  if(!is.null(ylim)){
    p <- p + ylim(ylim)
  }
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p)
}


#' Function for Rb analysis
#'
#' @param b1 Beta from first dataset.
#' @param se1 Standard error of beta from first dataset.
#' @param b2 Beta from second dataset.
#' @param se2 Standard error of beta from second dataset.
#' @param theta Variable representing sample overlap between two datasets. Should be set 0 if no sample overlap.
#'
#' @return Data frame with Rb, SE(Rb) and corresponding P-value.
#' @export
#'
#' @note This function is slightly adapted from the script shared by Ting Qi.
#'
#' @examples
calcu_cor_true <- function(b1, se1, b2, se2, theta=1) {
  idx <- which(is.infinite(b1) | is.infinite(b2) | is.infinite(se1) | is.infinite(se2))
  if (length(idx) > 0) {
    b1 <- b1[-idx]
    se1 <- se1[-idx]
    b2 <- b2[-idx]
    se2 <- se2[-idx]
    theta <- theta[-idx]
  }
  
  var_b1 <- var(b1, na.rm = T) - mean(se1^2, na.rm = T)
  var_b2 <- var(b2, na.rm = T) - mean(se2^2, na.rm = T)
  if (var_b1 < 0) {
    var_b1 <- var(b1, na.rm = T)
  }
  if (var_b2 < 0) {
    var_b2 <- var(b2, na.rm = T)
  }
  cov_b1_b2 <- cov(b1, b2, use = "complete.obs") - mean(theta, na.rm = T) * sqrt(mean(se1^2, na.rm = T) * mean(se2^2, na.rm = T))
  r <- cov_b1_b2 / sqrt(var_b1 * var_b2)
  
  r_jack <- c()
  n <- length(b1)
  for (k in 1:n) {
    b1_jack <- b1[-k]
    se1_jack <- se1[-k]
    var_b1_jack <- var(b1_jack, na.rm = T) - mean(se1_jack^2, na.rm = T)
    b2_jack <- b2[-k]
    se2_jack <- se2[-k]
    var_b2_jack <- var(b2_jack, na.rm = T) - mean(se2_jack^2, na.rm = T)
    if (var_b1_jack < 0) {
      var_b1_jack <- var(b1_jack, na.rm = T)
    }
    if (var_b2_jack < 0) {
      var_b2_jack <- var(b2_jack, na.rm = T)
    }
    theta_jack <- theta[-k]
    cov_e1_jack_e2_jack <- mean(theta_jack, na.rm = T) * sqrt(mean(se1_jack^2, na.rm = T) * mean(se2_jack^2, na.rm = T))
    cov_b1_b2_jack <- cov(b1_jack, b2_jack, use = "complete.obs") - cov_e1_jack_e2_jack
    r_tmp <- cov_b1_b2_jack / sqrt(var_b1_jack * var_b2_jack)
    r_jack <- c(r_jack, r_tmp)
  }
  r_mean <- mean(r_jack, na.rm = T)
  idx <- which(is.na(r_jack))
  if (length(idx) > 0) {
    se_r <- sqrt((n - 1) / n * sum((r_jack[-idx] - r_mean)^2))
  } else {
    se_r <- sqrt((n - 1) / n * sum((r_jack - r_mean)^2))
  }
  
  p <- pchisq((r / se_r)**2, df = 1, lower.tail = FALSE)
  
  res <- cbind(r, se_r, p)
  return(res)
}


match_qtl_tables <- function(qtl_replicating, 
                             qtl_replicating_in, 
                             qtl_replicating_name='sccre', 
                             qtl_replicating_in_name='pscre', 
                             cell_type_column='cell_type', 
                             feature_column='feature_id', 
                             variant_column='snp_id', 
                             emperical_p_column='empirical_feature_p_value', 
                             nominal_p_column='p_value', 
                             beta_column='beta', 
                             se_column='beta_se', 
                             zscore_column='zscore', 
                             allele_column='assessed_allele', 
                             significance_cutoff_emperical=0.05, 
                             significance_cutoff_nominal=0.05) {
  # we'll store the results in a list
  qtl_overlaps <- list()
  # convert to datatables
  qtl_replicating <- data.frame(qtl_replicating)
  qtl_replicating_in <- data.frame(qtl_replicating_in)
  # check each cell type
  cell_types_common <- as.character(intersect(unique(as.character(qtl_replicating[[cell_type_column]])), unique(as.character(qtl_replicating_in[[cell_type_column]]))))
  for (cell_type in cell_types_common) {
    qtl_replicating_ct <- qtl_replicating[!is.na(qtl_replicating[[cell_type_column]]) & qtl_replicating[[cell_type_column]] == cell_type, ]
    qtl_replicating_in_ct <- qtl_replicating_in[!is.na(qtl_replicating_in[[cell_type_column]]) & qtl_replicating_in[[cell_type_column]] == cell_type, ]
    # sort the replicating by the strongest effect for each feature
    qtl_replicating_ct <- qtl_replicating_ct[
      order(qtl_replicating_ct[[nominal_p_column]], decreasing = F), 
    ]
    # keep only the first entry
    qtl_replicating_ct_top <- qtl_replicating_ct[!duplicated(qtl_replicating_ct[[feature_column]]), ]
    # and keep only significant
    if (is.numeric(significance_cutoff_emperical)) {
      # if numeric
      qtl_replicating_ct_top <- qtl_replicating_ct_top[qtl_replicating_ct_top[[emperical_p_column]] < significance_cutoff_emperical, ]
    }
    else if(is.logical(significance_cutoff_emperical)) {
      # if logical
      qtl_replicating_ct_top <- qtl_replicating_ct_top[qtl_replicating_ct_top[[emperical_p_column]] == significance_cutoff_emperical, ]
    }
    else {
      stop('emperical cutoff must be logical or numeric')
    }
    # add QTL in both
    qtl_replicating_ct_top[['qtl']] <- paste(qtl_replicating_ct_top[[feature_column]], qtl_replicating_ct_top[[variant_column]], sep = '_')
    qtl_replicating_in_ct[['qtl']] <- paste(qtl_replicating_in_ct[[feature_column]], qtl_replicating_in_ct[[variant_column]], sep = '_')
    # subset both to columns we need
    qtl_replicating_ct_top_merge <- qtl_replicating_ct_top[, c('qtl', feature_column, variant_column, emperical_p_column, nominal_p_column, beta_column, se_column, zscore_column)]
    qtl_replicating_ct_in_merge <- qtl_replicating_in_ct[, c('qtl', emperical_p_column, nominal_p_column, beta_column, se_column, zscore_column)]
    # set the column names
    colnames(qtl_replicating_ct_top_merge) <- c('qtl', 'feature_id', 'variant_id', paste(qtl_replicating_name, c(emperical_p_column, nominal_p_column, beta_column, se_column, zscore_column), sep = '_'))
    colnames(qtl_replicating_ct_in_merge) <- c('qtl', paste(qtl_replicating_in_name, c(emperical_p_column, nominal_p_column, beta_column, se_column, zscore_column), sep = '_'))
    # add allele column if present
    if (!is.null(allele_column)) {
      qtl_replicating_ct_top_merge[[paste(qtl_replicating_name, allele_column, sep = '_')]] <- qtl_replicating_ct_top[[allele_column]]
      qtl_replicating_ct_in_merge[[paste(qtl_replicating_in_name, allele_column, sep = '_')]] <- qtl_replicating_in_ct[[allele_column]]
    }
    else {
      qtl_replicating_ct_top_merge[[paste(qtl_replicating_name, 'assessed_allele', sep = '_')]] <- 'X'
      qtl_replicating_ct_in_merge[[paste(qtl_replicating_in_name, 'assessed_allele', sep = '_')]] <- 'X'
    }
    # merge the two
    qtl_ct_both <- unique(merge(unique(qtl_replicating_ct_top_merge), unique(qtl_replicating_ct_in_merge), by = 'qtl'))
    # put in the list
    qtl_overlaps[[paste0(cell_type, '_vs_', cell_type)]] <- qtl_ct_both
  }
  return(qtl_overlaps)
}


get_stats_qtl_matches <- function(table_per_match, name_split='_vs_', 
                                  significance_column_replicating_emperical='fdr', 
                                  significance_column_replicating_nominal='p',
                                  significance_column_replicating_in_emperical='fdr',
                                  significance_column_replicating_in_nominal='p',
                                  significance_cutoff_replicating_emperical=0.05, 
                                  significance_cutoff_replicating_nominal=0.05, 
                                  significance_cutoff_replicating_in_emperical=0.05, 
                                  significance_cutoff_replicating_in_nominal=0.05, 
                                  score_column_replicating='zscore', 
                                  score_column_replicating_in='zscore', 
                                  beta_column_replicating='beta', 
                                  beta_column_replicating_in='beta', 
                                  se_column_replicating='se', 
                                  se_column_replicating_in='se', 
                                  allele_column_replicating='assessed_allele', 
                                  allele_column_replicating_in='assessed_allele', 
                                  cor_method='spearman') {
  # prepare a list of tables of all the comparisons
  stats_per_combination <- list()
  # check each of the tables
  for (combination in names(table_per_match)) {
    print(combination)
    # extract the two names
    replicating_name <- strsplit(combination, name_split)[[1]][1]
    replicating_in_name <- strsplit(combination, name_split)[[1]][2]
    # extract the table
    replication_table <- table_per_match[[combination]]
    # init variables
    replicated_emperical_emperical_in_replicating <- NA
    nominal_to_nominal_zscore_correlation <- NA
    emperical_to_nominal_zscore_correlation <- NA
    emperical_to_emperical_zscore_correlation <- NA
    nominal_to_nominal_concordance <- NA
    emperical_to_nominal_concordance <- NA
    emperical_to_emperical_concordance <- NA
    nominal_to_nominal_rb <- NA
    emperical_to_nominal_rb <- NA
    emperical_to_all_rb <- NA
    emperical_to_emperical_rb <- NA
    # if they are the same, we don't have to do much
    if (replicating_name == replicating_in_name) {
      # we check how many are significant
      n_sig <- NULL
      if (isTRUE(significance_cutoff_replicating_emperical)) {
        n_sig <- nrow(replication_table[
          replication_table[[paste0(replicating_name, '_', significance_column_replicating_emperical, '.x')]] == T, ])
      } else {
        n_sig <- nrow(replication_table[
          replication_table[[paste0(replicating_name, '_', significance_column_replicating_emperical, '.x')]] < significance_cutoff_replicating_emperical, ])
      }
      # and then all of these will be the same
      replicated_emperical_emperical_in_replicating <- n_sig
      replicated_emperical_nominal_in_replicating <- n_sig
      replicated_emperical_tested_in_replicating <- n_sig
      replicated_emperical_nominal_in_replicating <- n_sig
      # these will all be one
      nominal_to_nominal_zscore_correlation <- 1
      emperical_to_nominal_zscore_correlation <- 1
      emperical_to_emperical_zscore_correlation <- 1
      nominal_to_nominal_concordance <- 1
      emperical_to_nominal_concordance <- 1
      emperical_to_emperical_concordance <- 1
      nominal_to_nominal_rb <- 1
      emperical_to_nominal_rb <- 1
      emperical_to_all_rb <- 1
      emperical_to_emperical_rb <- 1
    }
    else {
      # check how many significant in replicating, are present in the other
      replicated_emperical_tested_in_replicating <- NULL
      # using the boolean value in the column if it was a boolean
      if (isTRUE(significance_cutoff_replicating_emperical)) {
        replicated_emperical_tested_in_replicating <- nrow(replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] == T, ])
      } else {
        # otherwise use it as a cutoff
        replicated_emperical_tested_in_replicating <- nrow(replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] < significance_cutoff_replicating_emperical, ])
      }
      
      # subset to emperical only in where we are replicating
      replication_table_emp_to_all <- NULL
      # if the significance cutoff is a T/F we don't use it as a cutoff
      if (isTRUE(significance_cutoff_replicating_emperical)) {
        replication_table_emp_to_all <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] == T, 
        ]
      } else {
        replication_table_emp_to_all <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] < significance_cutoff_replicating_emperical, 
        ]
      }
      # and the beta
      beta_replicating_in <- replication_table_emp_to_all[[paste(replicating_in_name, beta_column_replicating_in, sep = '_')]]
      # check which have the same allelic direction
      alleles_differ_locs <- which(replication_table_emp_to_all[[paste(replicating_name, allele_column_replicating, sep = '_')]] != replication_table_emp_to_all[[paste(replicating_in_name, allele_column_replicating_in, sep = '_')]])
      # and flip direction of betas to reflect this
      beta_replicating_in[alleles_differ_locs] <- -1 * beta_replicating_in[alleles_differ_locs]
      # calculate Rb
      emp_to_all_rb_table <- calcu_cor_true(b1 = replication_table_emp_to_all[[paste(replicating_name, beta_column_replicating, sep = '_')]], se1 = replication_table_emp_to_all[[paste(replicating_name, se_column_replicating, sep = '_')]], b2 = beta_replicating_in, se2 = replication_table_emp_to_all[[paste(replicating_in_name, se_column_replicating_in, sep = '_')]])
      emperical_to_all_rb <- emp_to_all_rb_table[1, 'r']
      
      # subset to nominal in both
      replication_table <- replication_table[
        replication_table[[paste(replicating_name, significance_column_replicating_nominal, sep = '_')]] < significance_cutoff_replicating_nominal &
          replication_table[[paste(replicating_in_name, significance_column_replicating_in_nominal, sep = '_')]] < significance_cutoff_replicating_in_nominal, 
      ]
      # get the Z scores
      zscore_replicating <- replication_table[[paste(replicating_name, score_column_replicating, sep = '_')]]
      zscore_replicating_in <- replication_table[[paste(replicating_in_name, score_column_replicating_in, sep = '_')]]
      # and the beta
      beta_replicating_in <- replication_table[[paste(replicating_in_name, beta_column_replicating_in, sep = '_')]]
      # check which have the same allelic direction
      alleles_differ_locs <- which(replication_table[[paste(replicating_name, allele_column_replicating, sep = '_')]] != replication_table[[paste(replicating_in_name, allele_column_replicating_in, sep = '_')]])
      # and flip direction of z score and betas to reflect this
      zscore_replicating_in[alleles_differ_locs] <- -1 * zscore_replicating_in[alleles_differ_locs]
      beta_replicating_in[alleles_differ_locs] <- -1 * beta_replicating_in[alleles_differ_locs]
      # get the correlation of the two
      nominal_to_nominal_zscore_correlation <- cor(zscore_replicating, zscore_replicating_in, method = cor_method)
      # and concordance
      nominal_to_nominal_concordance <- sum(sign(zscore_replicating) == sign(zscore_replicating_in)) / length(zscore_replicating)
      # and the rb
      #nominal_to_nominal_rb <- 1 - mean(abs(((zscore_replicating_in - zscore_replicating) / zscore_replicating)))
      if (nrow(replication_table) > 0) {
        nominal_to_nominal_rb_table <- calcu_cor_true(b1 = replication_table[[paste(replicating_name, beta_column_replicating, sep = '_')]], se1 = replication_table[[paste(replicating_name, se_column_replicating, sep = '_')]], b2 = beta_replicating_in, se2 = replication_table[[paste(replicating_in_name, se_column_replicating_in, sep = '_')]])
        nominal_to_nominal_rb <- nominal_to_nominal_rb_table[1, 'r']
      }
      else {
        nominal_to_nominal_rb <- NA
      }
      # subset the replicating to emperical for replicating
      if (isTRUE(significance_cutoff_replicating_emperical)) {
        # if the cutoff is just a 'T' value, we use that as the filter
        replication_table <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] == T &
            replication_table[[paste(replicating_in_name, significance_column_replicating_in_nominal, sep = '_')]] < significance_cutoff_replicating_in_nominal, 
        ]
      } else {
        replication_table <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] < significance_cutoff_replicating_emperical &
            replication_table[[paste(replicating_in_name, significance_column_replicating_in_nominal, sep = '_')]] < significance_cutoff_replicating_in_nominal, 
        ]
      }
      
      # check how many significant in replicating, are nominally significant in the other
      replicated_emperical_nominal_in_replicating <- nrow(replication_table)
      # get the Z scores
      zscore_replicating <- replication_table[[paste(replicating_name, score_column_replicating, sep = '_')]]
      zscore_replicating_in <- replication_table[[paste(replicating_in_name, score_column_replicating_in, sep = '_')]]
      # and the beta
      beta_replicating_in <- replication_table[[paste(replicating_in_name, beta_column_replicating_in, sep = '_')]]
      # check which have the same allelic direction
      alleles_differ_locs <- which(replication_table[[paste(replicating_name, allele_column_replicating, sep = '_')]] != replication_table[[paste(replicating_in_name, allele_column_replicating_in, sep = '_')]])
      # and flip direction of z score to reflect this
      zscore_replicating_in[alleles_differ_locs] <- -1 * zscore_replicating_in[alleles_differ_locs]
      beta_replicating_in[alleles_differ_locs] <- -1 * beta_replicating_in[alleles_differ_locs]
      # get the correlation of the two
      emperical_to_nominal_zscore_correlation <- cor(zscore_replicating, zscore_replicating_in, method = cor_method)
      # and concordance
      emperical_to_nominal_concordance <- sum(sign(zscore_replicating) == sign(zscore_replicating_in)) / length(zscore_replicating)
      # and the rb
      #emperical_to_nominal_rb <- 1 - mean(abs(((zscore_replicating_in - zscore_replicating) / zscore_replicating)))
      if (nrow(replication_table) > 0) {
        emperical_to_nominal_rb_table <- calcu_cor_true(b1 = replication_table[[paste(replicating_name, beta_column_replicating, sep = '_')]], se1 = replication_table[[paste(replicating_name, se_column_replicating, sep = '_')]], b2 = beta_replicating_in, se2 = replication_table[[paste(replicating_in_name, se_column_replicating_in, sep = '_')]])
        emperical_to_nominal_rb <- emperical_to_nominal_rb_table[1, 'r']
      }
      else {
        emperical_to_nominal_rb <- NA
      }
      
      # subset the replicating to emperical for both
      if (isTRUE(significance_cutoff_replicating_emperical) & isTRUE(significance_cutoff_replicating_in_emperical)) {
        # subset both just on that T value
        replication_table <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] == T &
            replication_table[[paste(replicating_in_name, significance_column_replicating_in_emperical, sep = '_')]] == T, 
        ]
      }
      else if(isTRUE(significance_cutoff_replicating_emperical)) {
        # or only one
        replication_table <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] == T &
            replication_table[[paste(replicating_in_name, significance_column_replicating_in_emperical, sep = '_')]] < significance_cutoff_replicating_in_emperical, 
        ]
      }
      else if(isTRUE(significance_cutoff_replicating_in_emperical)) {
        # or the other
        replication_table <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] < significance_cutoff_replicating_emperical &
            replication_table[[paste(replicating_in_name, significance_column_replicating_in_emperical, sep = '_')]] == T, 
        ]
      }
      else {
        # or on actual values
        replication_table <- replication_table[
          replication_table[[paste(replicating_name, significance_column_replicating_emperical, sep = '_')]] < significance_cutoff_replicating_emperical &
            replication_table[[paste(replicating_in_name, significance_column_replicating_in_emperical, sep = '_')]] < significance_cutoff_replicating_in_emperical, 
        ]
      }
      
      # check how many significant in replicating, are also significant in the other
      replicated_emperical_emperical_in_replicating <- nrow(replication_table)
      # get the Z scores
      zscore_replicating <- replication_table[[paste(replicating_name, score_column_replicating, sep = '_')]]
      zscore_replicating_in <- replication_table[[paste(replicating_in_name, score_column_replicating_in, sep = '_')]]
      # and the beta
      beta_replicating_in <- replication_table[[paste(replicating_in_name, beta_column_replicating_in, sep = '_')]]
      # check which have the same allelic direction
      alleles_differ_locs <- which(replication_table[[paste(replicating_name, allele_column_replicating, sep = '_')]] != replication_table[[paste(replicating_in_name, allele_column_replicating_in, sep = '_')]])
      # and flip direction of z score to reflect this
      zscore_replicating_in[alleles_differ_locs] <- -1 * zscore_replicating_in[alleles_differ_locs]
      beta_replicating_in[alleles_differ_locs] <- -1 * beta_replicating_in[alleles_differ_locs]
      # get the correlation of the two
      emperical_to_emperical_zscore_correlation <- cor(zscore_replicating, zscore_replicating_in, method = cor_method)
      # and concordance
      emperical_to_emperical_concordance <- sum(sign(zscore_replicating) == sign(zscore_replicating_in)) / length(zscore_replicating)
      # and the rb
      #emperical_to_emperical_rb <- 1 - mean(abs(((zscore_replicating_in - zscore_replicating) / zscore_replicating)))
      if (nrow(replication_table) > 0) {
        emperical_to_emperical_rb_table <- calcu_cor_true(b1 = replication_table[[paste(replicating_name, beta_column_replicating, sep = '_')]], se1 = replication_table[[paste(replicating_name, se_column_replicating, sep = '_')]], b2 = beta_replicating_in, se2 = replication_table[[paste(replicating_in_name, se_column_replicating_in, sep = '_')]])
        emperical_to_emperical_rb <- emperical_to_emperical_rb_table[1, 'r']
      }
      else {
        emperical_to_emperical_rb <- NA
      }
      
      
    }
    # put it all in a table
    stats_combination <- data.frame(
      'replicating' = c(replicating_name), 
      'replicating_in' = c(replicating_in_name), 
      'n_sig_tested' = c(replicated_emperical_tested_in_replicating), 
      'n_sig_nom' = c(replicated_emperical_nominal_in_replicating), 
      'n_sig_both' = c(replicated_emperical_emperical_in_replicating), 
      'cor_nom_nom' = c(nominal_to_nominal_zscore_correlation), 
      'cor_emp_nom' = c(emperical_to_nominal_zscore_correlation), 
      'cor_emp_emp' = c(emperical_to_emperical_zscore_correlation), 
      'con_nom_nom' = c(nominal_to_nominal_concordance), 
      'con_emp_nom' = c(emperical_to_nominal_concordance), 
      'con_emp_emp' = c(emperical_to_emperical_concordance), 
      'rb_nom_nom' = c(nominal_to_nominal_rb), 
      'rb_emp_nom' = c(emperical_to_nominal_rb), 
      'rb_emp_emp' = c(emperical_to_emperical_rb), 
      'rb_emp_all' = c(emperical_to_all_rb))
    # put in the list
    stats_per_combination[[combination]] <- stats_combination
  }
  # merge all tables
  stats_combination_all <- do.call('rbind', stats_per_combination)
  return(stats_combination_all)
}


plot_replication_stats <- function(replication_stats, replicating_column='replicating', replicating_in_column='replicating_in', matching_value_to_text=list('n_sig_tested' = 'N sig:'), comparing_value_to_text=list('n_sig_both' = 'N sig both:', 'n_sig_tested' = 'N test both:', 'cor_emp_nom' = 'nom cor:'), colour_column='cor_emp_nom', angle_labels=F, use_label_dict=T, legendless=F, paper_style=T, round_decimals=2, low_color='blue', mid_color='white', high_color='red', mid_value=0, category_plot_order=NULL, text_size=6) {
  # round values if required
  if (!is.null(round_decimals)) {
    for (col in colnames(replication_stats)) {
      if (is.numeric(replication_stats[[col]])) {
        replication_stats[[col]] <- round(replication_stats[[col]], digits = round_decimals)
      }
    }
  }
  # check each row of this table
  replication_stats[['display_text']] <- apply(replication_stats, 1, FUN = function(x) {
    # get the replicating cell type
    replicating <- x[[replicating_column]]
    # and what we are replicating in
    replicating_in <- x[[replicating_in_column]]
    # create the text
    final_text <- NULL
    # check if they are the same
    if (replicating == replicating_in) {
      # then we use this list
      list_to_check <- matching_value_to_text
    }
    else {
      # otherwise we use the other list
      list_to_check <- comparing_value_to_text
    }
    # check each column in the matching list
    for (column_to_check in names(list_to_check)) {
      # make the new piece of text by getting the matching text, and the value from the table
      display_text_part <- paste0(list_to_check[[column_to_check]], as.character(x[[column_to_check]]))
      # check if we have any text
      if (is.null(final_text)) {
        final_text <- display_text_part
      }
      else {
        # if we have, paste together
        final_text <- paste(final_text, display_text_part, sep = '\n')
      }
    }
    return(final_text)
  })
  # remap names if requested
  if (use_label_dict) {
    replication_stats[[replicating_column]] <- remap_with_label_dict(replication_stats[[replicating_column]])
    replication_stats[[replicating_in_column]] <- remap_with_label_dict(replication_stats[[replicating_in_column]])
  }
  # order labels if requested
  if (!is.null(category_plot_order)) {
    # double check if we have the labels
    all_labels <- unique(c(replication_stats[[replicating_column]], replication_stats[[replicating_in_column]]))
    # check if there is any difference between what we have and the categories
    labels_missing <- setdiff(all_labels, category_plot_order)
    # and stop it there are
    if (length(labels_missing) > 0) {
      stop(paste('category_plot_order does not contain all values in categories, missing are:', paste(labels_missing, sep = ',')))
    }
    replication_stats[[replicating_column]] <- factor(replication_stats[[replicating_column]], levels = category_plot_order)
    replication_stats[[replicating_in_column]] <- factor(replication_stats[[replicating_in_column]], levels = category_plot_order)
  }
  # now make the plot (use the !! rlang::sym(variable_name) synthax to be able to use a variable for ggplot)
  p <- ggplot(data=replication_stats, aes(x=!! rlang::sym(replicating_column), y=!! rlang::sym(replicating_in_column), fill=!! rlang::sym(colour_column))) + geom_tile() + scale_fill_gradient2(low=low_color, mid = mid_color, high=high_color, midpoint = mid_value) + geom_text(aes(label=display_text), size = text_size)
  # some options
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p)
}


get_overlapping_peakgene_links <- function(cre_tbl1, cre_tbl2, region_column1='snp_id', region_column2='snp_id', gene_column1='feature_id', gene_column2='feature_id', cell_type_column1='cell_type', cell_type_column2='cell_type', significance_cutoffs1=list('global_significance' = T, 'significant_nominal_cutoff' = T), significance_cutoffs2=list('global_significance' = T, 'significant_nominal_cutoff' = T), set_name1='sc', set_name2='ps') {
  # convert
  cre_tbl1 <- data.frame(cre_tbl1)
  cre_tbl2 <- data.frame(cre_tbl2)
  # save result in list first
  result_per_ct <- list()
  # also store for all cts
  r2g_all <- list('d1' = list(), 'd2' = list())
  # get the cell types
  cell_types_both <- unique(intersect(cre_tbl1[[cell_type_column1]], cre_tbl2[[cell_type_column2]]))
  # then check each of them
  for (ct in cell_types_both) {
    print(ct)
    # subset both on these cell types
    cre_tbl1_ct <- cre_tbl1[!is.na(cre_tbl1[[cell_type_column1]]) & cre_tbl1[[cell_type_column1]] == ct, ]
    cre_tbl2_ct <- cre_tbl2[!is.na(cre_tbl2[[cell_type_column2]]) & cre_tbl2[[cell_type_column2]] == ct, ]
    # now filter both on the significance cutoffs
    for (significance_column in names(significance_cutoffs1)) {
      # extract the cutoff
      significance_cutoff <- significance_cutoffs1[[significance_column]]
      # now based on the type, use to filter
      if (is.numeric(significance_cutoff)) {
        cre_tbl1_ct <- cre_tbl1_ct[cre_tbl1_ct[[significance_column]] < significance_cutoff, ]
      }
      else if (is.logical(significance_cutoff)) {
        cre_tbl1_ct <- cre_tbl1_ct[cre_tbl1_ct[[significance_column]] == significance_cutoff, ]
      }
    }
    # repeat for the other table
    for (significance_column in names(significance_cutoffs2)) {
      # extract the cutoff
      significance_cutoff <- significance_cutoffs2[[significance_column]]
      # now based on the type, use to filter
      if (is.numeric(significance_cutoff)) {
        cre_tbl2_ct <- cre_tbl2_ct[cre_tbl2_ct[[significance_column]] < significance_cutoff, ]
      }
      else if (is.logical(significance_cutoff)) {
        cre_tbl2_ct <- cre_tbl2_ct[cre_tbl2_ct[[significance_column]] == significance_cutoff, ]
      }
    }
    # extract the region-gene pairs for each of these
    r2g1 <- paste(cre_tbl1_ct[[region_column1]], cre_tbl1_ct[[gene_column1]])
    r2g2 <- paste(cre_tbl2_ct[[region_column2]], cre_tbl2_ct[[gene_column2]])
    # now get which are shared and exclusive
    r2g1_only <- setdiff(r2g1, r2g2)
    r2g2_only <- setdiff(r2g2, r2g1)
    r2g_shared <- intersect(r2g1, r2g2)
    # get the sizes
    r2g1_only_n <- length(r2g1_only)
    r2g2_only_n <- length(r2g2_only)
    r2g_shared_n <- length(r2g_shared)
    # put in list
    result_per_ct[[ct]] <- data.frame(
      'cell_type' = c(ct, ct, ct), 
      'shared' = c(set_name1, set_name2, 'shared'), 
      'n' = c(r2g1_only_n, r2g2_only_n, r2g_shared_n)
    )
    # put also in the total table
    r2g_all[['d1']][[ct]] <- r2g1
    r2g_all[['d2']][[ct]] <- r2g2
  }
  # merge all the d1 r2g
  r2g_all_d1 <- unique(do.call('c', r2g_all[['d1']]))
  r2g_all_d2 <- unique(do.call('c', r2g_all[['d2']]))
  # also get the sharedness here
  r2g1_only <- setdiff(r2g_all_d1, r2g_all_d2)
  r2g2_only <- setdiff(r2g_all_d2, r2g_all_d1)
  r2g_shared <- intersect(r2g_all_d1, r2g_all_d2)
  # get the sizes
  r2g1_only_n <- length(r2g1_only)
  r2g2_only_n <- length(r2g2_only)
  r2g_shared_n <- length(r2g_shared)
  # put in list
  result_per_ct[['all']] <- data.frame(
    'cell_type' = c('all', 'all', 'all'), 
    'shared' = c(set_name1, set_name2, 'shared'), 
    'n' = c(r2g1_only_n, r2g2_only_n, r2g_shared_n)
  )
  # merge all of them
  results_all <- do.call('rbind', result_per_ct)
  return(results_all)
}

# plot the concordance
plot_concordanace <- function(qtl_effect_table, ca_effect_column='ca_effect', e_effect_column='e_effect', main='Effect sizes of caQTLs versus eQTLs', xlab='caQTL effect size', ylab='eQTL effect size') {
  # set standardized columns
  qtl_effect_table[['ca_effect']] <- qtl_effect_table[[ca_effect_column]]
  qtl_effect_table[['e_effect']] <- qtl_effect_table[[e_effect_column]]
  
  # get the concordance
  con_qtl_effect_table <- sum(sign(qtl_effect_table[['e_effect']]) == sign(qtl_effect_table[['ca_effect']])) / nrow(qtl_effect_table)
  
  # get the min and max values on the axis
  max_beta1 <- max(abs(qtl_effect_table[['ca_effect']]))
  min_beta1 <- min(abs(qtl_effect_table[['ca_effect']]))
  max_beta2 <- max(abs(qtl_effect_table[['e_effect']]))
  min_beta2 <- min(abs(qtl_effect_table[['e_effect']]))
  
  # start making the plot
  plot(x = qtl_effect_table[['ca_effect']], y = qtl_effect_table[['e_effect']], 
       ylim = c(max_beta2 * -1, max_beta2),
       xlim = c(max_beta1 * -1, max_beta1),
       xlab = xlab,
       ylab = ylab,
       pch = 16,
       cex = 0.5,
       main = main
  )
  # the top to bottom rectangle of non-significant effects
  rect(xleft = -1 * min_beta1, xright = min_beta1, ybottom = -1 * max_beta2, ytop = max_beta2, border = NA, col = rgb(red = 1, green = 1, blue = 1, alpha = 0.5))
  # the left to right rectangle of non-significant effects
  rect(xleft = -1 * max_beta1, xright = max_beta1, ybottom = -1 * min_beta2, ytop = min_beta2, border = NA, col = rgb(red = 1, green = 1, blue = 1, alpha = 0.5))
  # bottom left concordant
  rect(xleft = -1 * max_beta1, xright = -1 *min_beta1, ybottom = -1 * max_beta2, ytop = -1 * min_beta2, border = NA, col = rgb(red = 0, green = .45, blue = .7, alpha = 0.5))
  # bottom right disconcordant
  rect(xleft = min_beta1, xright = max_beta1, ybottom = -1 * max_beta2, ytop = -1 * min_beta2, border = NA, col = rgb(red = .8, green = .4, blue = 0, alpha = 0.5))
  # top left disconcordant
  rect(xleft = -1 * max_beta1, xright = -1 *min_beta1, ybottom = min_beta2, ytop = max_beta2, border = NA, col = rgb(red = .8, green = .4, blue = 0, alpha = 0.5))
  # top right disconcordant
  rect(xleft = min_beta1, xright = max_beta1, ybottom = min_beta2, ytop = max_beta2, border = NA, col = rgb(red = 0, green = .45, blue = .7, alpha = 0.5))
  # box to put concordance label in
  # rect(xleft = max_beta2 * 0.60, xright = max_beta2 * 0.95, ybottom = max_beta2 * -0.8, ytop = max_beta2 * -0.70, col = 'white')
  rect(xleft = max_beta2 * 0.50, xright = max_beta2 * 1.15, ybottom = max_beta2 * -0.9, ytop = max_beta2 * -0.60, col = 'white')
  # add concordance label
  text(x = max_beta1 * 0.75 , y = max_beta2 * -0.75, labels = c(paste('concordance', round(con_qtl_effect_table, digits = 2), sep = ':\n')))
  # add the annotation of what if concordant and disconcordant
  #text(x = max_sig_z_venema * 0.75 , y = max_sig_z_eqtlgen * 0.75, labels = c('concordant'), col = rgb(red = 0, green = .45, blue = .7))
  #text(x = max_sig_z_venema * - 0.75 , y = max_sig_z_eqtlgen * 0.75, labels = c('disconcordant'), col = rgb(red = .8, green = .4, blue = 0))
  # add dashed lines
  lines(c(max_beta1 * -1, max_beta1), c(min_beta2 * -1, min_beta2 * -1), type = "l", lty = 2)
  lines(c(max_beta1 * -1, max_beta1), c(min_beta2, min_beta2), type = "l", lty = 2)
  lines(c(min_beta1 * -1, min_beta1 * -1), c(max_beta2 * -1, max_beta2), type = "l", lty = 2)
  lines(c(min_beta1, min_beta1), c(max_beta2 * -1, max_beta2), type = "l", lty = 2)
  # save the plot
  # plots/mo_top_eqtl_over_ct_to_caqtl_betas.pdf
  # 5,5
}


####################
# Settings         #
####################

# luck seed
set.seed(7777)


####################
# Main code        #
####################

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])
# # read the location of the UCSC annotations
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/genome_annotation.tsv'
gene_anno <- fread(gene_anno_loc, header = T, sep = '\t')
# rename columns to be the same as in limix
colnames(gene_anno) <-c('chrom', 'start', 'end', 'strand', 'gs','Transcription_Start_Site','Transcript_type')
# get extra annotations for the pseudobulk output
strand_information_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eQTA/LimixExpAnnotationFile.incStrand.txt'
strand_information <- fread(strand_information_loc, header = T, sep = '\t')


# location of the hybrid method
hybrid_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/'
# read hybrid method
hybrid_output_list <- read_pseudobulk_cre_output_per_celltype(hybrid_output_loc, add_mtc = F, filter_alpha = F, add_global_nominal_threshold = F, add_local_nominal_threshold = F, filename_output = 'qtl_results_annotated_all.txt', alpha_min = .8, alpha_max = 1.2, filter_significance = F)
# get the unique mappings
hybrid_mappings <- names(hybrid_output_list)
# extract the pseudobulk ones
pseudobulk_mappings <- hybrid_mappings[grep('_pb$', hybrid_mappings)]
# extract those
pseudobulk_output_list <- hybrid_output_list[pseudobulk_mappings]
# and remove the append of '_pb'
names(pseudobulk_output_list) <- gsub('_pb', '', names(pseudobulk_output_list))
# split those
hybrid_output_list <- hybrid_output_list[setdiff(hybrid_mappings, pseudobulk_mappings)]
# merge them
hybrid_output <- do.call('rbind', hybrid_output_list)
# add z score
hybrid_output[['zscore']] <- hybrid_output[['beta']] / hybrid_output[['beta_se']]
# add a correlation based on the Z score, taking sample size and removing 2 + 10 PCs to get the degrees of freedom
hybrid_output[['r']] <- hybrid_output[['zscore']] / sqrt(hybrid_output[['zscore']]^2 + (hybrid_output[['n_samples']][1] - 12))
# add a p based z
hybrid_output[['z_from_p']] <- qnorm(hybrid_output[['p_value']] / 2) * -1 * sign(hybrid_output[['beta']])

# merge pseudobulk results
pseudobulk_output <- do.call('rbind', pseudobulk_output_list)
# rename celltype
pseudobulk_output[['cell_type']] <- gsub('_pb$', '', pseudobulk_output[['cell_type']])
# add the location to the pseudobulk info
pseudobulk_output <- cbind(pseudobulk_output, cpeaks_anno[match(pseudobulk_output[['snp_id']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# get the distances
pseudobulk_distances <- get_closest_flanks(pseudobulk_output, 'start_hg38', 'end_hg38', 'feature_start', 'feature_end')
# add those distances
pseudobulk_output[['distance']] <- pseudobulk_distances[['min_dist']]
# and category
pseudobulk_output[['category']] <- 'pseudobulk'
# add to the pseudobulk
pseudobulk_output[['strand']] <- strand_information[match(pseudobulk_output[['feature_id']], strand_information[['feature_id']]), ][['strand']]
# add z score
pseudobulk_output[['zscore']] <- pseudobulk_output[['beta']] / pseudobulk_output[['beta_se']]


# and to the hybrid method
hybrid_output <- cbind(hybrid_output, cpeaks_anno[match(hybrid_output[['snp_id']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
hybrid_distances <- get_closest_flanks(hybrid_output, 'start_hg38', 'end_hg38', 'feature_start', 'feature_end')
hybrid_output[['distance']] <- hybrid_distances[['min_dist']]
hybrid_output[['category']] <- 'hybrid'
hybrid_output[['strand']] <- strand_information[match(hybrid_output[['feature_id']], strand_information[['feature_id']]), ][['strand']]

# add significance cutoff
hybrid_output <- add_nominal_p_value_cutoff(hybrid_output, emperical_p_column='global_significance', emperical_p_cutoff=T)
pseudobulk_output <- add_nominal_p_value_cutoff(pseudobulk_output, emperical_p_column='global_significance', emperical_p_cutoff=T)

# get this as numbers
hybrid_output_n_cre_gene <- data.frame(table(unique(hybrid_output[hybrid_output[['significant_nominal_cutoff']], c('snp_id', 'feature_id', 'cell_type')])[['cell_type']]))
colnames(hybrid_output_n_cre_gene) <- c('cell_type', 'n_cre_gene')
hybrid_output_n_cre <- data.frame(table(unique(hybrid_output[hybrid_output[['significant_nominal_cutoff']], c('snp_id', 'cell_type')])[['cell_type']]))
colnames(hybrid_output_n_cre) <- c('cell_type', 'n_cre')
hybrid_output_n_gene <- data.frame(table(unique(hybrid_output[hybrid_output[['significant_nominal_cutoff']], c('feature_id', 'cell_type')])[['cell_type']]))
colnames(hybrid_output_n_gene) <- c('cell_type', 'n_gene')
# now plot them
p_n_sccre_cre_gene <- egene_numbers_to_plot(hybrid_output_n_cre_gene, number_column = 'n_cre_gene', split_long_labels_to_lines = F, celltype_column = 'cell_type', use_label_dict = F, use_color_dict = T, legendless = T)
p_n_sccre_cre <- egene_numbers_to_plot(hybrid_output_n_cre, number_column = 'n_cre', split_long_labels_to_lines = F, celltype_column = 'cell_type', use_label_dict = F, use_color_dict = T, legendless = T)
p_n_sccre_gene <- egene_numbers_to_plot(hybrid_output_n_gene, number_column = 'n_gene', split_long_labels_to_lines = F, celltype_column = 'cell_type', use_label_dict = F, use_color_dict = T, legendless = T)

# get this as numbers
pseudo_output_n_cre_gene <- data.frame(table(unique(pseudobulk_output[pseudobulk_output[['significant_nominal_cutoff']], c('snp_id', 'feature_id', 'cell_type')])[['cell_type']]))
colnames(pseudo_output_n_cre_gene) <- c('cell_type', 'n_cre_gene')
pseudo_output_n_cre <- data.frame(table(unique(pseudobulk_output[pseudobulk_output[['significant_nominal_cutoff']], c('snp_id', 'cell_type')])[['cell_type']]))
colnames(pseudo_output_n_cre) <- c('cell_type', 'n_cre')
pseudo_output_n_gene <- data.frame(table(unique(pseudobulk_output[pseudobulk_output[['significant_nominal_cutoff']], c('feature_id', 'cell_type')])[['cell_type']]))
colnames(pseudo_output_n_gene) <- c('cell_type', 'n_gene')
# now plot them
p_n_pscre_cre_gene <- egene_numbers_to_plot(hybrid_output_n_cre_gene, number_column = 'n_cre_gene', split_long_labels_to_lines = F, celltype_column = 'cell_type', use_label_dict = F, use_color_dict = T, legendless = T)
p_n_pscre_cre <- egene_numbers_to_plot(hybrid_output_n_cre, number_column = 'n_cre', split_long_labels_to_lines = F, celltype_column = 'cell_type', use_label_dict = F, use_color_dict = T, legendless = T)
p_n_pscre_gene <- egene_numbers_to_plot(hybrid_output_n_gene, number_column = 'n_gene', split_long_labels_to_lines = F, celltype_column = 'cell_type', use_label_dict = F, use_color_dict = T, legendless = T)



# filter the outputs
hybrid_output_sig <- hybrid_output[hybrid_output[['significant_nominal_cutoff']] & hybrid_output[['global_significance']], ]
pseudobulk_output_sig <- pseudobulk_output[pseudobulk_output[['significant_nominal_cutoff']] & pseudobulk_output[['global_significance']], ]
# set the locations to share these
sccre_output_sig_loc <- '/groups/umcg-franke-scrna/tmp02/users/umcg-roelen/singularity/rstudio-server/simulated_home/tables/mo_sccre_significant.tsv.gz'
pscre_output_sig_loc <- '/groups/umcg-franke-scrna/tmp02/users/umcg-roelen/singularity/rstudio-server/simulated_home/tables/mo_pscre_significant.tsv.gz'
# and write these
write.table(hybrid_output_sig, gzfile(sccre_output_sig_loc), row.names = F, col.names = T, sep = '\t')
write.table(pseudobulk_output_sig, gzfile(pscre_output_sig_loc), row.names = F, col.names = T, sep = '\t')

# merge the sccre to pscre
ps_vs_sc_cre_rep <- match_qtl_tables(hybrid_output, pseudobulk_output, allele_column = NULL, emperical_p_column = 'global_significance', significance_cutoff_emperical = T)

# we'll store stats per cell type first
stats_per_ct_l <- list()
# check each cell type
for (ct in names(ps_vs_sc_cre_rep)) {
  # extract the table
  sc_vs_ps_ct <- ps_vs_sc_cre_rep[[ct]]
  # qvalue correct
  # sc_vs_ps_ct[['sccre_empirical_feature_q_value']] <- qvalue(sc_vs_ps_ct[['sccre_empirical_feature_p_value']], lambda = 0)$qvalues
  # # subset to what is significant in eQTLgen
  # sc_vs_ps_ct <- sc_vs_ps_ct[sc_vs_ps_ct[['sccre_empirical_feature_q_value']] < 0.05, ]
  sc_vs_ps_ct <- sc_vs_ps_ct[sc_vs_ps_ct[['sccre_global_significance']], ]
  # we only have significant effect now
  sc_vs_ps_ct[['sccre_significant']] <- T
  # we can only do this if there are significant effects
  if (nrow(sc_vs_ps_ct) > 0) {
    # then do FDR correction on mo
    sc_vs_ps_ct[['pscre_fdr']] <- p.adjust(sc_vs_ps_ct[['pscre_p_value']], method = 'BH')
    # and update the significance
    sc_vs_ps_ct[['pscre_significant']] <- sc_vs_ps_ct[['pscre_fdr']] < 0.05
    # put that into a list
    sc_vs_ps_ct_list <- list()
    sc_vs_ps_ct_list[[paste('sccre','pscre', sep = '_vs_')]] <- sc_vs_ps_ct
    # get the ct singular
    ct_singulary <- strsplit(ct, '_vs_')[[1]][[1]]
    # check replication in caQTLs across cell types
    sc_vs_ps_ct_stats <- get_stats_qtl_matches(sc_vs_ps_ct_list, 
                                               significance_column_replicating_nominal='p_value', 
                                               significance_column_replicating_in_nominal='p_value', 
                                               # significance_column_replicating_emperical='empirical_feature_p_value', 
                                               significance_column_replicating_emperical='p_value', # dummy column, because all are significant
                                               se_column_replicating='beta_se', 
                                               se_column_replicating_in='beta_se'
    )
    # now replace the name
    sc_vs_ps_ct_stats[['replicating_in']] <- ct_singulary
    # put in the list
    stats_per_ct_l[[ct]] <- sc_vs_ps_ct_stats
  } else {
    print(paste('skipping due to no overlap:', ct))
  }
}
# merge all now
stats_per_ct <- do.call('rbind', stats_per_ct_l)
# rename replicating, not because it is what we are replicating, but because it is what we are replicating in
stats_per_ct[['replicating']] <- 'pseudobulk'

# check replication in caQTLs across cell types
sccre_to_pscre_celltype_stat_p <- plot_replication_stats(stats_per_ct, mid_value = 0.5, comparing_value_to_text=list('n_sig_both' = 'N sig: ', 'n_sig_tested' = 'N test: ', 'rb_emp_all' = 'Rb: '), colour_column = 'rb_emp_all', category_plot_order = c('pseudobulk', 'B', 'CD4+ T', 'CD8+ T', 'NK', 'DC', 'Monocyte', 'monocyte'), text_size = 4, legendless = T) +
  xlab('Dataset') +
  ylab('Cell type replicating in') +
  ggtitle('Replication of CREs') +
  labs(fill="Rb") +
  theme(
    # X label font size
    axis.title.x = element_text(size = 20), 
    # Y label font size
    axis.title.y = element_text(size = 20), 
    # X tick font size
    axis.text.x = element_text(size = 16), 
    # Y tick font size
    axis.text.y = element_text(size = 16), 
    # title size
    title = element_text(size = 24)
  ) + 
  scale_y_discrete(position = "right")
# save the figure
ggsave('~/plots/mo_sccre_vs_pscre_celltype_replication_simple.pdf', plot = sccre_to_pscre_celltype_stat_p, width = 3, height = 8)
# display the figure
sccre_to_pscre_celltype_stat_p
# do concordance as well
sccre_to_pscre_celltype_stat_p_cond <- plot_replication_stats(stats_per_ct, mid_value = 0.5, comparing_value_to_text=list('con_emp_emp' = ''), colour_column = 'con_emp_emp', category_plot_order = c('pseudobulk', 'B', 'CD4+ T', 'CD8+ T', 'NK', 'DC', 'Monocyte', 'monocyte'), text_size = 4, legendless = T) +
  xlab('Dataset') +
  ylab('Cell type replicating in') +
  ggtitle('Replication of CREs') +
  labs(fill="Concordance") +
  theme(
    # X label font size
    axis.title.x = element_text(size = 20), 
    # Y label font size
    axis.title.y = element_text(size = 20), 
    # X tick font size
    axis.text.x = element_text(size = 16), 
    # Y tick font size
    axis.text.y = element_text(size = 16), 
    # title size
    title = element_text(size = 24)
  ) + 
  scale_y_discrete(position = "right")
# save the figure
ggsave('~/plots/mo_sccre_vs_pscre_celltype_replication_concordance.pdf', plot = sccre_to_pscre_celltype_stat_p_cond, width = 3, height = 8)
# display the figure
sccre_to_pscre_celltype_stat_p_cond

# get overlapping numbers
cre_overlapping <- get_overlapping_peakgene_links(hybrid_output, pseudobulk_output)
# replace cell type names
cre_overlapping[['cell_type_nice']] <- rename_labels(cre_overlapping[['cell_type']])
# add cell type with sharedness
cre_overlapping[['ct_sharedness']] <- paste(cre_overlapping[['cell_type_nice']], cre_overlapping[['shared']])
# and plot them
p_overlapping_r2g_ps_sc <- ggplot(data = cre_overlapping, mapping = aes(x = cell_type_nice, y = n, fill = ct_sharedness)) +
  geom_bar(stat = 'identity', position = 'stack') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  theme(legend.position = 'none') +
  xlab('Cell type') +
  ylab('N region-gene links') + 
  ggtitle('Number of region-gene links') +
  scale_fill_manual(values = get_color_coding_dict())
# show overlap
p_overlapping_r2g_ps_sc
# save the figure
ggsave('~/plots/mo_sccre_vs_pscre_celltype_overlap.pdf', plot = p_overlapping_r2g_ps_sc, width = 4, height = 4)

# sort by effect size
hybrid_output <- hybrid_output[order(abs(hybrid_output[['zscore']])), ]
pseudobulk_output <- pseudobulk_output[order(abs(pseudobulk_output[['zscore']])), ]
# merge significant results
overlap_top_effects <- merge(pseudobulk_output[!duplicated(paste(pseudobulk_output[['snp_id']], pseudobulk_output[['feature_id']])), ], hybrid_output[!duplicated(paste(hybrid_output[['snp_id']], hybrid_output[['feature_id']])), ], by = c('snp_id', 'feature_id'))
# plot this
plot_concordanace(overlap_top_effects, ca_effect_column='zscore.x', e_effect_column='zscore.y', main='Effect sizes of pbCREs vs scCREs', xlab = 'pbCRE Z-score', ylab = 'scCRE Z-score')

