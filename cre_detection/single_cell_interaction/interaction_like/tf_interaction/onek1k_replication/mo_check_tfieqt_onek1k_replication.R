#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_tfieqt_onek1k_replication.R
# Function: check TFa-i-eQTLs of OneK1K versus scMO
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(ggplot2)
library(qvalue)
library(stringr)


####################
# Functions        #
####################

#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @param gene_column which column to use as the gene identifier
#' @returns a list with the egenes per cell type
#' 
get_egenes_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all.txt.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, verbose=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  egenes_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(qtl_output_loc, '/', cell_type, '/', output_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    # filter the results on significance
    if (!is.null(significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[significance_column]]) &
          cell_type_output[[significance_column]] < significance_cutoff, 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
      }
    }
    # get the unique genes in this file
    unique_genes <- unique(cell_type_output[[gene_column]])
    # add to the list
    egenes_per_celltype[[cell_type]] <- unique_genes
  }
  # turn into a dataframe
  return(egenes_per_celltype)
}


get_qtls_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column='pval_nominal_threshold_global', nominal_significance_column='p_value', verbose=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  qtls_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(qtl_output_loc, '/', cell_type, '/', output_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    # filter the results on significance
    if (!is.null(significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by significance column', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[significance_column]]) &
          cell_type_output[[significance_column]] < significance_cutoff, 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by significance column', nrow(cell_type_output)))
      }
    }
    if (!is.null(nominal_cutoff_column) & !is.null(nominal_significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by nominal cutoff', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[nominal_cutoff_column]]) & 
          !is.na(cell_type_output[[nominal_significance_column]]) &
          cell_type_output[[nominal_significance_column]] <= cell_type_output[[nominal_cutoff_column]], 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by nominal cutoff', nrow(cell_type_output)))
      }
    }
    # add the celltype as a column
    cell_type_output[['cell_type']] <- cell_type
    # add to the list
    qtls_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(qtls_per_celltype)
}


add_top_effect_annotation <- function(qtls_per_celltype, feature_column='feature_id', variant_column='snp_id', significance_column='p_value', decreasing=F) {
  # go through each of the cell types
  for (cell_type in names(qtls_per_celltype)) {
    # extract that table
    qtls_celltype <- qtls_per_celltype[[cell_type]]
    # order by the significancce
    qtls_celltype_ordered <- qtls_celltype[order(qtls_celltype[[significance_column]], decreasing = decreasing), ]
    # keep only the top effect
    qtls_celltype_top <- qtls_celltype_ordered[!duplicated(qtls_celltype_ordered[[feature_column]]), ]
    # and annotate in the original table if the variant-feature combination was the top one
    qtls_celltype[['is_top_variant']] <- paste(qtls_celltype[[variant_column]], qtls_celltype[[feature_column]]) %in% paste(qtls_celltype_top[[variant_column]], qtls_celltype_top[[feature_column]])
    # put that back in the list
    qtls_per_celltype[[cell_type]] <- qtls_celltype
  }
  return(qtls_per_celltype)
}


plot_concondance <- function(dataset_to_compare, d1_effect_column='d1_zscore', d2_effect_column='d2_zscore') {
  # remove any NA
  dataset_to_compare <- dataset_to_compare[!is.na(dataset_to_compare[[d1_effect_column]]) & !is.na(dataset_to_compare[[d2_effect_column]]), ]
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


#' get the number eGenes per cell type from QTL output
#' 
#' @param unfiltered_loc base location of the folder containing files to filter
#' @param unfiltered_file which output file to read for the results
#' @param filtered_loc base location of the folder to put filtered output
#' @param filtered_file name of filtered file to create
#' @param significance_column column denoting significance
#' @param significance_cutoff cutoff for which to set significance
#' @param verbose print progress
#' @param add_mtc add multiple testing before filtering down
#' @param mtc_column the column of values to apply multiple testing on
#' @param feature_mtc_column the column that has the feature group to perform the multiple testing on
#' @param mtc_column_to_add the name of the column that has the mtc-corrected values
#' @param folders optional vector of folders to look at specifically
#' @returns 0 if success
#' 
filter_output_by_significance <- function(unfiltered_loc, unfiltered_file='qtl_results_all.txt.gz', filtered_loc=NULL, filtered_file=NULL, significance_column='p_value', significance_cutoff=0.05, verbose=T, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', folders=NULL) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(unfiltered_loc, full.names = F, recursive = F)
  # check if overlaps with the folders that we want to take a look at
  if (!is.null(folders)) {
    cell_types <- intersect(cell_types, folders)
  }
  # we will store the results in a list for now
  numbers_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(unfiltered_loc, '/', cell_type, '/', unfiltered_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- fread(full_cell_type_path, sep = '\t', header = T)
    
    # get the features and the emperical p value
    if (add_mtc) {
      # get just the two columns we care about
      cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column)]
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
    }
    
    # print progress if requested
    if (verbose) {
      print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
    }
    # filter
    cell_type_output <- cell_type_output[
      !is.na(cell_type_output[[significance_column]]) &
        cell_type_output[[significance_column]] < significance_cutoff, 
    ]
    if (verbose) {
      print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
    }
    # get the output location
    output_dir <- unfiltered_loc
    # if supplied, set the output directory
    if (!is.null(filtered_loc)) {
      output_dir <- filtered_loc
    }
    # set the output file location
    output_file <- NULL
    # if supplied set that
    if (!is.null(filtered_file)) {
      output_file <- filtered_file
    }
    else {
      output_file <- paste('filtered', '.', unfiltered_file, sep = '')
    }
    # make the full path
    full_output_loc <- paste(output_dir, '/', cell_type, '/', output_file, sep = '')
    # store the gz connection if we need it
    full_output_loc_wzip <- full_output_loc
    # gzip it if the extention ends on gz
    if (grepl('.gz$', full_output_loc)) {
      full_output_loc_wzip <- gzfile(full_output_loc)
    }
    # write result
    write.table(cell_type_output, full_output_loc_wzip, sep = '\t', row.names= F, col.names = T)
    # create md5
    mdfiver::create_md5_for_file(full_output_loc)
  }
  # return 0 upon success
  return(0)
}


#' Merge Chromosome Output
#'
#' This function merges chromosome-specific output files from different cell types into a single file.
#'
#' @param input_dir A character string specifying the directory containing the input files.
#' @param input_prepend A character string specifying the prefix of the input files. Default is 'qtl_results_all_qval_'.
#' @param input_append A character string specifying the suffix of the input files. Default is '_fdr01_significant.txt.gz'.
#' @param output_dir A character string specifying the directory to save the merged output file. Default is NULL.
#' @param output_file A character string specifying the name of the merged output file. Default is 'qtl_results_all_qval_allchroms_fdr01_significant.txt.gz'.
#' @param folders A character vector specifying the folders to include. Default is NULL.
#'
#' @return An integer value indicating the success of the operation. The function writes the merged data to the specified output location.
#' @export
#'
#' @examples
#' \dontrun{
#'   merge_chromosome_output(
#'     input_dir = "path/to/input/dir",
#'     input_prepend = 'qtl_results_all_qval_',
#'     input_append = '_fdr01_significant.txt.gz',
#'     output_dir = "path/to/output/dir",
#'     output_file = 'qtl_results_all_qval_allchroms_fdr01_significant.txt.gz',
#'     folders = c("folder1", "folder2")
#'   )
#' }
merge_chromosome_output <- function(input_dir, input_prepend='qtl_results_all_qval_', input_append='_fdr01_significant.txt.gz', output_dir=NULL, output_file='qtl_results_all_qval_allchroms_fdr01_significant.txt.gz', folders=NULL) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(input_dir, full.names = F, recursive = F)
  # subset to folders we are interested in, if supplied with that option
  if (!is.null(folders)) {
    cell_types <- intersect(cell_types, folders)
  }
  # we will store the results in a list for now
  numbers_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the cell type folder
    input_celltype_dir <- paste(input_dir, '/', cell_type, '/', sep = '')
    # list the files
    input_files_celltype <- list.files(input_celltype_dir)
    # now filter only for the ones we want
    input_files_celltype <- input_files_celltype[grep(paste(input_prepend, '\\d+', input_append, '$', sep  = ''), input_files_celltype)]
    # we'll store each file in a list
    input_files_celltype_list <- list()
    # and go through each file
    for (input_file_celltype in input_files_celltype) {
      # read the file
      input_celltype_chrom <- fread(paste(input_celltype_dir, '/', input_file_celltype, sep = ''), header = T, sep = '\t')
      # put in the list
      if (nrow(input_celltype_chrom) > 0) {
        input_files_celltype_list[[input_file_celltype]] <- input_celltype_chrom
      }
      else {
        warning(paste('skipping', paste(input_celltype_dir, '/', input_file_celltype, sep = ''), 'because it has no rows'))
      }
    }
    # now merge all of them
    input_celltypes_all <- do.call('rbind', input_files_celltype_list)
    # get the output location
    output_location <- input_dir
    # if supplied, set the output directory
    if (!is.null(output_dir)) {
      output_location <- output_dir
    }
    # make full output location
    full_output_loc <- paste(output_location, '/', cell_type, '/', output_file, sep = '')
    # store the gz connection if we need it
    full_output_loc_wzip <- full_output_loc
    # gzip it if the extention ends on gz
    if (grepl('.gz$', full_output_loc)) {
      full_output_loc_wzip <- gzfile(full_output_loc)
    }
    # write result
    write.table(input_celltypes_all, full_output_loc_wzip, sep = '\t', row.names= F, col.names = T)
    # create md5
    mdfiver::create_md5_for_file(full_output_loc)
  }
  return(0)
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
  color_coding_dict[['okada']] <- 'lightgreen'
  color_coding_dict[['onek1k']] <- 'lightpink'
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

# set location of both tables
tfa_ieqtls_scmo_loc <- paste('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/', 'results_fdr_with_replication.tsv.gz', sep = '/')
tfa_ieqtls_onek1k_loc <- paste('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/onek1k/L1/all/merged/', 'results_fdr.tsv.gz', sep = '/')
# read these
tfa_ieqtls <- fread(tfa_ieqtls_scmo_loc, header = T, sep = '\t')
tfa_ieqtls_onek1k <- fread(tfa_ieqtls_onek1k_loc, header = T, sep = '\t')
# rename columns
colnames(tfa_ieqtls) <- gsub('tf', 'tfa', colnames(tfa_ieqtls))
colnames(tfa_ieqtls_onek1k) <- gsub('region', 'tfa', colnames(tfa_ieqtls_onek1k))
# remove the i-egene from the eregulon (added in replication step)
tfa_ieqtls_onek1k[['tfa']] <- gsub('\\)_.*', ')', tfa_ieqtls_onek1k[['tfa']])
# add explicit tf column
tfa_ieqtls[['tf']] <- gsub("_(extended|direct).*", "", tfa_ieqtls[['tfa']])
tfa_ieqtls_onek1k[['tf']] <- gsub("_(extended|direct).*", "", tfa_ieqtls_onek1k[['tfa']])
# and add a directed tfa beta
tfa_ieqtls[['tfa_beta_directed']] <- tfa_ieqtls[['tfa_beta']]
tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa_beta_directed']] <- -1 * tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa_beta']]
tfa_ieqtls[['tfa:genotype_beta_directed']] <- tfa_ieqtls[['tfa:genotype_beta']]
tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta_directed']] <- -1 * tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta']]
tfa_ieqtls_onek1k[['tfa_beta_directed']] <- tfa_ieqtls_onek1k[['tfa_beta']]
tfa_ieqtls_onek1k[tfa_ieqtls_onek1k[['tfa_gene_direction']] == '-/+', ][['tfa_beta_directed']] <- -1 * tfa_ieqtls_onek1k[tfa_ieqtls_onek1k[['tfa_gene_direction']] == '-/+', ][['tfa_beta']]
tfa_ieqtls_onek1k[['tfa:genotype_beta_directed']] <- tfa_ieqtls_onek1k[['tfa:genotype_beta']]
tfa_ieqtls_onek1k[tfa_ieqtls_onek1k[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta_directed']] <- -1 * tfa_ieqtls_onek1k[tfa_ieqtls_onek1k[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta']]
# add z
tfa_ieqtls[['genotype_z']] <- tfa_ieqtls[['genotype_beta']] / tfa_ieqtls[['genotype_se']]
tfa_ieqtls[['tfa_z']] <- tfa_ieqtls[['tfa_beta']] / tfa_ieqtls[['tfa_se']]
tfa_ieqtls[['tfa:genotype_z']] <- tfa_ieqtls[['tfa:genotype_beta']] / tfa_ieqtls[['tfa:genotype_se']]
tfa_ieqtls[['tfa:genotype_rep_z']] <- tfa_ieqtls[['tfa:genotype_beta_rep']] / tfa_ieqtls[['tfa:genotype_se_rep']]
tfa_ieqtls[['tfa:genotype_z_directed']] <- tfa_ieqtls[['tfa:genotype_beta_directed']] / tfa_ieqtls[['tfa:genotype_se']]
tfa_ieqtls[['tfa_z_directed']] <- tfa_ieqtls[['tfa_beta_directed']] / tfa_ieqtls[['tfa_se']]
tfa_ieqtls_onek1k[['genotype_z']] <- tfa_ieqtls_onek1k[['genotype_beta']] / tfa_ieqtls_onek1k[['genotype_se']]
tfa_ieqtls_onek1k[['tfa:genotype_z']] <- tfa_ieqtls_onek1k[['tfa:genotype_beta']] / tfa_ieqtls_onek1k[['tfa:genotype_se']]
tfa_ieqtls_onek1k[['tfa_z']] <- tfa_ieqtls_onek1k[['tfa_beta']] / tfa_ieqtls_onek1k[['tfa_se']]
tfa_ieqtls_onek1k[['tfa:genotype_z_directed']] <- tfa_ieqtls_onek1k[['tfa:genotype_beta_directed']] / tfa_ieqtls_onek1k[['tfa:genotype_se']]
tfa_ieqtls_onek1k[['tfa_z_directed']] <- tfa_ieqtls_onek1k[['tfa_beta_directed']] / tfa_ieqtls_onek1k[['tfa_se']]


# read the eQTL output
eqtl_output_onek1k_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/output/L1/'
eqtl_output_onek1k_per_ct <- get_qtls_per_celltype_limix(eqtl_output_onek1k_loc, output_file = 'top_qtl_results_all.txt.gz', significance_column = NULL, nominal_significance_column = NULL)
# perform qvalue on these top effects
for (cell_type in names(eqtl_output_onek1k_per_ct)) {
  eqtl_output_onek1k_per_ct[[cell_type]][['feature_q_value']] <- qvalue::qvalue(eqtl_output_onek1k_per_ct[[cell_type]][['empirical_feature_p_value']])$qvalues
}
# merge all of them
eqtl_output_onek1k_top <- do.call(rbind, eqtl_output_onek1k_per_ct)
# then only keep what is significant at the qvalue level
eqtl_output_onek1k_top <- eqtl_output_onek1k_top[eqtl_output_onek1k_top[['feature_q_value']] < 0.05, ]

# location of the eQTL output
eqtl_output_scmo_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# get all the QTL output
eqtl_output_scmo <- get_qtls_per_celltype_limix(eqtl_output_scmo_loc)
# add the top effect information
eqtl_output_scmo <- add_top_effect_annotation(eqtl_output_scmo)
# merge all the results
eqtl_output_scmo_all <- do.call('rbind', eqtl_output_scmo)
# order by eQTL effect size
eqtl_output_scmo_all[['z']] <- eqtl_output_scmo_all[['beta']] / eqtl_output_scmo_all[['beta_se']]
eqtl_output_scmo_all <- eqtl_output_scmo_all[order(abs(eqtl_output_scmo_all[['z']])), ]
# add eQTL to the tf-i-eQTL table
tfa_ieqtls[['scmo_eqtl_z']] <- eqtl_output_scmo_all[match(paste(tfa_ieqtls[['variant']], tfa_ieqtls[['gene']]), paste(eqtl_output_scmo_all[['snp_id']], eqtl_output_scmo_all[['feature_id']])), ][['z']] * -1

# plot the single-cell vs bulk eQTL z-scores
p_sc_vs_bulk_eqtl <- plot_concondance(tfa_ieqtls[!duplicated(paste(tfa_ieqtls[['variant']], tfa_ieqtls[['gene']])), ], d1_effect_column = 'genotype_z', d2_effect_column = 'scmo_eqtl_z') + 
  xlab('scMO ps-eQTL Z-score') +
  ylab('scMO sc-eQTL Z-score')
# show plot
p_sc_vs_bulk_eqtl

# get all the QTL output
eqtl_output_onek1k <- get_qtls_per_celltype_limix(eqtl_output_onek1k_loc, output_file = 'qtl_results_all_qval_allchroms_fdr005_significant.txt.gz', significance_column = NULL, nominal_significance_column = NULL)
# add the top effect information
eqtl_output_onek1k <- add_top_effect_annotation(eqtl_output_onek1k)
# for each of these, add get the p value cutoff
for (cell_type in names(eqtl_output_onek1k)) {
  eqtl_output_onek1k[[cell_type]][['p_value_cutoff_global']] <- max(eqtl_output_onek1k[[cell_type]][eqtl_output_onek1k[[cell_type]][['is_top_variant']], ][['p_value']])
}
# merge all of them
eqtl_output_onek1k_all <- rbindlist(eqtl_output_onek1k[c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')], fill = T)
# order by eQTL effect size
eqtl_output_onek1k_all[['z']] <- eqtl_output_onek1k_all[['beta']] / eqtl_output_onek1k_all[['beta_se']]
eqtl_output_onek1k_all <- eqtl_output_onek1k_all[order(abs(eqtl_output_onek1k_all[['z']])), ]
# add eQTL to the tf-i-eQTL table
tfa_ieqtls_onek1k[['onek1k_eqtl_z']] <- eqtl_output_onek1k_all[match(paste(tfa_ieqtls_onek1k[['variant']], tfa_ieqtls_onek1k[['gene']]), paste(eqtl_output_onek1k_all[['snp_id']], eqtl_output_onek1k_all[['feature_id']])), ][['z']] * -1

# plot the single-cell vs bulk eQTL z-scores
p_sc_vs_bulk_eqtl_onek1k <- plot_concondance(tfa_ieqtls_onek1k[!duplicated(paste(tfa_ieqtls_onek1k[['variant']], tfa_ieqtls_onek1k[['gene']])), ], d1_effect_column = 'genotype_z', d2_effect_column = 'onek1k_eqtl_z') + 
  xlab('OneK1K ps-eQTL Z-score') +
  ylab('OneK1K sc-eQTL Z-score')
# show plot
p_sc_vs_bulk_eqtl_onek1k

# filter the onek1k ieqtls to only those that are significant in the eQTL analysis
tfa_ieqtls_onek1k_confined <- tfa_ieqtls_onek1k[
  tfa_ieqtls_onek1k[['gene']] %in% eqtl_output_onek1k_top[['feature_id']], 
]

# redo MTC for these
p_columns <- colnames(tfa_ieqtls_onek1k_confined)[grep('_p$', colnames(tfa_ieqtls_onek1k_confined))]
# do MTC for each of these columns
for (p_column in p_columns) {
  # check which rows have a value for this p
  p_valid_i <- !is.na(tfa_ieqtls_onek1k_confined[[p_column]])
  # add MTC column
  bh_column <- gsub('_p$', '_bh', p_column)
  tfa_ieqtls_onek1k_confined[[bh_column]] <- NA
  # extract valid p values, and calculate B&H
  bhs <- p.adjust(tfa_ieqtls_onek1k_confined[[p_column]][p_valid_i], method = 'BH')
  # then place those BHs
  tfa_ieqtls_onek1k_confined[[bh_column]][p_valid_i] <- bhs
  # add bf column as well
  bf_column <- gsub('_p$', '_bf', p_column)
  # adjust the valid p values
  bfs <- p.adjust(tfa_ieqtls_onek1k_confined[[p_column]][p_valid_i], method = 'bonferroni')
  # add the bonferroni column
  tfa_ieqtls_onek1k_confined[[bf_column]] <- NA
  # and place the bonferroni adjusted p values
  tfa_ieqtls_onek1k_confined[[bf_column]][p_valid_i] <- bfs
}


# merge the tf and region qtls
tfa_tfa_ieqtls_onek1k_overlap <- merge(
  tfa_ieqtls, 
  tfa_ieqtls_onek1k_confined, 
  by.x = c('variant', 'tfa', 'gene'), 
  by.y = c('variant', 'tfa', 'gene')
)
# get what is significant in both
tfa_tfa_ieqtls_onek1k_overlap_tfasig <- tfa_tfa_ieqtls_onek1k_overlap[
    tfa_tfa_ieqtls_onek1k_overlap[['anova_bh.x']] & 
    tfa_tfa_ieqtls_onek1k_overlap[['genotype_bh.x']] & 
    tfa_tfa_ieqtls_onek1k_overlap[['tfa:genotype_bh.x']] & 
    tfa_tfa_ieqtls_onek1k_overlap[['anova_bh.y']] < 0.05 & 
    tfa_tfa_ieqtls_onek1k_overlap[['genotype_bh.y']] < 0.05 &
    tfa_tfa_ieqtls_onek1k_overlap[['tfa:genotype_bh.y']] < 0.05, 
]
# show genotype level concordance
p_gt_vs_gt_ieqtls_onek1k_all <- plot_concondance(tfa_tfa_ieqtls_onek1k_overlap_tfasig, d1_effect_column = 'genotype_z.x', d2_effect_column = 'genotype_z.y') + 
  xlab('scMO eQTL Z-score') +
  ylab('OneK1K eQTL Z-score')
# show the plot
p_gt_vs_gt_ieqtls_onek1k_all

# also a pseudobulk level
p_gt_vs_gt_ps_ieqtls_onek1k_all <- plot_concondance(tfa_tfa_ieqtls_onek1k_overlap, d1_effect_column = 'scmo_eqtl_z', d2_effect_column = 'onek1k_eqtl_z') + 
  xlab('scMO ps-eQTL Z-score') +
  ylab('OneK1K ps-eQTL Z-score')
# show the plot
p_gt_vs_gt_ps_ieqtls_onek1k_all

# show TFa level concordance
p_tfa_vs_tfa_ieqtls_onek1k_all <- plot_concondance(tfa_tfa_ieqtls_onek1k_overlap_tfasig, d1_effect_column = 'tfa_z.x', d2_effect_column = 'tfa_z.y') + 
  xlab('scMO eQTL Z-score') +
  ylab('OneK1K eQTL Z-score')
# show the plot
p_tfa_vs_tfa_ieqtls_onek1k_all

# show concordance
p_tfa_vs_tfa_ieqtls_onek1k_all <- plot_concondance(tfa_tfa_ieqtls_onek1k_overlap_tfasig, d1_effect_column = 'tfa:genotype_z_directed.x', d2_effect_column = 'tfa:genotype_z_directed.y') + 
  xlab('scMO TFa-i-eQTL interaction Z-score') +
  ylab('OneK1K TFa-i-eQTL interaction Z-score')
# show the plot
p_tfa_vs_tfa_ieqtls_onek1k_all
# save plot
ggsave(file = '~/multiome/plots/mo_tfa_vs_tfa_ieqtls_onek1k_all.pdf', plot = p_tfa_vs_tfa_ieqtls_onek1k_all, width = 5, height = 5)

# read the TFe-TFa correlations
tfe_tfa_scmo_loc <- '~/multiome/tables/mo_tfe_vs_tfa_correlation_conditions.tsv.gz'
tfe_tfa_onek1k_loc <- '~/multiome/tables/mo_tfe_vs_tfa_correlation_onek1k.tsv.gz'
tfe_tfa_okada_loc <- '~/multiome/tables/mo_tfe_vs_tfa_correlation_okada.tsv.gz'
# read scMO
tfe_tfa_scmo <- fread(tfe_tfa_scmo_loc, header = T, sep = '\t')
# read onek1k
tfe_tfa_onek1k <- fread(tfe_tfa_onek1k_loc, header = T, sep = '\t')
# add gene signature direction
tfe_tfa_onek1k[['gene_signature_direction']] <- str_extract(tfe_tfa_onek1k[['eregulon']], '\\+\\/\\+|\\-\\/\\-|\\+\\/\\-|\\-\\/\\+')
# add corrected correlation
tfe_tfa_onek1k[['correlation_corrected']] <- tfe_tfa_onek1k[['correlation']]
tfe_tfa_onek1k[tfe_tfa_onek1k[['gene_signature_direction']] %in% c('-/-', '-/+'), ][['correlation_corrected']] <- -1 * tfe_tfa_onek1k[tfe_tfa_onek1k[['gene_signature_direction']] %in% c('-/-', '-/+'), ][['correlation']]
# add scMO type column
tfe_tfa_onek1k[['eregulon_scmo']] <- gsub('_', '-', tfe_tfa_onek1k[['eregulon']])
# get okada
tfe_tfa_okada <- fread(tfe_tfa_okada_loc, header = T, sep = '\t')
# add scMO type column
tfe_tfa_okada[['eregulon_scmo']] <- gsub('_', '-', tfe_tfa_okada[['eregulon']])
# rename columns
colnames(tfe_tfa_onek1k) <- c('tf','eregulon_iegene','correlation_uncorrected','dataset','eregulon_original','eregulon_nice', 'gene_signature_direction', 'correlation', 'eregulon')
colnames(tfe_tfa_okada) <- c('tf','eregulon_iegene','correlation','dataset','eregulon_original','eregulon_nice', 'eregulon')
# split scmo
tfe_tfa_scmo_split <- rbind(
  data.frame(
    eregulon = tfe_tfa_scmo[['eregulon']], 
    correlation = tfe_tfa_scmo[['rho_ut']], 
    dataset = rep('scMO UT', times = nrow(tfe_tfa_scmo))
  ), 
  data.frame(
    eregulon = tfe_tfa_scmo[['eregulon']], 
    correlation = tfe_tfa_scmo[['rho_24hca']], 
    dataset = rep('scMO 24hCA', times = nrow(tfe_tfa_scmo))
  )
)
# merge all onto the full table
tfe_tfa_all <- do.call('rbind', 
                       list(tfe_tfa_scmo_split, 
                            tfe_tfa_onek1k[, c('eregulon', 'correlation', 'dataset')], 
                            tfe_tfa_okada[, c('eregulon', 'correlation', 'dataset')])
)

# the non TF-i-eQTLs
ggplot(
  data = tfe_tfa_all[complete.cases(tfe_tfa_all), ],
  mapping = aes(
    x=dataset, 
    y=eregulon, 
    fill=correlation)
) + geom_tile() + scale_fill_gradient2(low='darkblue', mid = 'white', high='darkred', midpoint = 0)  + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('Dataset') + ylab('eRegulon') + 
  theme(legend.position = 'none')

