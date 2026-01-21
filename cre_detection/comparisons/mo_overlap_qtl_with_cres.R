#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_overlap_qtl_with_cres.R
# Function: overlap eQTLs and caQTLs with CREs detected in SCENIX+ or LIMIX
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(IRanges)
library(ggplot2)
library(cowplot)


####################
# Functions        #
####################


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


#' Merge QTL Data with Openness Regions Based on Genomic Overlaps
#'
#' This function merges QTL (Quantitative Trait Loci) data with genomic openness data
#' by identifying overlaps between QTL variant positions and openness regions on the same chromosome.
#' It returns the QTL data annotated with overlapping openness region identifiers and openness values.
#'
#' @param qtl_input A `data.table` or `data.frame` containing QTL data.
#' @param openness_input A `data.table` or `data.frame` containing openness region data.
#' @param variant_column_qtls Name of the column in `qtl_input` that contains variant identifiers. Default is `'snp_id'`.
#' @param chromosome_column_qtls Name of the column in `qtl_input` that contains chromosome identifiers. Default is `'snp_chromosome'`.
#' @param position_column_qtls Name of the column in `qtl_input` that contains variant positions. Default is `'snp_position'`.
#' @param openness_column_region Name of the column in `openness_input` that contains openness region identifiers. Default is `'name'`.
#' @param opennes_column_openness Name of the column in `openness_input` that contains openness values. Default is `'pct_exp'`.
#' @param openness_column_chromosome Name of the column in `openness_input` that contains chromosome identifiers. Default is `'#chrom'`.
#' @param openness_column_start Name of the column in `openness_input` that contains region start positions. Default is `'start'`.
#' @param openness_column_end Name of the column in `openness_input` that contains region end positions. Default is `'end'`.
#' @param overlapping_region_column Name of the output column for overlapping openness region identifiers. Default is `'openness_region'`.
#' @param overlapping_openness_column Name of the output column for overlapping openness values. Default is `'openness_openness'`.
#'
#' @return A `data.table` containing the original QTL data with additional columns for overlapping openness region and openness value.
#'
#' @import data.table
#' @importFrom IRanges IRanges findOverlaps pintersect
#' @export
#'
#' @examples
#' # Example usage:
#' # result <- qtl_merge_with_openness(qtl_data, openness_data)
qtl_merge_with_openness <- function(qtl_input, openness_input, variant_column_qtls='snp_id', chromosome_column_qtls='snp_chromosome', position_column_qtls='snp_position', openness_column_region='name', opennes_column_openness='pct_exp', openness_column_chromosome='#chrom', openness_column_start='start', openness_column_end='end', overlapping_region_column='openness_region', overlapping_openness_column='openness_openness', variant_window_left=0, variant_window_right=0) {
  # get the chromosomes in the qtl data
  qtl_chroms <- unique(qtl_input[[chromosome_column_qtls]])
  # and in the dars
  openness_chroms <- unique(openness_input[[openness_column_chromosome]])
  # only do the ones present in both
  chroms_both <- intersect(qtl_chroms, openness_chroms)
  # we'll save the results in a list
  overlaps_per_chrom <- list()
  # and check each chromosome
  for (chrom in chroms_both) {
    # subset to this chrom
    openness_input_regions_chromosome <- openness_input[!is.na(openness_input[[openness_column_chromosome]]) & openness_input[[openness_column_chromosome]] == chrom, ]
    qtl_regions_chromosome <- qtl_input[!is.na(qtl_input[[chromosome_column_qtls]]) & qtl_input[[chromosome_column_qtls]] == chrom, ]
    # then subset to unique variants for the QTLs
    qtl_regions_chromosome_variants <- unique(qtl_regions_chromosome[, c(variant_column_qtls, chromosome_column_qtls, position_column_qtls)])
    # turn into iranges objects
    openness_input_chromosome_iranges <- IRanges(start = openness_input_regions_chromosome[[openness_column_start]], end = openness_input_regions_chromosome[[openness_column_end]])
    # with the window supplied (making sure that we don't get a negative number)
    qtl_region_chromosome_variants_windowed_starts <- qtl_regions_chromosome_variants[[position_column_qtls]] - variant_window_left
    qtl_region_chromosome_variants_windowed_starts[qtl_region_chromosome_variants_windowed_starts < 1] <- 1
    qtl_chromosome_iranges <- IRanges(start = qtl_region_chromosome_variants_windowed_starts, end = qtl_regions_chromosome_variants[[position_column_qtls]] + variant_window_right)
    # find overlaps
    feature_chromosome_overlaps <- findOverlaps(openness_input_chromosome_iranges, qtl_chromosome_iranges)
    # extract overlapping ranges
    overlapping_ranges <- pintersect(qtl_chromosome_iranges[subjectHits(feature_chromosome_overlaps)], openness_input_chromosome_iranges[queryHits(feature_chromosome_overlaps)])
    # create a  table for the overlaps
    overlaps_table <- data.table(
      'variant_id' = qtl_regions_chromosome_variants[[variant_column_qtls]][subjectHits(feature_chromosome_overlaps)],
      'openness_region' = openness_input_regions_chromosome[[openness_column_region]][queryHits(feature_chromosome_overlaps)],
      'openness_openness' = openness_input_regions_chromosome[[opennes_column_openness]][queryHits(feature_chromosome_overlaps)]
    )
    # name them as we wanted
    colnames(overlaps_table) <- c('variant_id', overlapping_region_column, overlapping_openness_column)
    # add these positions
    qtl_regions_chromosome <- merge(qtl_regions_chromosome, overlaps_table, by.x = variant_column_qtls, by.y = 'variant_id', all.x = T, allow.cartesian=TRUE)
    # put in list
    overlaps_per_chrom[[as.character(chrom)]] <- qtl_regions_chromosome
  }
  overlaps_all <- do.call('rbind', overlaps_per_chrom)
  return(overlaps_all)
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
    # or when doing matching
    color_coding_dict[[paste(cell_type, 'matched')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'none')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, 'unmatched')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "black"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, 'matched only')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(50)[pct_whitening]
    # or when checking for matching CRE
    color_coding_dict[[paste(cell_type, 'is CRE')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'not CRE')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    # or if in CRE
    color_coding_dict[[paste(cell_type, 'in CRE')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'not in CRE')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    # if CRE in any tool
    color_coding_dict[[paste(cell_type, 'both')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'LIMIX')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, 'SCENIC+')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "black"))(100)[pct_whitening]
  }
  # general
  color_coding_dict[['AI']] <- 'darkblue'
  color_coding_dict[['NI']] <- 'darkred'
  color_coding_dict[['Actively Inflamed']] <- 'darkblue'
  color_coding_dict[['Non-Inflamed']] <- 'darkred'
  return(color_coding_dict)
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
  label_dict[['Plasma_cells']] <- 'Plasma cells'
  label_dict[['Stem_cells']] <- 'Stem cells'
  label_dict[['Stromal_cells']] <- 'Stromal cells'
  label_dict[['T_others']] <- 'other T cells'
  label_dict[['Transit_amplifying_cells']] <- 'Transit amplifying cells'
  label_dict[['AI']] <- 'Actively Inflamed'
  label_dict[['NI']] <- 'Non-Inflamed'
  return(label_dict)
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


####################
# Main code        #
####################

# location of the eQTL output
qtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# get all the QTL output
qtl_output <- get_qtls_per_celltype_limix(qtl_output_loc)
# add the top effect information
qtl_output <- add_top_effect_annotation(qtl_output)
# merge all the results
qtl_output_all <- do.call('rbind', qtl_output)


# # location of openness files
# openness_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/'
# # prepend and append
# openness_prepend <- 'mo_peaks_lane1to80_'
# openness_append <- '.bed'
# # and the openness cell types
# openness_cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
# 
# # we'll store openness data per cell type
# openness_table_per_celltype <- list()
# # read each of the openness tables
# for (cell_type in openness_cell_types) {
#   # we'll paste the path together
#   cell_type_openness_loc <- paste0(openness_output_loc, '/', openness_prepend, cell_type, openness_append)
#   # let the user know we are reading this data
#   message(paste('reading openness file at', cell_type_openness_loc))
#   # read the file
#   openness_table_per_celltype[[cell_type]] <- fread(cell_type_openness_loc, header = T, sep = '\t')
# }

# add 'chr' to the chromosome
qtl_output_all[['snp_chromosome']] <- paste0('chr', qtl_output_all[['snp_chromosome']])

# get the cpeaks overlaps for each variant
qtl_variants_all_cpeaks_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
qtl_variants_all_cpeaks <- fread(qtl_variants_all_cpeaks_loc, header = T, sep = '\t')

# # check each cell type, and add the openness for that cell type
# for (cell_type in names(openness_table_per_celltype)) {
#   # get openness for that cell type
#   openness_celltype <- openness_table_per_celltype[[cell_type]]
#   # keep only the things we need
#   openness_celltype <- openness_celltype[, c('name', 'pct_exp')]
#   # rename the column
#   colnames(openness_celltype) <- c('overlapping_feature', paste('openness', cell_type, sep = '_'))
#   # merge that onto the variant table
#   qtl_variants_all_ct_openness <- merge(qtl_variants_all_cpeaks, openness_celltype, by ='overlapping_feature', all.x = T)
#   # now sort this by openness
#   qtl_variants_all_ct_openness <- qtl_variants_all_ct_openness[order(qtl_variants_all_ct_openness[[paste('openness', cell_type, sep = '_')]], decreasing = T)]
#   # and keep the largest openness per variant
#   qtl_variants_all_ct_openness <- qtl_variants_all_ct_openness[!duplicated(qtl_variants_all_ct_openness[['snp_id']]), ]
#   # finally, add this information to the original QTL table
#   qtl_output_all[[paste('openness', cell_type, sep = '_')]] <- qtl_variants_all_ct_openness[match(qtl_output_all[['snp_id']], qtl_variants_all_ct_openness[['snp_id']]), ][[paste('openness', cell_type, sep = '_')]]
# }

# add overlapping feature to QTL
qtl_output_all[['region']] <- qtl_variants_all_cpeaks[match(qtl_output_all[['snp_id']], qtl_variants_all_cpeaks[['snp_id']]), ][['overlapping_feature']]
# filter on signifiacnce
qtl_output_all_sig <- qtl_output_all[qtl_output_all[['feature_q_value']] < 0.05 &
                                      qtl_output_all[['p_value']] < qtl_output_all[['pval_nominal_threshold_global']], ]

# get the region-gene pairs in the eqtls
r2g_eqtls <- unique(paste(qtl_output_all_sig$region, qtl_output_all_sig$feature_id))
# save these somewhere
r2g_eqtls_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/mo_variant_gene_region_significant.tsv.gz'
write.table(unique(qtl_output_all_sig[, c('snp_id', 'feature_id', 'region', 'cell_type', 'beta', 'beta_se')]), gzfile(r2g_eqtls_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# make checksum
mdfiver::create_sha256_for_file(r2g_eqtls_loc)

# colnames(ucsc_anno) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts')
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/genome_annotation.tsv'
gene_anno <- fread(gene_anno_loc, header = T, sep = '\t')
# rename columns to be the same as in limix
colnames(gene_anno) <-c('chrom', 'start', 'end', 'strand', 'gs','Transcription_Start_Site','Transcript_type')
# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# add location for the scenic table
scenic_output <- cbind(scenic_output, cpeaks_anno[match(scenic_output[['Region']], cpeaks_anno[['scenic_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# and the locations of the genes
scenic_output <- cbind(scenic_output, gene_anno[match(scenic_output[['Gene']], gene_anno[['gs']]), c('chrom', 'start', 'end')])
# get the distances again
scenic_distances <- get_closest_flanks(scenic_output, 'start_hg38', 'end_hg38', 'start', 'end')
# add that to the original table
scenic_output[['distance']] <- scenic_distances[['min_dist']]
# remove the entries that are more likely to be false positives
scenic_output_unfiltered <- scenic_output
scenic_output <- scenic_output_unfiltered[scenic_output_unfiltered[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# keep the first eRegulon if it is in both direct and extended
scenic_output <- scenic_output[order(scenic_output[['is_extended']]), ]
# grab the unique TFs
scenic_tfs <- unique(scenic_output[, c('TF', 'eRegulon_name')])
# remove duplicates, so that would be extended when we already have the direct
scenic_tfs <- scenic_tfs[!duplicated(scenic_tfs[['TF']]), ]
# now filter on thos eRegulons
scenic_output <- scenic_output[scenic_output[['eRegulon_name']] %in% scenic_tfs[['eRegulon_name']], ]
# remove the entries that are more likely to be false positives
scenic_output_unfiltered <- scenic_output
scenic_output <- scenic_output_unfiltered[scenic_output_unfiltered[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# and region-gene overlaps
scenic_output <- scenic_output[scenic_output[['distance']] > 0, ]
scenic_output <- scenic_output[scenic_output[['distance']] <= 150000, ]
# get the region-gene pairs in the SCENIC+ output
r2g_scenic <- unique(paste(gsub(':', '-', scenic_output[['Region']]), scenic_output[['Gene']]))

# check overlap with r2g of SCENIC+
length(intersect(r2g_eqtls, r2g_scenic))
# [1] 609

# annotate QTL output with the variant region being a CRE with the right gene
qtl_output_all_sig[['matching_cre']] <- !is.na(qtl_output_all_sig[['region']]) & paste(qtl_output_all_sig[['region']], qtl_output_all_sig[['feature_id']]) %in% r2g_scenic
# get number of QTLs where there is overlap
qtl_variant_cre <- data.frame(table(qtl_output_all_sig[, c('cell_type', 'matching_cre')]))
# get the number of top QTLs where there is overlap
qtl_variant_cre_top <- data.frame(table(qtl_output_all_sig[qtl_output_all_sig[['is_top_variant']], c('cell_type', 'matching_cre')]))
# order by overlap, to see if any overlap
qtl_output_all_sig <- qtl_output_all_sig[order(qtl_output_all_sig[['matching_cre']], decreasing = T), ]
# get the number of QTLs where there is any matching CRE
qtl_variant_cre_any <- data.frame(table(qtl_output_all_sig[!duplicated(paste(qtl_output_all_sig[['cell_type']], qtl_output_all_sig[['feature_id']])), c('cell_type', 'matching_cre')]))
# make string cre
qtl_variant_cre[['CRE']] <- ifelse(as.logical(qtl_variant_cre[['matching_cre']]), 'in CRE', 'not in CRE')
qtl_variant_cre_top[['CRE']] <- ifelse(as.logical(qtl_variant_cre_top[['matching_cre']]), 'in CRE', 'not in CRE')
qtl_variant_cre_any[['CRE']] <- ifelse(as.logical(qtl_variant_cre_any[['matching_cre']]), 'in CRE', 'not in CRE')
# add the cell type to the CRE annotation
qtl_variant_cre[['CRE_celltype']] <- paste(qtl_variant_cre[['cell_type']], qtl_variant_cre[['CRE']])
qtl_variant_cre_top[['CRE_celltype']] <- paste(qtl_variant_cre_top[['cell_type']], qtl_variant_cre_top[['CRE']])
qtl_variant_cre_any[['CRE_celltype']] <- paste(qtl_variant_cre_any[['cell_type']], qtl_variant_cre_any[['CRE']])
# make the plot
ggplot(data = qtl_variant_cre, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('All eQTL variants in correct CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data = qtl_variant_cre_top, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Top eQTL variant per gene in correct CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data = qtl_variant_cre_any, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Any eQTL variant per gene in correct CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# location of the eQTL output
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'
# get all the QTL output
caqtl_output <- get_qtls_per_celltype_limix(caqtl_output_loc)
# add the top effect information
caqtl_output <- add_top_effect_annotation(caqtl_output)
# merge all the results
caqtl_output_all <- do.call('rbind', caqtl_output)
# add 'chr' to the chromosome
caqtl_output_all[['snp_chromosome']] <- paste0('chr', caqtl_output_all[['snp_chromosome']])
# add overlapping feature to QTL
caqtl_output_all[['region']] <- qtl_variants_all_cpeaks[match(caqtl_output_all[['snp_id']], qtl_variants_all_cpeaks[['snp_id']]), ][['overlapping_feature']]
# filter on significance
caqtl_output_all_sig <- caqtl_output_all[caqtl_output_all[['feature_q_value']] < 0.05 &
                                       caqtl_output_all[['p_value']] < caqtl_output_all[['pval_nominal_threshold_global']], ]
# sort by the strongest effect size
caqtl_output_all_sig <- caqtl_output_all_sig[order(abs(caqtl_output_all_sig[['beta']]), decreasing = T), ]
# annotate QTL output with the variant region being a CRE with the right gene
caqtl_output_all_sig[['is_cre']] <- !is.na(caqtl_output_all_sig[['feature_id']]) & caqtl_output_all_sig[['feature_id']] %in% gsub(':', '-', scenic_output[['Region']])
# get number of QTLs where there is overlap
caqtl_variant_cre <- data.frame(table(caqtl_output_all_sig[!duplicated(paste(caqtl_output_all_sig[['cell_type']], caqtl_output_all_sig[['feature_id']])), c('cell_type', 'is_cre')]))
# make overlap into string
caqtl_variant_cre[['CRE']] <- ifelse(as.logical(caqtl_variant_cre[['is_cre']]), 'is CRE', 'not CRE')
# add cell type
caqtl_variant_cre[['CRE_celltype']] <- paste(caqtl_variant_cre[['cell_type']], caqtl_variant_cre[['CRE']])
# make the plot
ggplot(data = caqtl_variant_cre, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of chromatin peaks') +
  labs(fill = 'Chromatin peak is CRE') + 
  ggtitle('Chromatin peaks that are CREs') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# location of the hybrid method
hybrid_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/'
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
# add location info
hybrid_output <- cbind(hybrid_output, cpeaks_anno[match(hybrid_output[['snp_id']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# and distance info
hybrid_distances <- get_closest_flanks(hybrid_output, 'start_hg38', 'end_hg38', 'feature_start', 'feature_end')
hybrid_output[['distance']] <- hybrid_distances[['min_dist']]
# get extra annotations for the pseudobulk output
strand_information_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/LimixExpAnnotationFile.incStrand.txt'
strand_information <- fread(strand_information_loc, header = T, sep = '\t')
# add strand info to the hybrid output as well
hybrid_output[['strand']] <- strand_information[match(hybrid_output[['feature_id']], strand_information[['feature_id']]), ][['strand']]
# filter the hybrid one in the same way
hybrid_output_unfiltered <- hybrid_output
# check significance threshold
hybrid_output <- hybrid_output_unfiltered[hybrid_output_unfiltered[['global_significance']] == T, ]
# and filter on distance
hybrid_output <- hybrid_output[abs(hybrid_output[['distance']]) > 0, ]
hybrid_output <- hybrid_output[abs(hybrid_output[['distance']]) <= 150000, ]
# extract the significant r2g2c
hybrid_r2g2c <- unique(paste(hybrid_output[['snp_id']], hybrid_output[['feature_id']], hybrid_output[['cell_type']]))

# annotate QTL output with the variant region being a CRE with the right gene and cell type
qtl_output_all_sig[['matching_cre_limix']] <- !is.na(qtl_output_all_sig[['region']]) & paste(qtl_output_all_sig[['region']], qtl_output_all_sig[['feature_id']], qtl_output_all_sig[['cell_type']]) %in% hybrid_r2g2c
# get number of QTLs where there is overlap
qtl_variant_cre_limix <- data.frame(table(qtl_output_all_sig[, c('cell_type', 'matching_cre_limix')]))
# get the number of top QTLs where there is overlap
qtl_variant_cre_top_limix <- data.frame(table(qtl_output_all_sig[qtl_output_all_sig[['is_top_variant']], c('cell_type', 'matching_cre_limix')]))
# order by overlap, to see if any overlap
qtl_output_all_sig <- qtl_output_all_sig[order(qtl_output_all_sig[['matching_cre_limix']], decreasing = T), ]
# get the number of QTLs where there is any matching CRE
qtl_variant_cre_any_limix <- data.frame(table(qtl_output_all_sig[!duplicated(paste(qtl_output_all_sig[['cell_type']], qtl_output_all_sig[['feature_id']])), c('cell_type', 'matching_cre_limix')]))
# make string cre
qtl_variant_cre_limix[['CRE']] <- ifelse(as.logical(qtl_variant_cre_limix[['matching_cre_limix']]), 'in CRE', 'not in CRE')
qtl_variant_cre_top_limix[['CRE']] <- ifelse(as.logical(qtl_variant_cre_top_limix[['matching_cre_limix']]), 'in CRE', 'not in CRE')
qtl_variant_cre_any_limix[['CRE']] <- ifelse(as.logical(qtl_variant_cre_any_limix[['matching_cre_limix']]), 'in CRE', 'not in CRE')
# add the cell type to the CRE annotation
qtl_variant_cre_limix[['CRE_celltype']] <- paste(qtl_variant_cre_limix[['cell_type']], qtl_variant_cre_limix[['CRE']])
qtl_variant_cre_top_limix[['CRE_celltype']] <- paste(qtl_variant_cre_top_limix[['cell_type']], qtl_variant_cre_top_limix[['CRE']])
qtl_variant_cre_any_limix[['CRE_celltype']] <- paste(qtl_variant_cre_any_limix[['cell_type']], qtl_variant_cre_any_limix[['CRE']])
# make the plot
ggplot(data = qtl_variant_cre_limix, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('All eQTL variants in correct LIMIX CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data = qtl_variant_cre_top_limix, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Top eQTL variant per gene in correct LIMIX CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data = qtl_variant_cre_any_limix, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Any eQTL variant per gene in correct LIMIX CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# annotate QTL output with the variant region being a CRE with the right gene in both methods
qtl_output_all_sig[['matching_cre_both']] <- qtl_output_all_sig[['matching_cre']] & qtl_output_all_sig[['matching_cre_limix']]
# get number of QTLs where there is overlap
qtl_variant_cre_both <- data.frame(table(qtl_output_all_sig[, c('cell_type', 'matching_cre_both')]))
# get the number of top QTLs where there is overlap
qtl_variant_cre_top_both <- data.frame(table(qtl_output_all_sig[qtl_output_all_sig[['is_top_variant']], c('cell_type', 'matching_cre_both')]))
# order by overlap, to see if any overlap
qtl_output_all_sig <- qtl_output_all_sig[order(qtl_output_all_sig[['matching_cre_both']], decreasing = T), ]
# get the number of QTLs where there is any matching CRE
qtl_variant_cre_any_both <- data.frame(table(qtl_output_all_sig[!duplicated(paste(qtl_output_all_sig[['cell_type']], qtl_output_all_sig[['feature_id']])), c('cell_type', 'matching_cre_both')]))
# make string cre
qtl_variant_cre_both[['CRE']] <- ifelse(as.logical(qtl_variant_cre_both[['matching_cre_both']]), 'in CRE', 'not in CRE')
qtl_variant_cre_top_both[['CRE']] <- ifelse(as.logical(qtl_variant_cre_top_both[['matching_cre_both']]), 'in CRE', 'not in CRE')
qtl_variant_cre_any_both[['CRE']] <- ifelse(as.logical(qtl_variant_cre_any_both[['matching_cre_both']]), 'in CRE', 'not in CRE')
# add the cell type to the CRE annotation
qtl_variant_cre_both[['CRE_celltype']] <- paste(qtl_variant_cre_both[['cell_type']], qtl_variant_cre_both[['CRE']])
qtl_variant_cre_top_both[['CRE_celltype']] <- paste(qtl_variant_cre_top_both[['cell_type']], qtl_variant_cre_top_both[['CRE']])
qtl_variant_cre_any_both[['CRE_celltype']] <- paste(qtl_variant_cre_any_both[['cell_type']], qtl_variant_cre_any_both[['CRE']])
# make the plot
ggplot(data = qtl_variant_cre_both, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('All eQTL variants in correct both CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data = qtl_variant_cre_top_both, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Top eQTL variant per gene in correct both CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(data = qtl_variant_cre_any_both, mapping = aes(x = cell_type, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Any eQTL variant per gene in correct both CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# get which are significant in both
qtl_output_all_sig_both <- qtl_output_all_sig[qtl_output_all_sig[['matching_cre_both']], ]
# only SCENIC+
qtl_output_all_sig_limix_only <- qtl_output_all_sig[!(qtl_output_all_sig[['matching_cre_both']]) & qtl_output_all_sig[['matching_cre_limix']], ]
# only LIMIX
qtl_output_all_sig_scenic_only <- qtl_output_all_sig[!(qtl_output_all_sig[['matching_cre_both']]) & qtl_output_all_sig[['matching_cre']], ]
# get number of QTLs where there is overlap
qtl_variant_cre_both_scenic_limix <- data.frame(table(qtl_output_all_sig_both[, c('cell_type')]))
qtl_variant_cre_scenic_only <- data.frame(table(qtl_output_all_sig_scenic_only[, c('cell_type')]))
qtl_variant_cre_limix_only <- data.frame(table(qtl_output_all_sig_limix_only[, c('cell_type')]))
# add annotation where this came from
qtl_variant_cre_both_scenic_limix[['CRE']] <- 'both'
qtl_variant_cre_scenic_only[['CRE']] <- 'SCENIC+'
qtl_variant_cre_limix_only[['CRE']] <- 'LIMIX'
# merge them
qtl_variant_cre_both_all <- rbind(qtl_variant_cre_both_scenic_limix, qtl_variant_cre_scenic_only, qtl_variant_cre_limix_only)
# add cell type + method
qtl_variant_cre_both_all[['CRE_celltype']] <- paste(qtl_variant_cre_both_all[['Var1']], qtl_variant_cre_both_all[['CRE']])
# make a plot
ggplot(data = qtl_variant_cre_both_all, mapping = aes(x = Var1, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Any eQTL variant per gene in correct both CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# get which are significant in both
qtl_output_all_sig_both_top <- qtl_output_all_sig[qtl_output_all_sig[['is_top_variant']] & qtl_output_all_sig[['matching_cre_both']], ]
# only SCENIC+
qtl_output_all_sig_limix_only_top <- qtl_output_all_sig[qtl_output_all_sig[['is_top_variant']] & !(qtl_output_all_sig[['matching_cre_both']]) & qtl_output_all_sig[['matching_cre_limix']], ]
# only LIMIX
qtl_output_all_sig_screen_only_top <- qtl_output_all_sig[qtl_output_all_sig[['is_top_variant']] & !(qtl_output_all_sig[['matching_cre_both']]) & qtl_output_all_sig[['matching_cre']], ]
# get number of QTLs where there is overlap
qtl_variant_cre_both_scenic_limix_top <- data.frame(table(qtl_output_all_sig_both_top[, c('cell_type')]))
qtl_variant_cre_scenic_only_top <- data.frame(table(qtl_output_all_sig_scenic_only_top[, c('cell_type')]))
qtl_variant_cre_limix_only_top <- data.frame(table(qtl_output_all_sig_limix_only_top[, c('cell_type')]))
# add annotation where this came from
qtl_variant_cre_both_scenic_limix_top[['CRE']] <- 'both'
qtl_variant_cre_scenic_only_top[['CRE']] <- 'SCENIC+'
qtl_variant_cre_limix_only_top[['CRE']] <- 'LIMIX'
# merge them
qtl_variant_cre_both_all_top <- rbind(qtl_variant_cre_both_scenic_limix_top, qtl_variant_cre_scenic_only_top, qtl_variant_cre_limix_only_top)
# add cell type + method
qtl_variant_cre_both_all_top[['CRE_celltype']] <- paste(qtl_variant_cre_both_all_top[['Var1']], qtl_variant_cre_both_all_top[['CRE']])
# make a plot
ggplot(data = qtl_variant_cre_both_all_top, mapping = aes(x = Var1, y = Freq, fill = CRE_celltype)) +
  geom_bar(position = 'stack', stat = 'identity') + 
  scale_fill_manual(values = get_color_coding_dict()) +
  xlab('Cell type') +
  ylab('Number of variants') +
  labs(fill = 'Variant in CRE') + 
  ggtitle('Top eQTL variant per gene in correct both CRE') +
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
