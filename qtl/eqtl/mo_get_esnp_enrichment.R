#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_esnp_enrichment.R
# Function: check for overlap of eSNPs in open chromatin or SCREEN annotated cCREs
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(IRanges)
library(ggplot2)
library(cowplot)
library(progress)
library(pbapply)


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
  color_coding_dict[['CD4_T_cells']] <- '#264F88'
  color_coding_dict[['CD4T']] <- '#264F88'
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
  color_coding_dict[['CD4+ T cells']] <- '#264F88'
  color_coding_dict[['CD4+ T']] <- '#264F88'
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
  }
  # general
  color_coding_dict[['AI']] <- 'darkblue'
  color_coding_dict[['NI']] <- 'darkred'
  color_coding_dict[['Actively Inflamed']] <- 'darkblue'
  color_coding_dict[['Non-Inflamed']] <- 'darkred'
  return(color_coding_dict)
}


check_only_matching_openness <- function(qtl_output_all, cutoff=0.001) {
  # get the unique cell types
  cell_types <- unique(qtl_output_all[['cell_type']])
  # add the 'only' openness
  qtl_output_all[['only_matching']] <- apply(qtl_output_all, 1, function(x, cutoff) {
    # make sure to only keep the columns that mention the cell types and the column that has the cell type
    x <- x[grepl(paste('cell_type', paste(cell_types, collapse='|'), sep = '|'), names(x))]
    # get the openness columns
    matching_openness <- x[[paste('openness', x[['cell_type']], sep = '_')]]
    # make sure we are looking at numbers
    matching_openness <- as.numeric(matching_openness)
    # if the value is already NA or doesn't fit the cutoff, we don't have to do anything else
    if (is.na(matching_openness) | matching_openness < cutoff) {
      return(F)
    }
    # get all the other openness columns
    openness_columns <- names(x)[grepl('openness_', names(x))]
    # remove the matching one
    openness_columns <- setdiff(openness_columns, paste('openness', x[['cell_type']], sep = '_'))
    # get the openess values
    any_openness <- x[openness_columns]
    # make sure we are looking at numbers
    any_openness <- as.numeric(any_openness)
    # and then anything that is not NA
    any_openness <- any_openness[!is.na(any_openness)]
    # set the max openness
    max_any_openness <- 0
    if (length(any_openness) > 0) {
      max_any_openness <- max(any_openness)
    }
    # now check if any of the others are bigger than cutoff, and if the cell type one is bigger than the cutoff
    if (matching_openness < cutoff) {
      return(F)
    } else if (matching_openness >= cutoff & max_any_openness >= cutoff) {
      return(F)
    } else if (matching_openness >= cutoff & max_any_openness < cutoff) {
      return(T)
    }
  }, cutoff)
  return(qtl_output_all)
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


add_significant_celltypes_as_strings_vectorised <- function(overlap_table, variant1_col='hit1', variant2_col='hit2', trait1_col='trait1', trait2_col='trait2', cell_type_column='cell_type') {
  # convert to dataframe format
  overlap_table_df <- data.frame(overlap_table)
  # add overlap 
  overlap_table_df[['full_overlap']] <- paste(overlap_table_df[[variant1_col]], overlap_table_df[[variant2_col]],
                                              overlap_table_df[[trait1_col]], overlap_table_df[[trait2_col]])
  
  # split into groups
  split_cts <- split(overlap_table_df[[cell_type_column]], overlap_table_df[['full_overlap']])
  
  # summarise to string per group
  ct_strings <- sapply(split_cts, function(x) {
    ct_string <- paste(sort(unique(x)), collapse = ",")
    return(ct_string)
  })
  
  # map back
  overlap_table_df[['other_ct']] <- ct_strings[overlap_table_df[['full_overlap']]]
  # back to data table format
  return(data.table(overlap_table_df))
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


# location of openness files
openness_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/'
# prepend and append
openness_prepend <- 'mo_peaks_lane1to80_'
openness_append <- '.bed'
# and the openness cell types
openness_cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

# we'll store openness data per cell type
openness_table_per_celltype <- list()
# read each of the openness tables
for (cell_type in openness_cell_types) {
  # we'll paste the path together
  cell_type_openness_loc <- paste0(openness_output_loc, '/', openness_prepend, cell_type, openness_append)
  # let the user know we are reading this data
  message(paste('reading openness file at', cell_type_openness_loc))
  # read the file
  openness_table_per_celltype[[cell_type]] <- fread(cell_type_openness_loc, header = T, sep = '\t')
}

# add 'chr' to the chromosome
qtl_output_all[['snp_chromosome']] <- paste0('chr', qtl_output_all[['snp_chromosome']])

# get the cpeaks overlaps for each variant
qtl_variants_all_cpeaks_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
qtl_variants_all_cpeaks <- fread(qtl_variants_all_cpeaks_loc, header = T, sep = '\t')

# check each cell type, and add the openness for that cell type
for (cell_type in names(openness_table_per_celltype)) {
  # get openness for that cell type
  openness_celltype <- openness_table_per_celltype[[cell_type]]
  # keep only the things we need
  openness_celltype <- openness_celltype[, c('name', 'pct_exp')]
  # rename the column
  colnames(openness_celltype) <- c('overlapping_feature', paste('openness', cell_type, sep = '_'))
  # merge that onto the variant table
  qtl_variants_all_ct_openness <- merge(qtl_variants_all_cpeaks, openness_celltype, by ='overlapping_feature', all.x = T)
  # now sort this by openness
  qtl_variants_all_ct_openness <- qtl_variants_all_ct_openness[order(qtl_variants_all_ct_openness[[paste('openness', cell_type, sep = '_')]], decreasing = T)]
  # and keep the largest openness per variant
  qtl_variants_all_ct_openness <- qtl_variants_all_ct_openness[!duplicated(qtl_variants_all_ct_openness[['snp_id']]), ]
  # finally, add this information to the original QTL table
  qtl_output_all[[paste('openness', cell_type, sep = '_')]] <- qtl_variants_all_ct_openness[match(qtl_output_all[['snp_id']], qtl_variants_all_ct_openness[['snp_id']]), ][[paste('openness', cell_type, sep = '_')]]
}

# add the max openness of any region if present
qtl_output_all[['openness_max']] <- apply(qtl_output_all, 1, function(x) {
  # get the openness columns
  any_openness <- x[grepl('openness_', names(x))]
  # make sure we are looking at numbers
  any_openness <- as.numeric(any_openness)
  # and then anything that is not NA
  any_openness <- any_openness[!is.na(any_openness)]
  # if any values are not NA, we'll return that
  if (length(any_openness) > 0) {
    return(max(any_openness))
  }
  else {
    return(NA)
  }
})

# add the matching openness of any region if present
qtl_output_all[['openness_matching']] <- apply(qtl_output_all, 1, function(x) {
  # get the openness columns
  matching_openness <- x[[paste('openness', x[['cell_type']], sep = '_')]]
  # make sure we are looking at numbers
  matching_openness <- as.numeric(matching_openness)
  # return this
  return(matching_openness)
})

# add the only-matching openness
qtl_output_all <- check_only_matching_openness(qtl_output_all)

# subset to what is significant also at the snp-level
qtl_output_all_snpsig <- qtl_output_all[!is.na(qtl_output_all[['p_value']]) & 
                                          !is.na(qtl_output_all[['pval_nominal_threshold_global']]) &
                                          qtl_output_all[['p_value']] <= qtl_output_all[['pval_nominal_threshold_global']], ]

# add information about if the variant is in open chromatin, set to none first
qtl_output_all_snpsig[['in_open_chromatin']] <- 'none'
# update where there is any openness
qtl_output_all_snpsig[!is.na(qtl_output_all_snpsig[['openness_max']]) & qtl_output_all_snpsig[['openness_max']] >= 0.001, 'in_open_chromatin'] <- 'unmatched'
# or if it matches the cell type
qtl_output_all_snpsig[!is.na(qtl_output_all_snpsig[['openness_matching']]) & qtl_output_all_snpsig[['openness_matching']] >= 0.001, 'in_open_chromatin'] <- 'matched'
# also add information if it is only matched
qtl_output_all_snpsig[['in_open_chromatin_only']] <- qtl_output_all_snpsig[['in_open_chromatin']]
# and add that as added information
qtl_output_all_snpsig[qtl_output_all_snpsig[['only_matching']], 'in_open_chromatin_only'] <- 'matched only'
# set levels
qtl_output_all_snpsig[['in_open_chromatin']] <- factor(qtl_output_all_snpsig[['in_open_chromatin']], levels = c('none', 'unmatched', 'matched'))
qtl_output_all_snpsig[['in_open_chromatin_only']] <- factor(qtl_output_all_snpsig[['in_open_chromatin_only']], levels = c('none', 'unmatched', 'matched', 'matched only'))
# now make this into a table per cell type
qtl_output_all_snpsig_openatac_occurences <- data.frame(table(qtl_output_all_snpsig[, c('cell_type', 'in_open_chromatin')]))
qtl_output_all_snpsig_openatac_occurences_only <- data.frame(table(qtl_output_all_snpsig[, c('cell_type', 'in_open_chromatin_only')]))
# also for just the top effects
qtl_output_all_snpsig_openatac_occurences_leads <- data.frame(table(qtl_output_all_snpsig[qtl_output_all_snpsig[['is_top_variant']] == T, c('cell_type', 'in_open_chromatin')]))
qtl_output_all_snpsig_openatac_occurences_leads_only <- data.frame(table(qtl_output_all_snpsig[qtl_output_all_snpsig[['is_top_variant']] == T, c('cell_type', 'in_open_chromatin_only')]))
# order by chromatin state, so that when we plot them later, this is the order they occur in
# qtl_output_all_snpsig_openatac_occurences[['in_open_chromatin']] <- factor(qtl_output_all_snpsig_openatac_occurences[['in_open_chromatin']], levels = c('none', 'unmatched', 'matched'))
# qtl_output_all_snpsig_openatac_occurences_leads[['in_open_chromatin']] <- factor(qtl_output_all_snpsig_openatac_occurences_leads[['in_open_chromatin']], levels = c('none', 'unmatched', 'matched'))
qtl_output_all_snpsig_openatac_occurences <- qtl_output_all_snpsig_openatac_occurences[order(qtl_output_all_snpsig_openatac_occurences[['in_open_chromatin']]), ]
qtl_output_all_snpsig_openatac_occurences_leads <- qtl_output_all_snpsig_openatac_occurences_leads[order(qtl_output_all_snpsig_openatac_occurences_leads[['in_open_chromatin']]), ]
qtl_output_all_snpsig_openatac_occurences_only <- qtl_output_all_snpsig_openatac_occurences_only[order(qtl_output_all_snpsig_openatac_occurences_only[['in_open_chromatin_only']]), ]
qtl_output_all_snpsig_openatac_occurences_leads_only <- qtl_output_all_snpsig_openatac_occurences_leads_only[order(qtl_output_all_snpsig_openatac_occurences_leads_only[['in_open_chromatin_only']]), ]
# add nicer cell type name
qtl_output_all_snpsig_openatac_occurences[['cell_type_nice']] <- remap_with_label_dict(as.character(qtl_output_all_snpsig_openatac_occurences[['cell_type']]))
qtl_output_all_snpsig_openatac_occurences_leads[['cell_type_nice']] <- remap_with_label_dict(as.character(qtl_output_all_snpsig_openatac_occurences_leads[['cell_type']]))
qtl_output_all_snpsig_openatac_occurences_only[['cell_type_nice']] <- remap_with_label_dict(as.character(qtl_output_all_snpsig_openatac_occurences_only[['cell_type']]))
qtl_output_all_snpsig_openatac_occurences_leads_only[['cell_type_nice']] <- remap_with_label_dict(as.character(qtl_output_all_snpsig_openatac_occurences_leads_only[['cell_type']]))
# add annotation for cell type and chromatin status
qtl_output_all_snpsig_openatac_occurences[['celltype_in_open_chromatin']] <- paste(qtl_output_all_snpsig_openatac_occurences[['cell_type_nice']], qtl_output_all_snpsig_openatac_occurences[['in_open_chromatin']])
qtl_output_all_snpsig_openatac_occurences_leads[['celltype_in_open_chromatin']] <- paste(qtl_output_all_snpsig_openatac_occurences_leads[['cell_type_nice']], qtl_output_all_snpsig_openatac_occurences_leads[['in_open_chromatin']])
qtl_output_all_snpsig_openatac_occurences_only[['celltype_in_open_chromatin']] <- paste(qtl_output_all_snpsig_openatac_occurences_only[['cell_type_nice']], qtl_output_all_snpsig_openatac_occurences_only[['in_open_chromatin_only']])
qtl_output_all_snpsig_openatac_occurences_leads_only[['celltype_in_open_chromatin']] <- paste(qtl_output_all_snpsig_openatac_occurences_leads_only[['cell_type_nice']], qtl_output_all_snpsig_openatac_occurences_leads_only[['in_open_chromatin_only']])
# because we ordered by in_open_chromatin before, if we order based on the current order, the chromatin state with the cell type should follow the same order
qtl_output_all_snpsig_openatac_occurences[['celltype_in_open_chromatin']] <- factor(qtl_output_all_snpsig_openatac_occurences[['celltype_in_open_chromatin']], levels = qtl_output_all_snpsig_openatac_occurences[['celltype_in_open_chromatin']])
qtl_output_all_snpsig_openatac_occurences_leads[['celltype_in_open_chromatin']] <- factor(qtl_output_all_snpsig_openatac_occurences_leads[['celltype_in_open_chromatin']], levels = qtl_output_all_snpsig_openatac_occurences_leads[['celltype_in_open_chromatin']])
qtl_output_all_snpsig_openatac_occurences_only[['celltype_in_open_chromatin']] <- factor(qtl_output_all_snpsig_openatac_occurences_only[['celltype_in_open_chromatin']], levels = qtl_output_all_snpsig_openatac_occurences_only[['celltype_in_open_chromatin']])
qtl_output_all_snpsig_openatac_occurences_leads_only[['celltype_in_open_chromatin']] <- factor(qtl_output_all_snpsig_openatac_occurences_leads_only[['celltype_in_open_chromatin']], levels = qtl_output_all_snpsig_openatac_occurences_leads_only[['celltype_in_open_chromatin']])
# make these into plots
p_all_variant_openatac_overlap <- ggplot(data = qtl_output_all_snpsig_openatac_occurences, mapping = aes(x = cell_type_nice, y = Freq, fill = celltype_in_open_chromatin)) +
  # barplots specifically
  geom_bar(stat = 'identity', position = 'stack') +
  # with manual colors
  scale_fill_manual(values = get_color_coding_dict()) +
  # labels
  xlab('Cell type') + 
  ylab('Number of variants') + 
  ggtitle('Open chromatin state of eSNPs') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # and set the label for the SCREEN annotation
  labs(fill = "Open chromatin state")
p_all_variant_openatac_overlap_lead <- ggplot(data = qtl_output_all_snpsig_openatac_occurences_leads, mapping = aes(x = cell_type_nice, y = Freq, fill = celltype_in_open_chromatin)) +
  # barplots specifically
  geom_bar(stat = 'identity', position = 'stack') +
  # with manual colors
  scale_fill_manual(values = get_color_coding_dict()) +
  # labels
  xlab('Cell type') + 
  ylab('Number of variants') + 
  ggtitle('Open chromatin state of lead eSNPs') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # and set the label for the SCREEN annotation
  labs(fill = "Open chromatin state")
# make these into plots
p_all_variant_openatac_overlap_only <- ggplot(data = qtl_output_all_snpsig_openatac_occurences_only, mapping = aes(x = cell_type_nice, y = Freq, fill = celltype_in_open_chromatin)) +
  # barplots specifically
  geom_bar(stat = 'identity', position = 'stack') +
  # with manual colors
  scale_fill_manual(values = get_color_coding_dict()) +
  # labels
  xlab('Cell type') + 
  ylab('Number of variants') + 
  ggtitle('Open chromatin state of eSNPs') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # and set the label for the SCREEN annotation
  labs(fill = "Open chromatin state")
p_all_variant_openatac_overlap_lead_only <- ggplot(data = qtl_output_all_snpsig_openatac_occurences_leads_only, mapping = aes(x = cell_type_nice, y = Freq, fill = celltype_in_open_chromatin)) +
  # barplots specifically
  geom_bar(stat = 'identity', position = 'stack') +
  # with manual colors
  scale_fill_manual(values = get_color_coding_dict()) +
  # labels
  xlab('Cell type') + 
  ylab('Number of variants') + 
  ggtitle('Open chromatin state of lead eSNPs') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # and set the label for the SCREEN annotation
  labs(fill = "Open chromatin state")

# make with a nicer legend
p_all_variant_openatac_overlap_lead_seplegend <- plot_grid(
  # plot without a legend
  p_all_variant_openatac_overlap_lead + theme(legend.position='none') +
    # label sizes
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 20)
    ), 
  # with a dummy legend
  plot_grid(
    cowplot::get_legend(
      ggplot(
        data = data.frame('open_chromatin_state' = factor(c('none', 'unmatched', 'matched'), levels = c('none', 'unmatched', 'matched')), 'y' = c(1,2,3)), mapping = aes(x = open_chromatin_state, fill = open_chromatin_state, y = y)
      ) + 
        geom_bar(stat = 'identity') + 
        scale_fill_manual(values = list('none' = 'lightgray', 'unmatched' = 'black', 'matched' = 'darkgray')) + 
        labs(fill = "Chromatin\nstate") + 
        theme(legend.title = element_text(size = 18), legend.text = element_text(size = 16))
    ), 
    nrow = 2
  ), 
  ncol = 2, 
  rel_widths = c(0.8, 0.2)
  )
# save the plot
ggsave('~/plots/mo_lead_esnp_chromatin_overlap.pdf', width = 8, height = 8, plot = p_all_variant_openatac_overlap_lead_seplegend)

# redo this, but add the total frequences
qtl_output_all_snpsig_openatac_occurences_total <- aggregate(Freq ~ cell_type, qtl_output_all_snpsig_openatac_occurences[, c('cell_type', 'Freq')], sum)
# add that frequency to the table
qtl_output_all_snpsig_openatac_occurences[['total']] <- qtl_output_all_snpsig_openatac_occurences_total[match(qtl_output_all_snpsig_openatac_occurences[['cell_type']], qtl_output_all_snpsig_openatac_occurences_total[['cell_type']]), ][['Freq']]
# the add the fraction
qtl_output_all_snpsig_openatac_occurences[['frac']] <- qtl_output_all_snpsig_openatac_occurences[['Freq']] / qtl_output_all_snpsig_openatac_occurences[['total']]
# same for the leads
qtl_output_all_snpsig_openatac_occurences_leads_total <- aggregate(Freq ~ cell_type, qtl_output_all_snpsig_openatac_occurences_leads[, c('cell_type', 'Freq')], sum)
qtl_output_all_snpsig_openatac_occurences_leads[['total']] <- qtl_output_all_snpsig_openatac_occurences_leads_total[match(qtl_output_all_snpsig_openatac_occurences_leads[['cell_type']], qtl_output_all_snpsig_openatac_occurences_leads_total[['cell_type']]), ][['Freq']]
qtl_output_all_snpsig_openatac_occurences_leads[['frac']] <- qtl_output_all_snpsig_openatac_occurences_leads[['Freq']] / qtl_output_all_snpsig_openatac_occurences_leads[['total']]
# make the figure again, but plot the fractions instead
p_all_variant_openatac_overlap_lead_frac <- ggplot(data = qtl_output_all_snpsig_openatac_occurences_leads[qtl_output_all_snpsig_openatac_occurences_leads[['in_open_chromatin']] != 'none', ], mapping = aes(x = cell_type_nice, y = frac, fill = celltype_in_open_chromatin)) +
  # barplots specifically
  geom_bar(stat = 'identity', position = 'stack') +
  # with manual colors
  scale_fill_manual(values = get_color_coding_dict()) +
  # labels
  xlab('Cell type') + 
  ylab('Fraction of variants') + 
  ggtitle('Open chromatin state of lead eSNPs') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # and set the label for the SCREEN annotation
  labs(fill = "Open chromatin state") +
  # and set the y axis limits to be between 0 and 1
  ylim(0, 1)
# make with a nicer legend
p_all_variant_openatac_overlap_lead_seplegend_frac <- plot_grid(
  # plot without a legend
  p_all_variant_openatac_overlap_lead_frac + theme(legend.position='none') +
    # label sizes
    theme(
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 20)
    ), 
  # with a dummy legend
  plot_grid(
    cowplot::get_legend(
      ggplot(
        data = data.frame('open_chromatin_state' = factor(c('unmatched', 'matched'), levels = c('unmatched', 'matched')), 'y' = c(1,2)), mapping = aes(x = open_chromatin_state, fill = open_chromatin_state, y = y)
      ) + 
        geom_bar(stat = 'identity') + 
        scale_fill_manual(values = list('unmatched' = 'black', 'matched' = 'darkgray')) + 
        labs(fill = "Chromatin\nstate") + 
        theme(legend.title = element_text(size = 18), legend.text = element_text(size = 16))
    ), 
    nrow = 2
  ), 
  ncol = 2, 
  rel_widths = c(0.8, 0.2)
)
# save the plot
ggsave('~/plots/mo_lead_esnp_chromatin_overlap_frac.pdf', width = 8, height = 8, plot = p_all_variant_openatac_overlap_lead_seplegend_frac)



# get the screen matches for each variant
qtl_variants_all_screen_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_screen_overlap.tsv.gz'
qtl_variants_all_screen <- fread(qtl_variants_all_screen_loc, header = T, sep = '\t')
# add the matching screen region to the variants
qtl_output_all_snpsig[['screen_region']] <- qtl_variants_all_screen[match(qtl_output_all_snpsig[['snp_id']], qtl_variants_all_screen[['snp_id']]), ][['overlapping_feature']]
# get the screen annotations as well
screen_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/GRCh38-cCREs.bed.gz'
screen_anno <- fread(screen_anno_loc, header = F, sep = '\t')
colnames(screen_anno) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
# add name based on the location
screen_anno[['signac_hg38']] <- paste(screen_anno[['chromosome']], screen_anno[['start']], screen_anno[['end']], sep = '-')
# and add the type information onto the QTL output
qtl_output_all_snpsig[['screen_annotation']] <- screen_anno[match(qtl_output_all_snpsig[['screen_region']], screen_anno[['signac_hg38']]), ][['type']]
# make the NA for the annotations into 'none'
qtl_output_all_snpsig[is.na(qtl_output_all_snpsig[['screen_annotation']]), 'screen_annotation'] <- 'none'
# set levels
qtl_output_all_snpsig[['screen_annotation']] <- factor(qtl_output_all_snpsig[['screen_annotation']], levels = c('none', setdiff(unique(qtl_output_all_snpsig[['screen_annotation']]), 'none')))
# now make this into a table per cell type
qtl_output_all_snpsig_screen_occurences <- data.frame(table(qtl_output_all_snpsig[, c('cell_type', 'screen_annotation')]))
# also for just the top effects
qtl_output_all_snpsig_screen_occurences_leads <- data.frame(table(qtl_output_all_snpsig[qtl_output_all_snpsig[['is_top_variant']] == T, c('cell_type', 'screen_annotation')]))
# get colours for each category
qtl_output_all_snpsig_screen_occurences_colours <- roycols::get_color_list(unique(qtl_output_all_snpsig_screen_occurences[['screen_annotation']]))
# but make 'none' gray
qtl_output_all_snpsig_screen_occurences_colours[['none']] <- 'gray'

# make these into plots
p_all_variant_screen_overlap <- ggplot(data = qtl_output_all_snpsig_screen_occurences, mapping = aes(x = cell_type, y = Freq, fill = screen_annotation)) +
  # barplots specifically
  geom_bar(stat = 'identity', position = 'stack') +
  # with manual colors
  scale_fill_manual(values = qtl_output_all_snpsig_screen_occurences_colours) +
  # labels
  xlab('Cell type') + 
  ylab('Number of variants') + 
  ggtitle('Categories of elements from ENCODE SCREEN v4\noverlapping with eSNPs') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # and set the label for the SCREEN annotation
  labs(fill = "SCREEN annotation")
p_all_variant_screen_overlap_lead <- ggplot(data = qtl_output_all_snpsig_screen_occurences_leads, mapping = aes(x = cell_type, y = Freq, fill = screen_annotation)) +
  # barplots specifically
  geom_bar(stat = 'identity', position = 'stack') +
  # with manual colors
  scale_fill_manual(values = qtl_output_all_snpsig_screen_occurences_colours) +
  # labels
  xlab('Cell type') + 
  ylab('Number of variants') + 
  ggtitle('Categories of elements from ENCODE SCREEN v4\noverlapping with lead eSNPs') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  # and set the label for the SCREEN annotation
  labs(fill = "SCREEN annotation")

# add the lineage of each cell
qtl_output_all_snpsig[['lineage']] <- 'none'
qtl_output_all_snpsig[qtl_output_all_snpsig[['cell_type']] %in% c('B', 'CD4T', 'CD8T', 'NK'), ][['lineage']] <- 'lymphoid'
qtl_output_all_snpsig[qtl_output_all_snpsig[['cell_type']] %in% c('DC', 'monocyte'), ][['lineage']] <- 'myeloid'
# get number of values
n_qtl <- nrow(qtl_output_all_snpsig)
# store the matching and non-matching openness
lineage_openness <- data.frame(
  matrix(NA, nrow = n_qtl, ncol = 2, dimnames = list(NULL, c('matching_lineage_openness', 'nonmatching_lineage_openness')))
)
# # set up progress bar
# pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = n_qtl)
# # initialize
# pb$tick(0)
# by checking each row
for (row_i in 1 : n_qtl) {
  # update pb
  # pb$tick(row_i)
  # get the cell type
  ct <- qtl_output_all_snpsig[row_i, 'cell_type']
  
  # get the openness of lymphoid
  lymphoid_openness <- as.vector(unlist(qtl_output_all_snpsig[row_i, paste('openness', c('B', 'CD4T', 'CD8T', 'NK'), sep = '_')]))
  # change each NA to 0
  lymphoid_openness[is.na(lymphoid_openness)] <- 0
  # then get the average
  lymphoid_avg_openness <- mean(lymphoid_openness)
  
  # get the openness of myeloid
  myeloid_openness <- as.vector(unlist(qtl_output_all_snpsig[row_i, paste('openness', c('DC', 'monocyte'), sep = '_')]))
  # change each NA to 0
  myeloid_openness[is.na(myeloid_openness)] <- 0
  # then get the average
  myeloid_avg_openness <- mean(myeloid_openness)
  
  # depending on the cell type, we'll set columns
  if (cell_type %in% c('B', 'CD4T', 'CD8T', 'NK')) {
    lineage_openness[row_i, 'matching_lineage_openness'] <- lymphoid_avg_openness
    lineage_openness[row_i, 'nonmatching_lineage_openness'] <- myeloid_avg_openness
  } else if(cell_type %in% c('DC', 'monocyte')) {
    lineage_openness[row_i, 'matching_lineage_openness'] <- myeloid_avg_openness
    lineage_openness[row_i, 'nonmatching_lineage_openness'] <- lymphoid_avg_openness
  }
}
# add this info
qtl_output_all_snpsig <- cbind(qtl_output_all_snpsig, lineage_openness)
# add other cell types
qtl_output_all_snpsig <- add_significant_celltypes_as_strings_vectorised(qtl_output_all_snpsig, trait2_col = 'feature_id', trait1_col = 'feature_id', variant1_col = 'snp_id', variant2_col = 'snp_id')
# check each to see if they are lineage only
qtl_output_all_snpsig[['lineage_specific']] <- pbapply(qtl_output_all_snpsig, 1, function(x) {
  # get the cell type
  ct <- x[['cell_type']]
  # get the vectorised cts
  other_cts_string <- x[['other_ct']]
  # split to cts
  other_cts <- strsplit(other_cts_string, split = ',')[[1]]
  # check these
  if (ct %in% c('B', 'CD4T', 'CD8T', 'NK') & (
    length(
      setdiff(other_cts, c('B', 'CD4T', 'CD8T', 'NK'))
    ) == 0
  )) {
    return(T)
  } else if (ct %in% c('DC', 'monocyte') & (
    length(
      setdiff(other_cts, c('DC', 'monocyte'))
    ) == 0
  )) {
    return(T)
  } else {
    return(F)
  }
})
# get top variants
qtl_output_all_snpsig_top <- qtl_output_all_snpsig[qtl_output_all_snpsig[['is_top_variant']], ]
# add lineage openness
qtl_output_all_snpsig_top[['nonmatching_lineage_closed']] <- qtl_output_all_snpsig_top[['nonmatching_lineage_openness']] < 0.001
qtl_output_all_snpsig_top[['matching_lineage_open']] <- qtl_output_all_snpsig_top[['matching_lineage_openness']] >= 0.001

