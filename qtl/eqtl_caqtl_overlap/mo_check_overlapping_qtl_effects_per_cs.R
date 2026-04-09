#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_overlapping_qtl_effects_per_cs.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)
library(ggplot2)


####################
# Functions        #
####################

get_coloc_output <- function(coloc_output_loc, append='_eqtl_caqtl_coloc.tsv.gz') {
  # list all of the files in the directory
  coloc_files <- list.files(coloc_output_loc, recursive = F, full.names = F)
  # filter to files that fit the pattern
  coloc_files <- coloc_files[grepl(paste0('.*', append, '$'), coloc_files)]
  # now store each resulting file into a list
  coloc_per_ct <- list()
  # check each file
  for (coloc_file in coloc_files) {
    # read the file
    coloc_file_loc <- paste(coloc_output_loc, coloc_file, sep = '/')
    coloc_file_output <- fread(coloc_file_loc, header = T, sep = '\t')
    # grab the cell type from the filename
    ct <- gsub(append, '', coloc_file)
    # put that into the table
    coloc_file_output <- cbind(data.table('cell_type' = rep(ct, times = nrow(coloc_file_output))), coloc_file_output)
    # put in the list
    coloc_per_ct[[ct]] <- coloc_file_output
  }
  # merge all
  coloc_all <- do.call('rbind', coloc_per_ct)
  return(coloc_all)
}


get_coloc_variants <- function(finemapped_colocs, cell_type_column='ct', feature_column='Gene', variant_column='hit1', PP_H4_abf_column='PP.H4.abf', PP_H4_abf_column_cutoff=0.8, link_split_char='_') {
  # subset to significant
  finemapped_colocs_colocing <- finemapped_colocs[finemapped_colocs[[PP_H4_abf_column]] >= PP_H4_abf_column_cutoff, ]
  # extract the variant
  finemapped_colocs_colocing_short <- data.frame('variant' = finemapped_colocs_colocing[[variant_column]])
  # extract the feature link
  finemapped_colocs_colocing_short[, c('gene', 'region')] <- str_split_fixed(finemapped_colocs_colocing[[feature_column]], link_split_char, 2)
  # and the cell type
  finemapped_colocs_colocing_short[['cell_type']] <- finemapped_colocs_colocing[[cell_type_column]]
  return(finemapped_colocs_colocing_short)
}


#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @param gene_column which column to use as the gene identifier
#' @returns a list with the output tables per cell type
#' 
get_output_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all_qval_allchroms_fdr005_significant_cs.tsv.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, verbose=T, add_nominal_cutoff=T, nominal_p_value_column='p_value') {
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
      # add nominal cutoff by first getting the top effects
      cell_type_output_top <- cell_type_output[order(cell_type_output[[nominal_p_value_column]]), ]
      cell_type_output_top <- cell_type_output_top[!duplicated(cell_type_output_top[[gene_column]]), ]
      # get the highest still significant p value
      max_sig_p <- max(cell_type_output_top[[nominal_p_value_column]])
      # and add that to the output
      cell_type_output[['nominal_p_value_cutoff']] <- max_sig_p
    }
    # add to the list
    egenes_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(egenes_per_celltype)
}


add_credible_sets_output <- function(output_table_per_celltype, feature_column='feature_id', cs_column='CS', item_column='cell_type') {
  # filter on pval_nominal_threshold_global
  eqtl_outputs_filtered <- list()
  for (ct in names(output_table_per_celltype)) {
    # extract this eQTL output
    ct_eqtl_output <- output_table_per_celltype[[ct]]
    # # filter on the pval threshold
    # ct_eqtl_output <- ct_eqtl_output[ct_eqtl_output[['p_value']] < ct_eqtl_output[['pval_nominal_threshold_global']], ]
    # get the unique combinations of features and credible sets
    ct_eqtl_feature_cs <- unique(ct_eqtl_output[, c(feature_column, cs_column)])
    # get features with NA values
    ct_eqtl_feature_cs_na <- ct_eqtl_feature_cs[is.na(ct_eqtl_feature_cs[[cs_column]]), feature_column]
    # and the ones without
    ct_eqtl_feature_cs_nona <- ct_eqtl_feature_cs[!is.na(ct_eqtl_feature_cs[[cs_column]]), feature_column]
    # get the ones which are only with na
    ct_eqtl_feature_cs_onlyna <- setdiff(ct_eqtl_feature_cs_na, ct_eqtl_feature_cs_nona)
    if (length(ct_eqtl_feature_cs_onlyna)) {
      # for those, set the CS to simply be L1
      ct_eqtl_output[ct_eqtl_output[[feature_column]] %in% ct_eqtl_feature_cs_onlyna, ][[cs_column]] <- 'L1'
    }
    # get the ones which are partly na
    ct_eqtl_feature_cs_somena <- intersect(ct_eqtl_feature_cs_nona, ct_eqtl_feature_cs_na)
    if (length(ct_eqtl_feature_cs_somena) > 0) {
      ct_eqtl_output[(ct_eqtl_output[[feature_column]] %in% ct_eqtl_feature_cs_somena) & is.na(ct_eqtl_output[[cs_column]]), ][[cs_column]] <- 'L0'
    }
    # add ct
    ct_eqtl_output[[item_column]] <- ct
    # put back in list
    eqtl_outputs_filtered[[ct]] <- ct_eqtl_output
  }
  return(eqtl_outputs_filtered)
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
    egene_numbers[[celltype_column]] <- remap_with_label_dict(as.character(egene_numbers[[celltype_column]]))
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

# location of the finemapped results
finemapped_colocs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/colocalization/eqtl_caqtl/ut_and_24hca_significant/'
# read the tabel
finemapped_colocs <- get_coloc_output(finemapped_colocs_loc)
# add the way that we merged these
finemapped_colocs[['method']] <- 'coloc'
# filter the colocs
finemapped_colocs <- finemapped_colocs[!is.na(finemapped_colocs[['PP.H4.abf']]) & finemapped_colocs[['PP.H4.abf']] >= 0.75, ]

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
eqtl_outputs <- get_output_per_celltype_limix(eqtl_output_loc)
# location of the QTL outputs
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/caqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
caqtl_outputs <- get_output_per_celltype_limix(caqtl_output_loc)

# add cs info
eqtl_outputs_cs <- add_credible_sets_output(eqtl_outputs)
caqtl_outputs_cs <- add_credible_sets_output(caqtl_outputs)

# merge all of them
eqtl_outputs_cs_all <- rbindlist(eqtl_outputs_cs, fill = T)
caqtl_outputs_cs_all <- rbindlist(caqtl_outputs_cs, fill = T)
# remove variants not at threshold
eqtl_outputs_cs_all <- eqtl_outputs_cs_all[eqtl_outputs_cs_all[['p_value']] <= eqtl_outputs_cs_all[['pval_nominal_threshold_global']], ]
# fix significance column
caqtl_outputs_cs_all[is.na(caqtl_outputs_cs_all[['pval_nominal_threshold_global']]), ][['pval_nominal_threshold_global']] <- caqtl_outputs_cs_all[is.na(caqtl_outputs_cs_all[['pval_nominal_threshold_global']]), ][['nominal_p_value_cutoff']]
caqtl_outputs_cs_all <- caqtl_outputs_cs_all[caqtl_outputs_cs_all[['p_value']] <= caqtl_outputs_cs_all[['pval_nominal_threshold_global']], ]
# get variants for each
eqtl_outputs_cs_all_variants <- eqtl_outputs_cs_all[['snp_id']]
caqtl_outputs_cs_all_variants <- caqtl_outputs_cs_all[['snp_id']]
# and cell types
eqtl_outputs_cs_all_cts <- eqtl_outputs_cs_all[['cell_type']]
caqtl_outputs_cs_all_cts <- caqtl_outputs_cs_all[['cell_type']]
# remove the columns
eqtl_outputs_cs_all[['snp_id']] <- NULL
caqtl_outputs_cs_all[['snp_id']] <- NULL
eqtl_outputs_cs_all[['cell_type']] <- NULL
caqtl_outputs_cs_all[['cell_type']] <- NULL
# rename so that the e and ca are added to the names
colnames(eqtl_outputs_cs_all) <- paste('e', colnames(eqtl_outputs_cs_all), sep = '_')
colnames(caqtl_outputs_cs_all) <- paste('ca', colnames(caqtl_outputs_cs_all), sep = '_')
# add variants back
eqtl_outputs_cs_all <- cbind(data.table('cell_type' = eqtl_outputs_cs_all_cts, 'variant_id' = eqtl_outputs_cs_all_variants), eqtl_outputs_cs_all)
caqtl_outputs_cs_all <- cbind(data.table('cell_type' = caqtl_outputs_cs_all_cts, 'variant_id' = caqtl_outputs_cs_all_variants), caqtl_outputs_cs_all)
# merge the outputs based on overlapping variant (for that cell type)
qtl_outputs_cs_overlapping <- merge(x = caqtl_outputs_cs_all, y = eqtl_outputs_cs_all, by = c('cell_type', 'variant_id'))

# add the CS info
finemapped_colocs[['e_CS']] <- eqtl_outputs_cs_all[
  match(
    paste(finemapped_colocs[['hit2']], finemapped_colocs[['trait2']]), 
    paste(eqtl_outputs_cs_all[['variant_id']], eqtl_outputs_cs_all[['e_feature_id']])
  ), 
][['e_CS']]

# get the unique CS-region-gene
finemapped_colocs_unique_overlap <- unique(finemapped_colocs[, c('trait1', 'trait2', 'e_CS')])
qtl_outputs_cs_unique_overlapping <- unique(qtl_outputs_cs_overlapping[, c('ca_feature_id', 'e_feature_id', 'e_CS')])
# set the same column names
colnames(finemapped_colocs_unique_overlap) <- c('region', 'gene', 'cs')
colnames(qtl_outputs_cs_unique_overlapping) <- c('region', 'gene', 'cs')
# rbind
colocs_unique_overlap_both <- unique(rbind(
  finemapped_colocs_unique_overlap, 
  qtl_outputs_cs_unique_overlapping
))
# for now, just keep the L0 as L1
colocs_unique_overlap_both[['cs']] <- gsub('L0', 'L1', colocs_unique_overlap_both[['cs']])
# do unique again
colocs_unique_overlap_both <- unique(colocs_unique_overlap_both)

# get unique overlaps
n_overlapping <- nrow(colocs_unique_overlap_both)

# get unique gene and CS combinations that are in the overlap
length(unique(paste(colocs_unique_overlap_both$gene, colocs_unique_overlap_both$cs)))
# 1753
# get the number of unique CS-gene pairs
eqtl_outputs_cs_all_nol0 <- eqtl_outputs_cs_all
eqtl_outputs_cs_all_nol0[['e_CS']] <- gsub('L0', 'L1', eqtl_outputs_cs_all_nol0[['e_CS']])
length(unique(paste(eqtl_outputs_cs_all_nol0$e_feature_id, eqtl_outputs_cs_all_nol0$e_CS)))
# 560

# add new credible set column
qtl_outputs_cs_overlapping[['e_CS_noL0']] <- gsub('L0', 'L1', qtl_outputs_cs_overlapping[['e_CS']])
# get the number of unique region-gene-cs combinations in the overlap
n_gene_cs_region_tbl <- data.frame(table(unique(qtl_outputs_cs_overlapping[, c('cell_type', 'ca_feature_id', 'e_feature_id', 'e_CS_noL0')])[['cell_type']]))
# set columns
colnames(n_gene_cs_region_tbl) <- c('cell_type', 'n_overlapping')
# plot these numbers
p_ngene_and_cs_to_peaks <- egene_numbers_to_plot(n_gene_cs_region_tbl, celltype_column='cell_type', number_column='n_overlapping', use_label_dict = T, use_color_dict = T, split_long_labels_to_lines = F, use_distinct_colours = F, legendless = T) +
  xlab('Cell type') +
  ylab('Number of independent eQTL effects and caQTL overlaps')
# show plot
p_ngene_and_cs_to_peaks
# and save
ggsave(filename = '~/multiome/plots/mo_egene_cs_caqtl_overlap_numbers.pdf', plot = p_ngene_and_cs_to_peaks, width = 5, height = 5)

# add credible set column also to original overlap table
eqtl_outputs_cs_all[['e_CS_noL0']] <- gsub('L0', 'L1', eqtl_outputs_cs_all[['e_CS']])
# get number of unique gene-CS combinations in the original eQTL output
n_gene_cs_tbl <- data.frame(table(unique(eqtl_outputs_cs_all[, c('cell_type', 'e_feature_id', 'e_CS_noL0')])[['cell_type']]))
# set columns
colnames(n_gene_cs_tbl) <- c('cell_type', 'n_total')
# get unique gene-CS combinations in the overlap
n_gene_cs_overlap_tbl <- data.frame(table(unique(qtl_outputs_cs_overlapping[, c('cell_type', 'e_feature_id', 'e_CS_noL0')])[['cell_type']]))
# set columns
colnames(n_gene_cs_overlap_tbl) <- c('cell_type', 'n_overlapping')
# now merge both
n_gene_cs_merged <- merge(x = n_gene_cs_tbl, y = n_gene_cs_overlap_tbl, by = 'cell_type')
# calculate percentage of gene-CS combinations that are in the overlap
n_gene_cs_merged[['frac_with_caqtl']] <- n_gene_cs_merged[['n_overlapping']] / n_gene_cs_merged[['n_total']]
# plot these numbers
p_frac_gene_and_cs_wcaqtl <- egene_numbers_to_plot(n_gene_cs_merged, celltype_column='cell_type', number_column='frac_with_caqtl', use_label_dict = T, use_color_dict = T, split_long_labels_to_lines = F, use_distinct_colours = F, legendless = T) + 
  xlab('Cell type') + 
  ylab('Fraction of independent eQTL effects with caQTL overlap')
# show the plot
p_frac_gene_and_cs_wcaqtl
# and save
ggsave(filename = '~/multiome/plots/mo_egene_cs_with_caqtl_overlap_fraction.pdf', plot = p_frac_gene_and_cs_wcaqtl, width = 5, height = 5)

# get number of unique gene-CS combinations in the original eQTL output
n_gene_tbl <- data.frame(table(unique(eqtl_outputs_cs_all[, c('cell_type', 'e_feature_id')])[['cell_type']]))
# set columns
colnames(n_gene_tbl) <- c('cell_type', 'n_total')
# get unique gene-CS combinations in the overlap
n_gene_overlap_tbl <- data.frame(table(unique(qtl_outputs_cs_overlapping[, c('cell_type', 'e_feature_id')])[['cell_type']]))
# set columns
colnames(n_gene_overlap_tbl) <- c('cell_type', 'n_overlapping')
# now merge both
n_gene_merged <- merge(x = n_gene_tbl, y = n_gene_overlap_tbl, by = 'cell_type')
# calculate percentage of gene-CS combinations that are in the overlap
n_gene_merged[['frac_with_caqtl']] <- n_gene_merged[['n_overlapping']] / n_gene_merged[['n_total']]
# plot these numbers
p_frac_gene_wcaqtl <- egene_numbers_to_plot(n_gene_merged, celltype_column='cell_type', number_column='frac_with_caqtl', use_label_dict = T, use_color_dict = T, split_long_labels_to_lines = F, use_distinct_colours = F, legendless = T) + 
  xlab('Cell type') + 
  ylab('Fraction of eGenes with caQTL overlap')
# show the plot
p_frac_gene_wcaqtl
# and save
ggsave(filename = '~/multiome/plots/mo_egene_with_caqtl_overlap_fraction.pdf', plot = p_frac_gene_wcaqtl, width = 5, height = 5)

# get the number of unique eGenes
egenes_all <- unique(eqtl_outputs_cs_all[['e_feature_id']])
# get which are are at any point overlapping with a caQTL
egenes_all_caqtl_overlap <- egenes_all[egenes_all %in% qtl_outputs_cs_overlapping[['e_feature_id']]]
# and what the fraction is
length(egenes_all_caqtl_overlap) / length(egenes_all)
# 0.8523234


# add zscores
qtl_outputs_cs_overlapping[['e_z']] <- qtl_outputs_cs_overlapping[['e_beta']] / qtl_outputs_cs_overlapping[['e_beta_se']]
qtl_outputs_cs_overlapping[['ca_z']] <- qtl_outputs_cs_overlapping[['ca_beta']] / qtl_outputs_cs_overlapping[['ca_beta_se']]
# and absolute z score to make things easier for Monique
qtl_outputs_cs_overlapping[['e_z_abs']] <- abs(qtl_outputs_cs_overlapping[['e_z']])
qtl_outputs_cs_overlapping[['ca_z_abs']] <- abs(qtl_outputs_cs_overlapping[['ca_z']])

# export specifically the LY86 examples Monique is interested in
ly86_overlaps <- qtl_outputs_cs_overlapping[qtl_outputs_cs_overlapping[['e_feature_id']] == 'LY86' & 
                                              qtl_outputs_cs_overlapping[['cell_type']] %in% c('B', 'monocyte'), 
                                            c('cell_type', 'variant_id', 'ca_assessed_allele', 
                                              'ca_feature_id', 'ca_CS', 'ca_pip', 'ca_p_value', 'ca_beta', 'ca_beta_se', 'ca_z', 'ca_z_abs',
                                              'e_feature_id', 'e_CS', 'e_pip', 'e_p_value', 'e_beta', 'e_beta_se', 'e_z', 'e_z_abs')]
# save that specifically
ly86_overlaps_loc <- '~/multiome/tables/mo_ly86_caqtl_eqtl_pairs_b_monocyte.tsv.gz'
write.table(ly86_overlaps, gzfile(ly86_overlaps_loc), sep = '\t', quote = F, row.names = F, col.names = T)
mdfiver::create_sha256_for_file(ly86_overlaps_loc)

# location of the i-eqtl output
ieqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# get all the QTL output
ieqtl_output <- get_qtls_per_celltype_limix(ieqtl_output_loc, output_file = 'inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', gene_column='feature', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column=NULL, nominal_significance_column=NULL)
# merge all the results
ieqtl_output_all <- do.call('rbind', ieqtl_output)
# filter on signifiacnce
ieqtl_output_all_sig <- ieqtl_output_all[ieqtl_output_all[['feature_q_value']] < 0.05 &
                                         ieqtl_output_all[['feature_bf_eigen']] < 0.05, ]
# rename the columns
colnames(ieqtl_output_all_sig) <- c('i_beta', 'i_beta_se', 'i_empirical_feature_p_value', 'i_p_value', 'snp_id', 'feature_id', 'i_n_tests_feature', 'i_feature_bf_eigen', 'i_total_bf_eigen', 'i_feature_q_value', 'i_cell_type')
# merge the iqtl to the qtl output
eqtl_outputs_cs_all <- merge(eqtl_outputs_cs_all, ieqtl_output_all_sig, all.x = T, by.x = c('variant_id', 'e_feature_id', 'cell_type'), by.y = c('snp_id', 'feature_id', 'i_cell_type'))

# location of the i-eqtl output
icaqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# get all the QTL output
icaqtl_output <- get_qtls_per_celltype_limix(icaqtl_output_loc, output_file = '/inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', gene_column='feature', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column=NULL, nominal_significance_column=NULL)
# merge all the results
icaqtl_output_all <- do.call('rbind', icaqtl_output)
# filter on signifiacnce
icaqtl_output_all_sig <- icaqtl_output_all[icaqtl_output_all[['feature_q_value']] < 0.05 &
                                         icaqtl_output_all[['feature_bf_eigen']] < 0.05, ]
# rename the columns
colnames(icaqtl_output_all_sig) <- c('i_beta', 'i_beta_se', 'i_empirical_feature_p_value', 'i_p_value', 'snp_id', 'feature_id', 'i_n_tests_feature', 'i_feature_bf_eigen', 'i_total_bf_eigen', 'i_feature_q_value', 'i_cell_type')
# merge the iqtl to the qtl output
caqtl_outputs_cs_all <- merge(caqtl_outputs_cs_all, icaqtl_output_all_sig, all.x = T, by.x = c('variant_id', 'ca_feature_id', 'cell_type'), by.y = c('snp_id', 'feature_id', 'i_cell_type'))

# add Zs here
eqtl_outputs_cs_all[['e_z']] <- eqtl_outputs_cs_all[['e_beta']] / eqtl_outputs_cs_all[['e_beta_se']]
caqtl_outputs_cs_all[['ca_z']] <- caqtl_outputs_cs_all[['ca_beta']] / caqtl_outputs_cs_all[['ca_beta_se']]
# and absolute z score to make things easier for Monique
eqtl_outputs_cs_all[['e_z_abs']] <- abs(eqtl_outputs_cs_all[['e_z']])
caqtl_outputs_cs_all[['ca_z_abs']] <- abs(caqtl_outputs_cs_all[['ca_z']])

# export specifically the LY86 examples Monique is interested in
eqtl_outputs_cs_ly86 <- eqtl_outputs_cs_all[eqtl_outputs_cs_all[['e_feature_id']] == 'LY86' & 
                                              eqtl_outputs_cs_all[['cell_type']] %in% c('B', 'monocyte'), 
                                            c('cell_type', 'variant_id', 'e_assessed_allele', 
                                              'e_feature_id', 'e_CS', 'e_pip', 'e_p_value', 'e_beta', 'e_beta_se', 'e_z', 'e_z_abs', 
                                              'i_p_value', 'i_beta', 'i_beta_se', 'i_feature_bf_eigen', 'i_feature_q_value')]
# save that specifically
eqtl_outputs_cs_ly86_loc <- '~/multiome/tables/mo_ly86_eqtl_b_monocyte.tsv.gz'
write.table(eqtl_outputs_cs_ly86, gzfile(eqtl_outputs_cs_ly86_loc), sep = '\t', quote = F, row.names = F, col.names = T)
mdfiver::create_sha256_for_file(eqtl_outputs_cs_ly86_loc)
