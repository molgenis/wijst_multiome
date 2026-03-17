#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cre_sccre_interaction_confinement.R
# Function: create variant-TF-gene confinements for TF interaction-eQTL mapping
# 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)


####################
# Functions        #
####################

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


get_credible_sets_output <- function(output_table_per_celltype, feature_column='feature_id', cs_column='CS', item_column='cell_type') {
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
    # for the only NA ones, do L1
    ct_eqtl_feature_cs[ct_eqtl_feature_cs[[feature_column]] %in% ct_eqtl_feature_cs_onlyna, cs_column] <- 'L1'
    # then remove the NA ones
    ct_eqtl_feature_cs <- ct_eqtl_feature_cs[!is.na(ct_eqtl_feature_cs[[cs_column]]), ]
    # add the cell type
    ct_eqtl_feature_cs[[item_column]] <- ct
    # put it back in the list
    eqtl_outputs_filtered[[ct]] <- ct_eqtl_feature_cs
  }
  return(eqtl_outputs_filtered)
}


credible_sets_to_numbers_primary_and_other <- function(output_table_per_celltype, feature_column='feature_id', cs_column='CS') {
  # create a table to put results in
  nrs_table <- data.frame(
    'celltype' = rep(NA, times = length(names(output_table_per_celltype)) * 2), 
    'effect' = rep(NA, times = length(names(output_table_per_celltype)) * 2), 
    'number' = rep(NA, times = length(names(output_table_per_celltype)) * 2)
  )
  # check each cell type
  for (i in 1 : length(names(output_table_per_celltype))) {
    # extract cell type
    cell_type <- names(output_table_per_celltype)[i]
    # get the table
    ct_table <- output_table_per_celltype[[cell_type]]
    # get the unique features
    unique_features_ct <- unique(ct_table[[feature_column]])
    # get how many effects
    n_effects <- nrow(ct_table)
    # and how many unique effects
    n_unique_effects <- length(unique_features_ct)
    # which means there are non-primary effects
    n_non_primary_effects <- n_effects - n_unique_effects
    # put that into the table
    nrs_table[(i * 2 - 1), 'celltype'] <- cell_type
    nrs_table[(i * 2), 'celltype'] <- cell_type
    nrs_table[(i * 2 - 1), 'effect'] <- 'primary'
    nrs_table[(i * 2), 'effect'] <- 'additional'
    nrs_table[(i * 2 - 1), 'number'] <- n_unique_effects
    nrs_table[(i * 2), 'number'] <- n_non_primary_effects
    print(nrs_table)
  }
  return(nrs_table)
}


cs_number_to_name <- function(cs_number_or_name) {
  # this is our list with the mapping
  name_mapping <- list()
  name_mapping[['1']] <- 'primary'
  name_mapping[['2']] <- 'secondary'
  name_mapping[['3']] <- 'tertiary'
  name_mapping[['4']] <- 'quaternary'
  name_mapping[['5']] <- 'quinary'
  name_mapping[['6']] <- 'senary'
  name_mapping[['7']] <- 'septenary'
  name_mapping[['8']] <- 'octonary'
  name_mapping[['9']] <- 'nonary'
  name_mapping[['10']] <- 'senary'
  # now also do them with the 'L' added
  for (effect in names(name_mapping)) {
    name_mapping[[paste0('l', effect)]] <- name_mapping[[effect]]
    name_mapping[[paste0('L', effect)]] <- name_mapping[[effect]]
  }
  # check if the number or name was in the list
  if (as.character(cs_number_or_name) %in% names(name_mapping)) {
    return(name_mapping[[as.character(cs_number_or_name)]])
  }
  # otherwise just return the original
  else {
    # with a warning
    warning(paste('supplied credible set number or name does not have a string representation, value is returned as unchanged string'))
    return(as.character(cs_number_or_name))
  }
}


credible_sets_to_effect_number <- function(output_table_per_celltype, feature_column='feature_id', cs_column='CS') {
  # we'll store per cell type first
  effect_numbers_per_celltype <- list()
  # we'll check each cell type
  for (cell_type in names(output_table_per_celltype)) {
    # extract for that cell type
    ct_table <- output_table_per_celltype[[cell_type]]
    # now get the number for each credible set
    ct_cs_numbers <- data.frame(table(ct_table[[cs_column]]))
    # order that by the number
    ct_cs_numbers <- ct_cs_numbers[order(ct_cs_numbers[['Freq']], decreasing = T), ]
    # sometimes we'll have L4, but not L3. To account for that, we'll just make a new ordering
    ct_cs_numbers[['effect_nr']] <- 1:nrow(ct_cs_numbers)
    # and then give that a nice name
    ct_cs_numbers[['effect']] <- apply(ct_cs_numbers, 1, function(x){cs_number_to_name(x['effect_nr'])})
    # add the cell type
    ct_cs_numbers[['celltype']] <- cell_type
    # now filter to just what we care about
    ct_cs_numbers <- ct_cs_numbers[, c('celltype', 'effect', 'Freq')]
    # and rename columns to be how we like
    colnames(ct_cs_numbers) <- c('celltype', 'effect', 'number')
    # put in the list
    effect_numbers_per_celltype[[cell_type]] <- ct_cs_numbers
  }
  # merge all tables
  effect_numbers_all <- do.call('rbind', effect_numbers_per_celltype)
  return(effect_numbers_all)
}

credible_sets_to_numbers <- function(output_table_per_celltype, feature_column='feature_id', cs_column='CS', split_effects=F) {
  # initialize variable
  nrs_table <- NULL
  # we might care about primary, secondary, tertiary etc.
  if (split_effects) {
    nrs_table <- credible_sets_to_effect_number(output_table_per_celltype, feature_column = feature_column, cs_column = cs_column)
  }
  # or we just want primary and the rest
  else {
    nrs_table <- credible_sets_to_numbers_primary_and_other(output_table_per_celltype, feature_column = feature_column, cs_column = cs_column)
  }
  return(nrs_table)
}


#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @returns a list with the output tables per cell type
#' 
get_output_per_celltype_coloc <- function(coloc_output_loc, file_prepend='', file_append='_eqtl_caqtl_coloc.tsv.gz', h4_column='PP.H4.abf') {
  # get the folders in the directory, which should be the cell types
  cell_type_files <- list.files(coloc_output_loc, full.names = F, recursive = F)
  # filter these on the pattern
  cell_type_pattern <- paste0(file_prepend, '.*', file_append, '$')
  cell_type_files <- cell_type_files[grep(cell_type_pattern, cell_type_files)]
  # store the result per cell type
  coloc_per_celltype <- list()
  # check each file
  for (cell_type_file in cell_type_files) {
    # paste the full path together
    cell_type_file_full <- paste0(coloc_output_loc, '/', cell_type_file)
    # read the file
    coloc_celltype <- fread(cell_type_file_full, header = T, sep = '\t')
    # filter on existing h4
    coloc_celltype <- coloc_celltype[!is.na(coloc_celltype[[h4_column]]), ]
    # remove the prepend from the file
    cell_type <- gsub(file_prepend, '', cell_type_file)
    # and the append
    cell_type <- gsub(file_append, '', cell_type)
    # add the cell type to the table
    coloc_celltype[['cell_type']] <- cell_type
    # and put in list
    coloc_per_celltype[[cell_type]] <- coloc_celltype
  }
  return(coloc_per_celltype)
}


get_top_effect_per_cs <- function(cs_output_per_ct, feature_column='feature_id', variant_column='snp_id', sort_column='p_value', cs_column='CS', decreasing_sort=F) {
  # store the variant-feature confinement per ct
  var_feature_per_ct <- list()
  # check each cell type
  for (cell_type in names(cs_output_per_ct)) {
    # get output for this cell type
    cs_output_ct <- cs_output_per_ct[[cell_type]]
    # order
    cs_output_ct <- cs_output_ct[order(cs_output_ct[[sort_column]], decreasing = decreasing_sort), ]
    # check if there are entries where the CS is not known
    cs_output_ct_cs_na <- is.na(cs_output_ct[[cs_column]])
    # if there are more than zero
    if (length(cs_output_ct_cs_na) > 0) {
      # set the CS as 'U'
      cs_output_ct[cs_output_ct_cs_na, ][[cs_column]] <- 'U'
      # get the top variant per feature and cs
      cs_output_ct <- cs_output_ct[!duplicated(paste(cs_output_ct[[feature_column]], cs_output_ct[[cs_column]])), ]
      # get the entries where the CS is 'U'
      cs_output_ct_u <- cs_output_ct[cs_output_ct[[cs_column]] == 'U', ]
      # count the credible sets with unknown credible sets
      cs_output_ct_u_ct_counts <- data.frame(table(cs_output_ct_u[[feature_column]]))
      # get the ones that have a credible set in addition to the U one
      cs_output_ct_u_ct_also_not_u <- cs_output_ct_u_ct_counts[cs_output_ct_u_ct_counts[['Freq']] > 1, ]
      # check if there are any
      if (nrow(cs_output_ct_u_ct_also_not_u) > 1) {
        # get the features
        cs_output_ct_u_ct_also_not_u_features <- cs_output_ct_u_ct_also_not_u[[feature_id]]
        # then keep only what is not U and at the same time also not-U
        cs_output_ct <- cs_output_ct[
          !(cs_output_ct[[cs_column]] == 'U' & cs_output_ct[[feature_column]] %in% cs_output_ct_u_ct_also_not_u_features), 
        ]
      }
    }
    # get the variant and feature
    top_var_feature_per_cs <- data.frame('cell_type' = rep(cell_type, times = nrow(cs_output_ct)), 'variant' = cs_output_ct[[variant_column]], 'feature' = cs_output_ct[[feature_column]])
    # put in the list
    var_feature_per_ct[[cell_type]] <- top_var_feature_per_cs
  }
  # merge all
  var_feature_all <- do.call('rbind', var_feature_per_ct)
  return(var_feature_all)
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

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
eqtl_outputs <- get_output_per_celltype_limix(eqtl_output_loc)
# get the top effects per credible set
eqtl_top_var_feature_per_cs <- get_top_effect_per_cs(eqtl_outputs)

# location of the QTL outputs
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/caqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
caqtl_outputs <- get_output_per_celltype_limix(caqtl_output_loc)
# get the top effects per credible set
caqtl_top_var_feature_per_cs <- get_top_effect_per_cs(caqtl_outputs)

# the location of SCENIC output
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# remove what SCENIC thinks is less likely
scenic_output <- scenic_output[scenic_output[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# order by the extended
scenic_output <- scenic_output[order(scenic_output[['is_extended']]), ]
# get which are non-extended if it was both extended and non-extended
scenic_eregs <- unique(scenic_output[, c('TF', 'Gene_signature_direction', 'Gene_signature_name', 'source')])
scenic_eregs <-scenic_eregs[order(scenic_eregs[['source']]), ]
scenic_eregs_to_keep <- scenic_eregs[!duplicated(paste(scenic_eregs[['TF']], scenic_eregs[['Gene_signature_direction']])), ]
# then use that to filer
scenic_output <- scenic_output[scenic_output[['Gene_signature_name']] %in% scenic_eregs_to_keep[['Gene_signature_name']], ]
# rename the regions
scenic_output[['region_cpeaks']] <- gsub(':', '-', scenic_output[['Region']])

# add the TF to the top QTLs
eqtl_top_var_feature_per_cs_tf <- merge(eqtl_top_var_feature_per_cs, scenic_output[, c('Gene', 'Gene_signature_name')], by.x = 'feature', by.y = 'Gene')
# keep across all cell types
eqtl_top_var_feature_per_cs_tf_all <- unique(eqtl_top_var_feature_per_cs_tf[, c('variant', 'Gene_signature_name', 'feature')])
# set column names
colnames(eqtl_top_var_feature_per_cs_tf_all) <- c('variant', 'eregulon', 'feature')

# add the TF to the top caQTLs
caqtl_top_var_feature_per_cs_tf <- merge(caqtl_top_var_feature_per_cs, scenic_output[, c('region_cpeaks', 'Gene', 'Gene_signature_name')], by.x = 'feature', by.y = 'region_cpeaks')
# we will ook at the gene-TF relations, so we need to go from region to gene here
caqtl_top_var_feature_per_cs_tf_all <- unique(caqtl_top_var_feature_per_cs_tf[, c('variant', 'Gene_signature_name', 'Gene')])
# set column names
colnames(caqtl_top_var_feature_per_cs_tf_all) <- c('variant', 'eregulon', 'feature')

# merge the gene-TF results from the eQTLs and caQTLs
qtl_top_var_feature_per_cs_tf_all_both <- unique(rbind(
  eqtl_top_var_feature_per_cs_tf_all, 
  caqtl_top_var_feature_per_cs_tf_all
))
# sort to make comparisons easier
qtl_top_var_feature_per_cs_tf_all_both <- qtl_top_var_feature_per_cs_tf_all_both[order(qtl_top_var_feature_per_cs_tf_all_both[['variant']], qtl_top_var_feature_per_cs_tf_all_both[['eregulon']], qtl_top_var_feature_per_cs_tf_all_both[['feature']]), ]

# write this file
tf_to_gene_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement_inclcaqtls.tsv.gz'
write.table(qtl_top_var_feature_per_cs_tf_all_both, gzfile(tf_to_gene_output_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_sha256_for_file(tf_to_gene_output_loc)

# add scenic region info
eqtl_outputs_scenic_regions <- add_associated_region(eqtl_outputs, scenic_output, scenic_region_column = 'region_cpeaks')
# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])
# set names
colnames(cpeaks_anno) <- paste('atac', colnames(cpeaks_anno), sep = '_')
# store filtered regions
eqtl_outputs_scenic_regions_filtered <- list()
# the region location for each cpeaks region
for (ct in names(eqtl_outputs_scenic_regions)) {
  # get the table
  eqtl_outputs_scenic_regions_ct <- eqtl_outputs_scenic_regions[[ct]]
  # add to the overlap table
  eqtl_outputs_scenic_regions_ct <- cbind(eqtl_outputs_scenic_regions_ct, cpeaks_anno[match(eqtl_outputs_scenic_regions_ct[['region']], cpeaks_anno[['atac_signac_hg38']]), c('atac_chr_hg38', 'atac_start_hg38', 'atac_end_hg38', 'atac_housekeeping', 'atac_screen_all')])
  # get the distances
  eqtl_outputs_scenic_regions_ct_gene_distances <- get_closest_flanks(eqtl_outputs_scenic_regions_ct, 'atac_start_hg38', 'atac_end_hg38', 'feature_start', 'feature_end')
  # add that to the original table
  eqtl_outputs_scenic_regions_ct[['region_gene_distance']] <- eqtl_outputs_scenic_regions_ct_gene_distances[['min_dist']]
  # get the distances
  eqtl_outputs_scenic_regions_ct_esnp_distances <- get_closest_flanks(eqtl_outputs_scenic_regions_ct, 'atac_start_hg38', 'atac_end_hg38', 'snp_position', 'snp_position')
  # add that to the original table
  eqtl_outputs_scenic_regions_ct[['region_esnp_distance']] <- eqtl_outputs_scenic_regions_ct_esnp_distances[['min_dist']]
  # now keep only where the variant is in the chromatin region, and the chromatin region does not overlap the gene
  eqtl_outputs_scenic_regions_ct <- eqtl_outputs_scenic_regions_ct[
    !is.na(eqtl_outputs_scenic_regions_ct[['region_gene_distance']]) &
      !is.na(eqtl_outputs_scenic_regions_ct[['region_esnp_distance']]) &
      eqtl_outputs_scenic_regions_ct[['region_gene_distance']] > 0 &
      eqtl_outputs_scenic_regions_ct[['region_esnp_distance']] == 0,  
  ]
  # put that back into the list
  eqtl_outputs_scenic_regions_filtered[[ct]] <- eqtl_outputs_scenic_regions_ct
}
# get the top effects per credible set again, after filtering
eqtl_top_var_feature_per_cs_filtered <- get_top_effect_per_cs(eqtl_outputs_scenic_regions_filtered)
# merge the output
eqtl_outputs_scenic_regions_filtered_merged <- rbindlist(eqtl_outputs_scenic_regions_filtered)
# and add the region to the top effects
eqtl_top_var_feature_per_cs_filtered[['region']] <- eqtl_outputs_scenic_regions_filtered_merged[
  match(paste(eqtl_top_var_feature_per_cs_filtered[['variant']], eqtl_top_var_feature_per_cs_filtered[['feature']]), 
        paste(eqtl_outputs_scenic_regions_filtered_merged[['snp_id']], eqtl_outputs_scenic_regions_filtered_merged[['feature_id']])
        ), 
][['region']]
# add the TF to the top QTLs
eqtl_top_var_feature_per_cs_tf <- merge(eqtl_top_var_feature_per_cs_filtered, scenic_output[, c('region_cpeaks', 'Gene', 'Gene_signature_name')], by.x = c('region', 'feature'), by.y = c('region_cpeaks', 'Gene'))
# keep across all cell types
eqtl_top_var_feature_per_cs_tf_all <- unique(eqtl_top_var_feature_per_cs_tf[, c('variant', 'Gene_signature_name', 'feature')])
# set column names
colnames(eqtl_top_var_feature_per_cs_tf_all) <- c('variant', 'eregulon', 'feature')
# merge the gene-TF results from the eQTLs and caQTLs
qtl_top_var_feature_per_cs_tf_all_both <- unique(rbind(
  eqtl_top_var_feature_per_cs_tf_all, 
  caqtl_top_var_feature_per_cs_tf_all
))
# sort to make comparisons easier
qtl_top_var_feature_per_cs_tf_all_both <- qtl_top_var_feature_per_cs_tf_all_both[order(qtl_top_var_feature_per_cs_tf_all_both[['variant']], qtl_top_var_feature_per_cs_tf_all_both[['eregulon']], qtl_top_var_feature_per_cs_tf_all_both[['feature']]), ]
# write this file
tf_to_gene_filtered_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement_inclcaqtls_varinregion.tsv.gz'
write.table(qtl_top_var_feature_per_cs_tf_all_both, gzfile(tf_to_gene_filtered_output_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_sha256_for_file(tf_to_gene_filtered_output_loc)
