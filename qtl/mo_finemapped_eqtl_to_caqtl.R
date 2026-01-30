#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_finemapped_eqtl_to_caqtl.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(stringr)
library(UpSetR)


####################
# Settings        #
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

get_betas_file <- function(cell_type_file_loc, confinement_table, cell_type_name, confinement_variant_column='variant', confinement_feature_column='gene', confinement_celltype_column='cell_type', qtl_variant_column='snp_id', qtl_feature_column='feature_id', qtl_beta_column='beta', qtl_effect_column='assessed_allele', qtl_se_column='beta_se', qtl_pvalue_column='p_value', qtl_pvalue_cutoff_column='pval_nominal_threshold_global', qtl_fdr_column='feature_q_value') {
  # read the file
  cell_type_table <- fread(cell_type_file_loc, header = T, sep = '\t')
  # get this for the qtl file as well
  var_features_qtls <- paste(cell_type_table[[qtl_variant_column]], cell_type_table[[qtl_feature_column]], sep = '_')
  # filter the cell type table for these
  if (!is.null(confinement_table)) {
    # check which var-features we want to get
    var_features_celltype <- confinement_table[confinement_table[[confinement_celltype_column]] == cell_type_name, ][['var_feature']]
    # and filter on what we have
    cell_type_table <- cell_type_table[var_features_qtls %in% var_features_celltype, ]
  }
  # extract the data we care about
  cell_type_table_slim <- data.frame('cell_type' = rep(cell_type_name, times = nrow(cell_type_table)), 
                                     'variant' = cell_type_table[[qtl_variant_column]], 
                                     'feature' = cell_type_table[[qtl_feature_column]], 
                                     'effect' = cell_type_table[[qtl_beta_column]], 
                                     'allele' = cell_type_table[[qtl_effect_column]], 
                                     'effect_se' = cell_type_table[[qtl_se_column]], 
                                     'pvalue' = cell_type_table[[qtl_pvalue_column]], 
                                     'pvalue_cutoff' = cell_type_table[[qtl_pvalue_cutoff_column]], 
                                     'qvalue' = cell_type_table[[qtl_fdr_column]])
  return(cell_type_table_slim)
}


get_betas_per_celltype <- function(qtl_output_loc, confinement_table, cell_types=NULL, confinement_variant_column='variant', confinement_feature_column='gene', confinement_celltype_column='cell_type', qtl_file_name='qtl_results_all.txt.gz', qtl_variant_column='snp_id', qtl_feature_column='feature_id', qtl_beta_column='beta', qtl_effect_column='assessed_allele', qtl_se_column='beta_se', qtl_pvalue_column='p_value', qtl_pvalue_cutoff_column='pval_nominal_threshold_global', qtl_fdr_column='feature_q_value') {
  # get the var-gene combinations in the confinement
  if (!is.null(confinement_table)) {
    confinement_table[['var_feature']] <- paste(confinement_table[[confinement_variant_column]], confinement_table[[confinement_feature_column]], sep = '_')
  }
  # list all cell types
  cell_type_folders <- list.dirs(qtl_output_loc, recursive = F, full.names = F)
  # overlap with cell types supplied if done so
  if (!is.null(cell_types)) {
    cell_type_folders <- intersect(cell_type_folders, cell_types)
  }
  # store the cell type result in a list
  cell_types_list <- list()
  # check cell type
  for (cell_type in cell_type_folders) {
    # paste together the full path
    cell_type_file_loc <- paste0(qtl_output_loc, '/', cell_type, '/', qtl_file_name)
    message(paste('reading', cell_type_file_loc, '\n'))
    # get the betas
    cell_type_table_slim <- get_betas_file(cell_type_file_loc, 
                                           confinement_table, 
                                           cell_type_name=cell_type, 
                                           confinement_variant_column=confinement_variant_column, 
                                           confinement_feature_column=confinement_feature_column, 
                                           confinement_celltype_column=confinement_celltype_column, 
                                           qtl_variant_column=qtl_variant_column, 
                                           qtl_feature_column=qtl_feature_column, 
                                           qtl_beta_column=qtl_beta_column, 
                                           qtl_effect_column=qtl_effect_column, 
                                           qtl_se_column=qtl_se_column, 
                                           qtl_pvalue_column=qtl_pvalue_column, 
                                           qtl_pvalue_cutoff_column=qtl_pvalue_cutoff_column, 
                                           qtl_fdr_column=qtl_fdr_column)
    # put in the list
    cell_types_list[[cell_type]] <- cell_type_table_slim
  }
  cell_types_all <- do.call('rbind', cell_types_list)
  return(cell_types_all)
}


merge_beta_files <- function(
    matching_table,
    beta_table1, 
    beta_table2, 
    matching_variant_column='variant', 
    matching_feature1_column='region', 
    matching_feature2_column='gene', 
    matching_celltype_column='cell_type', 
    variant_column1='variant', 
    variant_column2='variant', 
    feature_column1='feature', 
    feature_column2='feature', 
    cell_type_column1='cell_type', 
    cell_type_column2='cell_type', 
    beta_column1='effect', 
    beta_column2='effect', 
    allele_column1='allele', 
    allele_column2='allele') {
  # get the variant-gene-celltype
  matching_set1 <- paste(matching_table[[matching_celltype_column]], matching_table[[matching_variant_column]], matching_table[[matching_feature1_column]])
  # the same for the first beta table
  beta_table1_set1 <- paste(beta_table1[[cell_type_column1]], beta_table1[[variant_column1]], beta_table1[[feature_column1]])
  # get in the right order from the set
  beta_table1_ordered <- beta_table1[match(matching_set1, beta_table1_set1), ]
  
  # get the variant-gene-celltype
  matching_set2 <- paste(matching_table[[matching_celltype_column]], matching_table[[matching_variant_column]], matching_table[[matching_feature2_column]])
  # the same for the decond beta table
  beta_table2_set2 <- paste(beta_table2[[cell_type_column2]], beta_table2[[variant_column2]], beta_table2[[feature_column2]])
  # get in the right order from the set
  beta_table2_ordered <- beta_table2[match(matching_set2, beta_table2_set2), ]
  
  # get the cell type
  cell_types <- matching_table[[matching_celltype_column]]
  # get the betas for the first
  betas1 <- beta_table1_ordered[[beta_column1]]
  # and the second one
  betas2 <- beta_table2_ordered[[beta_column2]]
  
  # drop empty entries
  indices_empty_entries <- is.na(beta_table1_ordered[[beta_column1]]) | is.na(beta_table2_ordered[[beta_column2]])
  if (length(indices_empty_entries) > 0) {
    warning(paste('removing', as.character(length(indices_empty_entries)), 'entries due to the beta being absent\n'))
    matching_table <- matching_table[!indices_empty_entries, ]
    beta_table1_ordered <- beta_table1_ordered[!indices_empty_entries, ]
    beta_table2_ordered <- beta_table2_ordered[!indices_empty_entries, ]
    cell_types <- cell_types[!indices_empty_entries]
    betas1 <- betas1[!indices_empty_entries]
    betas2 <- betas2[!indices_empty_entries]
  }
  
  # and convert beta2 if the alleles don't match
  betas2[beta_table1_ordered[[allele_column1]] != beta_table2_ordered[[allele_column2]]] <- betas2[beta_table1_ordered[[allele_column1]] != beta_table2_ordered[[allele_column2]]] * -1
  
  # make one table
  beta_comparison_table <- data.frame('effect1' = betas1, 'effect2' = betas2)
  beta_comparison_table <- cbind(matching_table, beta_comparison_table)
  return(beta_comparison_table)
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
  # rename if possible
  if (use_label_dict) {
    names(genes_per_ct) <- rename_labels(names(genes_per_ct))
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


# plot the concordance
plot_concordanace <- function(qtl_effect_table, ca_effect_column='ca_effect', e_effect_column='e_effect', main='Effect sizes of caQTLs versus eQTLs') {
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
       xlab = 'caQTL effect size',
       ylab = 'eQTL effect size',
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
# Main Code        #
####################

# location of the finemapped results
finemapped_colocs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/colocalization/eqtl_caqtl/ut_and_24hca_significant/'
# read the tabel
finemapped_colocs <- get_coloc_output(finemapped_colocs_loc)
# location of the eQTL output
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
eqtl_file_name <- 'qtl_results_all_qval_allchroms_fdr005_significant.txt.gz'
# location of the caQTL output
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'
# get the colocs
# significant_finemapped_colocs <- get_coloc_variants(finemapped_colocs)
finemapped_colocs <- finemapped_colocs[finemapped_colocs[['PP.H4.abf']] >= .8, ]
# get the eqtl betas
eqtl_betas <- get_betas_per_celltype(eqtl_output_loc, finemapped_colocs, confinement_variant_column = 'hit2', confinement_feature_column = 'trait2', confinement_celltype_column = 'cell_type', qtl_file_name = eqtl_file_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# and caqtl betas
caqtl_betas <- get_betas_per_celltype(caqtl_output_loc, finemapped_colocs, confinement_variant_column = 'hit1', confinement_feature_column = 'trait1', confinement_celltype_column = 'cell_type', qtl_file_name = eqtl_file_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# filter 
eqtl_betas_unique <- eqtl_betas[
  eqtl_betas[['pvalue']] <= eqtl_betas[['pvalue_cutoff']] &
    eqtl_betas[['qvalue']] < .05, 
]
caqtl_betas_unique <- caqtl_betas[
  caqtl_betas[['pvalue']] <= caqtl_betas[['pvalue_cutoff']] &
    caqtl_betas[['qvalue']] < .05, 
]

# remove empty entries
# significant_finemapped_colocs <- significant_finemapped_colocs[significant_finemapped_colocs[['variant']] != '-', ]
# get eqtl and caqtl betas
# betas_both <- merge_beta_files(significant_finemapped_colocs, caqtl_betas, eqtl_betas)
# set column names for betas
colnames(eqtl_betas_unique) <- c('cell_type', 'e_variant', 'e_feature', 'e_effect', 'e_allele', 'e_effect_se', 'e_pvalue', 'e_pvalue_cutoff', 'e_qvalue')
colnames(caqtl_betas_unique) <- c('cell_type', 'ca_variant', 'ca_feature', 'ca_effect', 'ca_allele', 'ca_effect_se', 'ca_pvalue', 'ca_pvalue_cutoff', 'ca_qvalue')
# merge the eQTLs onto the colocing signals
finemapped_colocs <- merge(finemapped_colocs, eqtl_betas_unique, by.x = c('cell_type', 'hit2', 'trait2'), by.y = c('cell_type', 'e_variant', 'e_feature'))
finemapped_colocs <- merge(finemapped_colocs, caqtl_betas_unique, by.x = c('cell_type', 'hit1', 'trait1'), by.y = c('cell_type', 'ca_variant', 'ca_feature'))
# add the way that we merged these
finemapped_colocs[['method']] <- 'coloc'
#betas_both <- merge(caqtl_betas_unique, eqtl_betas_unique, by = c('cell_type', 'variant'))

# get all betas, regardless of coloc
eqtl_betas_nc <- get_betas_per_celltype(eqtl_output_loc, confinement_table = NULL, confinement_variant_column = 'hit2', confinement_feature_column = 'trait2', confinement_celltype_column = 'cell_type', qtl_file_name = eqtl_file_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
caqtl_betas_nc <- get_betas_per_celltype(caqtl_output_loc, confinement_table = NULL, confinement_variant_column = 'hit1', confinement_feature_column = 'trait1', confinement_celltype_column = 'cell_type', qtl_file_name = eqtl_file_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# filter 
eqtl_betas_nc <- eqtl_betas_nc[
  eqtl_betas_nc[['pvalue']] <= eqtl_betas_nc[['pvalue_cutoff']] &
    eqtl_betas_nc[['qvalue']] < .05, 
]
caqtl_betas_nc <- caqtl_betas_nc[
  caqtl_betas_nc[['pvalue']] <= caqtl_betas_nc[['pvalue_cutoff']] &
    caqtl_betas_nc[['qvalue']] < .05, 
]
# set column names for betas
colnames(eqtl_betas_nc) <- c('cell_type', 'variant2', 'trait2', 'e_effect', 'e_allele', 'e_effect_se', 'e_pvalue', 'e_pvalue_cutoff', 'e_qvalue')
colnames(caqtl_betas_nc) <- c('cell_type', 'variant1', 'trait1', 'ca_effect', 'ca_allele', 'ca_effect_se', 'ca_pvalue', 'ca_pvalue_cutoff', 'ca_qvalue')
# merge these
betas_both_nc <- merge(eqtl_betas_nc, caqtl_betas_nc, by.x = c('cell_type', 'variant2'), by.y = c('cell_type', 'variant1'))
# add the method
betas_both_nc[['method']] <- 'overlap'
# add some missing columns
betas_both_nc[['variant1']] <- betas_both_nc[['variant2']]
betas_both_nc[['hit1']] <- betas_both_nc[['variant1']]
betas_both_nc[['hit2']] <- betas_both_nc[['variant1']]
betas_both_nc[['dataset1']] <- 'caQTLs'
betas_both_nc[['dataset2']] <- 'eQTLs'
# add other columns that don't have a value
for (missing_column in setdiff(colnames(finemapped_colocs), colnames(betas_both_nc))) {
  betas_both_nc[[missing_column]] <- NA
}
# now order the same
betas_both_nc <- betas_both_nc[, colnames(finemapped_colocs)]
# and merge all
overlap_complete <- rbind(finemapped_colocs, betas_both_nc)

# sort by eQTL, then caQTL strength
overlap_complete <- overlap_complete[order(abs(overlap_complete[['e_effect']]), abs(overlap_complete[['ca_effect']]), decreasing = T), ]
# now by method
overlap_complete <- overlap_complete[order(overlap_complete[['method']], decreasing = T), ]
# then keep the first entry
overlap_complete_unique <- overlap_complete[!duplicated(paste(overlap_complete[['cell_type']], overlap_complete[['hit1']], overlap_complete[['hit2']], overlap_complete[['trait1']], overlap_complete[['trait2']])), ]

# correlation of everything
sum(sign(overlap_complete_unique$ca_effect) == sign(overlap_complete_unique$e_effect)) / nrow(overlap_complete_unique)
# 0.7673073
cor(overlap_complete_unique$ca_effect, overlap_complete_unique$e_effect)
# 0.5794656

# correlation of best chromatin region and variant per gene
sum(sign(overlap_complete_unique[
  !duplicated(
    paste(
      #overlap_complete_unique[['hit1']], 
      #overlap_complete_unique[['hit2']], 
      overlap_complete_unique[['trait2']]
      )
    ), ]$ca_effect) 
  == sign(overlap_complete_unique[
    !duplicated(
      paste(
        #overlap_complete_unique[['hit1']], 
        #overlap_complete_unique[['hit2']], 
        overlap_complete_unique[['trait2']])), ]
    $e_effect)) / 
  nrow(overlap_complete_unique[
    !duplicated(
      paste(
        #overlap_complete_unique[['hit1']], 
        #overlap_complete_unique[['hit2']], 
        overlap_complete_unique[['trait2']])), ])
# 0.8087472
cor(
  overlap_complete_unique[
    !duplicated(
      paste(
        #overlap_complete_unique[['hit1']], 
        #overlap_complete_unique[['hit2']], 
        overlap_complete_unique[['trait2']])), ]$ca_effect, 
  overlap_complete_unique[
    !duplicated(
      paste(
        #overlap_complete_unique[['hit1']], 
        #overlap_complete_unique[['hit2']], 
        overlap_complete_unique[['trait2']])), ]$e_effect)
# 0.6506549

# write the result
overlap_complete_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/colocalization/eqtl_caqtl/ut_and_24hca_significant/mo_eqtl_cqtl_coloc_and_overlapping.tsv.gz'
write.table(overlap_complete, gzfile(overlap_complete_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# with a checksum
mdfiver::create_sha256_for_file(overlap_complete_loc)

# make a list to keep the significant caQTL-eQTL pairs
pairs_per_ct <- list()
genes_per_ct <- list()
regions_per_ct <- list()
overlapping_snps_per_ct <- list()
# check each of the cell types present
for (ct in unique(overlap_complete_unique[['cell_type']])) {
  # subset the rows for this cell type
  overlap_complete_unique_ct <- overlap_complete_unique[
    !is.na(overlap_complete_unique[['cell_type']]) &
      overlap_complete_unique[['cell_type']] == ct, 
  ]
  # get the unique region-gene pairs
  pairs_per_ct[[ct]] <- unique(paste(overlap_complete_unique_ct[['trait1']], overlap_complete_unique_ct[['trait2']]))
  genes_per_ct[[ct]] <- unique(paste(overlap_complete_unique_ct[['trait2']]))
  regions_per_ct[[ct]] <- unique(paste(overlap_complete_unique_ct[['trait1']]))
  # for the variants we'll only use overlapping effects
  overlap_complete_unique_ct_variant <- overlap_complete_unique_ct[
    !is.na(overlap_complete_unique_ct[['method']]) &
      overlap_complete_unique_ct[['method']] == 'overlap', 
  ]
  overlapping_snps_per_ct[[ct]] <- unique(overlap_complete_unique_ct_variant[['variant1']])
}
# plot this overlap
plot_sharing_per_celltype(pairs_per_ct, use_label_dict = T, use_color_dict = T)
plot_sharing_per_celltype(genes_per_ct, use_label_dict = T, use_color_dict = T)
plot_sharing_per_celltype(regions_per_ct, use_label_dict = T, use_color_dict = T)
plot_sharing_per_celltype(overlapping_snps_per_ct, use_label_dict = T, use_color_dict = T)


# read the file again
overlap_complete <- read.table(overlap_complete_loc, header = T, sep = '\t')
# make sure the order is the same again
overlap_complete <- overlap_complete[order(abs(overlap_complete[['e_effect']]), abs(overlap_complete[['ca_effect']]), decreasing = T), ]

# get extra annotations for the eQTLs output
strand_information_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/LimixExpAnnotationFile.incStrand.txt'
strand_information <- fread(strand_information_loc, header = T, sep = '\t')
# set names
colnames(strand_information) <- paste('gene', colnames(strand_information), sep = '_')
# add to the overlap table
overlap_complete <- cbind(overlap_complete, strand_information[match(overlap_complete[['trait2']], strand_information[['gene_feature_id']]), c('gene_chromosome', 'gene_start', 'gene_end', 'gene_strand')])
# set the TSS
overlap_complete[['TSS']] <- overlap_complete[['gene_start']]
# to the end if the strand was negative
overlap_complete[!is.na(overlap_complete[['gene_strand']]) & overlap_complete[['gene_strand']] == -1, ][['TSS']] <- overlap_complete[!is.na(overlap_complete[['gene_strand']]) & overlap_complete[['gene_strand']] == -1, ][['gene_end']]
# set to NA if strand info was NA
overlap_complete[is.na(overlap_complete[['gene_strand']]), ][['TSS']] <- NA
# get information on the QTL variants
qtl_variants_all_screen_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_screen_overlap.tsv.gz'
qtl_variants_all_screen <- fread(qtl_variants_all_screen_loc, header = T, sep = '\t')
# get the SNP position for the ATAC
overlap_complete[['atac_snp_chromosome']] <- qtl_variants_all_screen[match(overlap_complete[['hit1']], qtl_variants_all_screen[['snp_id']])][['snp_chromosome']]
overlap_complete[['atac_snp_position']] <- qtl_variants_all_screen[match(overlap_complete[['hit1']], qtl_variants_all_screen[['snp_id']])][['snp_position']]
# and for the gene
overlap_complete[['gene_snp_chromosome']] <- qtl_variants_all_screen[match(overlap_complete[['hit2']], qtl_variants_all_screen[['snp_id']])][['snp_chromosome']]
overlap_complete[['gene_snp_position']] <- qtl_variants_all_screen[match(overlap_complete[['hit2']], qtl_variants_all_screen[['snp_id']])][['snp_position']]
# calculate distance to tss
overlap_complete[['tss_dist']] <- overlap_complete[['TSS']] - overlap_complete[['gene_snp_position']]
# where the it was on the negative strand, the distance is in the other direction
overlap_complete[!is.na(overlap_complete[['gene_strand']]) & overlap_complete[['gene_strand']] == -1, ][['tss_dist']] <- -1 * overlap_complete[!is.na(overlap_complete[['gene_strand']]) & overlap_complete[['gene_strand']] == -1, ][['tss_dist']]
# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])
# set names
colnames(cpeaks_anno) <- paste('atac', colnames(cpeaks_anno), sep = '_')
# add to the overlap table
overlap_complete <- cbind(overlap_complete, cpeaks_anno[match(overlap_complete[['trait1']], cpeaks_anno[['atac_signac_hg38']]), c('atac_chr_hg38', 'atac_start_hg38', 'atac_end_hg38', 'atac_housekeeping', 'atac_screen_all')])
# get the distances
overlap_distances <- get_closest_flanks(overlap_complete, 'atac_start_hg38', 'atac_end_hg38', 'gene_start', 'gene_end')
# add that to the original table
overlap_complete[['distance']] <- overlap_distances[['min_dist']]
# make a version where we remove entries where the gene and region physically overlap
overlap_complete_distbiggerzero <- overlap_complete[overlap_complete[['distance']] > 0, ]

# keep the first entry
overlap_complete_unique_r2g <- overlap_complete[!duplicated(paste(overlap_complete[['trait1']], overlap_complete[['trait2']])), ]
# # get the concordance
# con_overlap_complete_unique_r2g <- sum(sign(overlap_complete_unique_r2g[['e_effect']]) == sign(overlap_complete_unique_r2g[['ca_effect']])) / nrow(overlap_complete_unique_r2g)
# check again, but where there is no physical ATAC and gene overlap
overlap_complete_unique_distbiggerzero_r2g <- overlap_complete_distbiggerzero[!duplicated(paste(overlap_complete_distbiggerzero[['trait1']], overlap_complete_distbiggerzero[['trait2']])), ]
# # get the concordance again
# con_overlap_complete_unique_distbiggerzero_r2g <- sum(sign(overlap_complete_unique_distbiggerzero_r2g[['e_effect']]) == sign(overlap_complete_unique_distbiggerzero_r2g[['ca_effect']])) / nrow(overlap_complete_unique_distbiggerzero_r2g)
# now just for mono
overlap_complete_unique_r2g_mono <- overlap_complete[overlap_complete[['cell_type']] == 'monocyte', ]
overlap_complete_unique_r2g_mono <- overlap_complete_unique_r2g_mono[!duplicated(paste(overlap_complete_unique_r2g_mono[['trait1']], overlap_complete_unique_r2g_mono[['trait2']])), ]
# and CD4T
overlap_complete_unique_r2g_cd4t <- overlap_complete[overlap_complete[['cell_type']] == 'CD4T', ]
overlap_complete_unique_r2g_cd4t <- overlap_complete_unique_r2g_cd4t[!duplicated(paste(overlap_complete_unique_r2g_cd4t[['trait1']], overlap_complete_unique_r2g_cd4t[['trait2']])), ]

# plot these
plot_concordanace(overlap_complete_unique_r2g)
plot_concordanace(overlap_complete_unique_distbiggerzero_r2g, main = 'Effect sizes of caQTLs versus eQTLs\n(no region/gene overlap)')
plot_concordanace(overlap_complete_unique_r2g_mono, main = 'Effect sizes of caQTLs versus eQTLs for monocytes')
plot_concordanace(overlap_complete_unique_r2g_cd4t, main = 'Effect sizes of caQTLs versus eQTLs for CD4+ T')
plot_concordanace(overlap_complete_unique_distbiggerzero_r2g[!duplicated(overlap_complete_unique_distbiggerzero_r2g[['trait2']]), ], main = 'Effect sizes of caQTLs versus eQTLs\n(no region/gene overlap, top caQTL effect per gene)')
plot_concordanace(overlap_complete_unique_r2g_mono[!duplicated(overlap_complete_unique_r2g_mono[['trait2']]), ], main = 'Effect sizes of caQTLs versus eQTLs for monocytes\n(top caQTL effect per gene)')
plot_concordanace(overlap_complete_unique_r2g_cd4t[!duplicated(overlap_complete_unique_r2g_cd4t[['trait2']]), ], main = 'Effect sizes of caQTLs versus eQTLs for CD4+ T(top caQTL effect per gene)')

# add the Z score instead
overlap_complete_unique_r2g_mono[['ca_z']] <- overlap_complete_unique_r2g_mono[['ca_effect']] / overlap_complete_unique_r2g_mono[['ca_effect_se']]
overlap_complete_unique_r2g_mono[['e_z']] <- overlap_complete_unique_r2g_mono[['e_effect']] / overlap_complete_unique_r2g_mono[['e_effect_se']]
# plot the one we'll use in the end
plot_concordanace(overlap_complete_unique_r2g_mono[!duplicated(overlap_complete_unique_r2g_mono[['trait2']]), ], main = 'Effect sizes of caQTLs versus eQTLs for monocytes\n(top caQTL effect per gene)', ca_effect_column = 'ca_z', e_effect_column = 'e_z')
# save this plot
pdf(file = '~/plots/mo_scatter_mo_caqtl_vs_eqtl_top_caqtl_per_gene_mono_z.pdf', width=5, height=5)
plot_concordanace(overlap_complete_unique_r2g_mono[!duplicated(overlap_complete_unique_r2g_mono[['trait2']]), ], main = 'Effect sizes of caQTLs versus eQTLs for monocytes\n(top caQTL effect per gene)', ca_effect_column = 'ca_z', e_effect_column = 'e_z')
dev.off()

