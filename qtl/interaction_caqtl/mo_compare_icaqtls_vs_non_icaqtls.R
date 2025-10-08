#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_compare_icaqtls_vs_non_icaqtls.R
# Function: compare the characteristics of the caQTLs that are also i-caQTLs versus those that are not
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


####################
# Main code        #
####################

# location of the eQTL output
qtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'
# get all the QTL output
qtl_output <- get_qtls_per_celltype_limix(qtl_output_loc)
# add the top effect information
qtl_output <- add_top_effect_annotation(qtl_output)
# merge all the results
qtl_output_all <- do.call('rbind', qtl_output)
# add 'chr' to the chromosome
qtl_output_all[['snp_chromosome']] <- paste0('chr', qtl_output_all[['snp_chromosome']])

# get the cpeaks overlaps for each variant
qtl_variants_all_cpeaks_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
qtl_variants_all_cpeaks <- fread(qtl_variants_all_cpeaks_loc, header = T, sep = '\t')
# add overlapping feature to QTL
qtl_output_all[['region']] <- qtl_variants_all_cpeaks[match(qtl_output_all[['snp_id']], qtl_variants_all_cpeaks[['snp_id']]), ][['overlapping_feature']]
# filter on signifiacnce
qtl_output_all_sig <- qtl_output_all[qtl_output_all[['feature_q_value']] < 0.05 &
                                       qtl_output_all[['p_value']] < qtl_output_all[['pval_nominal_threshold_global']], ]

# location of the i-eqtl output
iqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# get all the QTL output
iqtl_output <- get_qtls_per_celltype_limix(iqtl_output_loc, output_file = 'inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', gene_column='feature', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column=NULL, nominal_significance_column=NULL)
# merge all the results
iqtl_output_all <- do.call('rbind', iqtl_output)
# filter on signifiacnce
iqtl_output_all_sig <- iqtl_output_all[iqtl_output_all[['feature_q_value']] < 0.05 &
                                         iqtl_output_all[['feature_bf_eigen']] < 0.05, ]
# rename the columns
colnames(iqtl_output_all_sig) <- c('i_beta', 'i_beta_se', 'i_empirical_feature_p_value', 'i_p_value', 'snp_id', 'feature_id', 'i_n_tests_feature', 'i_feature_bf_eigen', 'i_total_bf_eigen', 'i_feature_q_value', 'i_cell_type')
# merge the iqtl to the qtl output
qtl_output_all_sig_wi <- merge(qtl_output_all_sig, iqtl_output_all_sig, on = c('snp_id', 'feature_id', 'cell_type'), all.x = T)

# get the distances
qtl_distances <- get_closest_flanks(qtl_output_all_sig_wi, 'feature_start', 'feature_end', 'snp_position', 'snp_position')
# add those distances
qtl_output_all_sig_wi[['distance']] <- qtl_distances[['min_dist']]
# add column denoting the interaction effect
qtl_output_all_sig_wi[['interaction_direction']] <- 'none'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & sign(qtl_output_all_sig_wi[['i_beta']]) == 1, ][['interaction_direction']] <- 'positive'
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['i_beta']]) & sign(qtl_output_all_sig_wi[['i_beta']]) == -1, ][['interaction_direction']] <- 'negative'

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])

# add screen annotation to the QTLs
qtl_output_all_sig_wi[['screen']] <- cpeaks_anno[match(qtl_output_all_sig_wi[['feature_id']], cpeaks_anno[['signac_hg38']]), ][['screen_all']]
# make list with tables depending on the effect strength
qtl_output_all_sig_wi_direction <- list(
  'none' = unique(qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['interaction_direction']]) & qtl_output_all_sig_wi[['interaction_direction']] == 'none', c('feature_id', 'screen', 'cell_type')]), 
     'positive' = unique(qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['interaction_direction']]) & qtl_output_all_sig_wi[['interaction_direction']] == 'positive', c('feature_id', 'screen', 'cell_type')]), 
     'negative' = unique(qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['interaction_direction']]) & qtl_output_all_sig_wi[['interaction_direction']] == 'negative', c('feature_id', 'screen', 'cell_type')]))
# get the percentage of each screen group
fracs_screen_is_ca <- get_group_proportions(qtl_output_all_sig_wi_direction, specific_trait_expression = 'CA', specific_trait_name = 'CA')
# make into a plot
p_screen_fractions_ca <- ggplot(data = fracs_screen_is_ca, mapping = aes(x = group, y = frac, fill = category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('i-caQTL direction') +
  ylab('fraction of CA screen regions') +
  ggtitle('CRE method result CA screen regions') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = roycols::get_color_list(c('is CA', 'is not CA')))
# get the percentage of each screen group
fracs_screen_is_tf <- get_group_proportions(qtl_output_all_sig_wi_direction, specific_trait_expression = 'TF', specific_trait_name = 'TF')
# make into a plot
p_screen_fractions_tf <- ggplot(data = fracs_screen_is_tf, mapping = aes(x = group, y = frac, fill = category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('i-caQTL direction') +
  ylab('fraction of TF screen regions') +
  ggtitle('CRE method result TF screen regions') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = roycols::get_color_list(c('is TF', 'is not TF')))
# for CTCF
fracs_screen_is_ctcf <- get_group_proportions(qtl_output_all_sig_wi_direction, specific_trait_expression = 'CTCF', specific_trait_name = 'CTCF')
# make into a plot
p_screen_fractions_ctcf <- ggplot(data = fracs_screen_is_ctcf, mapping = aes(x = group, y = frac, fill = category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('i-caQTL direction') +
  ylab('fraction of CTCF screen regions') +
  ggtitle('CRE method result CTCF screen regions') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = roycols::get_color_list(c('is CTCF', 'is not CTCF')))
# get the percentage of each screen group
fracs_screen_is_pls <- get_group_proportions(qtl_output_all_sig_wi_direction, specific_trait_expression = 'PLS', specific_trait_name = 'PLS')
# make into a plot
p_screen_fractions_pls <- ggplot(data = fracs_screen_is_pls, mapping = aes(x = group, y = frac, fill = category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('i-caQTL direction') +
  ylab('fraction of PLS screen regions') +
  ggtitle('CRE method result PLS screen regions') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = roycols::get_color_list(c('is PLS', 'is not PLS')))
# get the percentage of each screen group
fracs_screen_is_els <- get_group_proportions(qtl_output_all_sig_wi_direction, specific_trait_expression = 'ELS', specific_trait_name = 'ELS')
# make into a plot
p_screen_fractions_els <- ggplot(data = fracs_screen_is_els, mapping = aes(x = group, y = frac, fill = category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('i-caQTL direction') +
  ylab('fraction of ELS screen regions') +
  ggtitle('CRE method result ELS screen regions') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = roycols::get_color_list(c('is ELS', 'is not ELS')))
# get the percentage of each screen group
fracs_screen_is_pels <- get_group_proportions(qtl_output_all_sig_wi_direction, specific_trait_expression = 'pELS', specific_trait_name = 'pELS')
# make into a plot
p_screen_fractions_pels <- ggplot(data = fracs_screen_is_pels, mapping = aes(x = group, y = frac, fill = category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('i-caQTL direction') +
  ylab('fraction of pELS screen regions') +
  ggtitle('CRE method result pELS screen regions') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = roycols::get_color_list(c('is pELS', 'is not pELS')))
# get the percentage of each screen group
fracs_screen_is_dels <- get_group_proportions(qtl_output_all_sig_wi_direction, specific_trait_expression = 'dELS', specific_trait_name = 'dELS')
# make into a plot
p_screen_fractions_dels <- ggplot(data = fracs_screen_is_dels, mapping = aes(x = group, y = frac, fill = category)) + 
  geom_bar(stat = 'identity', position = 'stack') +
  xlab('i-caQTL direction') +
  ylab('fraction of dELS screen regions') +
  ggtitle('CRE method result dELS screen regions') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = roycols::get_color_list(c('is dELS', 'is not dELS')))
# show in one plot
plot_grid(
  p_screen_fractions_ca, 
  p_screen_fractions_tf, 
  p_screen_fractions_ctcf, 
  p_screen_fractions_pls, 
  p_screen_fractions_els, 
  p_screen_fractions_pels, 
  p_screen_fractions_dels, 
  nrow = 3, 
  ncol = 3
)
