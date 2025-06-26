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


####################
# Functions        #
####################


read_pseudobulk_cre_output_per_celltype <- function(pseudobulk_output_folder, cell_types=NULL, filename_output='qtl_results_all.txt.gz', significance_column='empirical_feature_p_value', significance_cutoff=0.05) {
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
    # filter the file if requested
    if (!is.null(significance_column)) {
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
  output_per_celltype <- read_pseudobulk_cre_output_per_celltype(binomial_output_folder, cell_types = cell_types, filename_output = filename_output, significance_column = significance_column, significance_cutoff = significance_cutoff)
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


####################
# Settings         #
####################

# luck seed
set.seed(7777)


####################
# Main code        #
####################

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# location of the overlapping caQTLs and eQTLs
qtl_overlap_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl_caqtl_overlap/combined/L1/all/eqtl_caqtl_overlapping_variants.tsv.gz'
# location of the pseudobulk CRE mapping
pseudobulk_output_folder <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/eQTA/L1/'
# location of the binomial method
binomial_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/combined/betas_ps/'

# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# read the qtl overlap
qtl_overlap <- fread(qtl_overlap_loc, header = T, sep = '\t')

# read the pseudobulk outputs
pseudobulk_output_ut_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, 'UT', sep = '/'))
pseudobulk_output_24hca_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, '24hCA', sep = '/'))
# merge the cell types
pseudobulk_output_ut <- do.call('rbind', pseudobulk_output_ut_list)
pseudobulk_output_24hca <- do.call('rbind', pseudobulk_output_24hca_list)
# add the condition
pseudobulk_output_ut[['condition']] <- 'UT'
pseudobulk_output_24hca[['condition']] <- '24hCA'
# merge them
pseudobulk_output <- do.call('rbind', list(pseudobulk_output_ut, pseudobulk_output_24hca))

# read the binomial results
binomial_output_list <- read_binomial_output_per_celltype(binomial_output_loc)
# merge them
binomial_output <- do.call('rbind', binomial_output_list)

# add the Z to the pseudobulk analysis
pseudobulk_output[['zscore']] <- pseudobulk_output[['beta']] / pseudobulk_output[['beta_se']]
# and the qtl overlap
qtl_overlap[['z_caqtl']] <- qtl_overlap[['beta_caqtl']] / qtl_overlap[['se_caqtl']]
qtl_overlap[['z_eqtl']] <- qtl_overlap[['beta_eqtl']] / qtl_overlap[['se_eqtl']]
# get the sign overlap
qtl_overlap[['sign']] <- sign(qtl_overlap[['z_caqtl']]) * sign(qtl_overlap[['z_eqtl']])


# sort all of them by the Z
pseudobulk_output <- pseudobulk_output[order(abs(pseudobulk_output[['zscore']])), ]
binomial_output <- binomial_output[order(abs(binomial_output[['meta_z']])), ]
scenic_output <- scenic_output[order(scenic_output[['rho_R2G']]), ]
qtl_overlap <- qtl_overlap[order(qtl_overlap[['z_eqtl']], qtl_overlap[['z_caqtl']]), ]

# get unique ones
pseudobulk_output_unique <- pseudobulk_output[!duplicated(paste(pseudobulk_output[['snp_id']], pseudobulk_output[['feature_id']])), ]
binomial_output_unique <- binomial_output[!duplicated(paste(binomial_output[['region']], binomial_output[['gene']])), ]
scenic_output_unique <- scenic_output[!duplicated(paste(scenic_output[['Region']], scenic_output[['Gene']])), ]
qtl_overlap_unique <- qtl_overlap[!duplicated(paste(qtl_overlap[['feature_caqtl']], qtl_overlap[['feature_eqtl']])), ]

# get the percentage positive
nrow(pseudobulk_output_unique[sign(pseudobulk_output_unique[['zscore']]) == 1, ]) / nrow(pseudobulk_output_unique)
# [1] 0.8557264
nrow(binomial_output_unique[sign(binomial_output_unique[['meta_z']]) == 1, ]) / nrow(binomial_output_unique)
# [1] 0.9983478
nrow(scenic_output_unique[sign(scenic_output_unique[['rho_R2G']]) == 1, ]) / nrow(scenic_output_unique)
# [1] 0.6687718
nrow(qtl_overlap_unique[sign(qtl_overlap_unique[['sign']]) == 1, ]) / nrow(qtl_overlap_unique)
# [1] 0.7630585

# add region to gene column
pseudobulk_output_unique[['r2g']] <- paste(pseudobulk_output_unique[['snp_id']], pseudobulk_output_unique[['feature_id']])
binomial_output_unique[['r2g']] <- paste(binomial_output_unique[['region']], binomial_output_unique[['gene']])
scenic_output_unique[['r2g']] <- paste(gsub(':', '-', scenic_output_unique[['Region']]), scenic_output_unique[['Gene']])
qtl_overlap_unique[['r2g']] <- paste(qtl_overlap_unique[['feature_caqtl']], qtl_overlap_unique[['feature_eqtl']])


# merge pseudobulk and binominal
pseudobulk_vs_binomial <- merge(x = pseudobulk_output_unique, y = binomial_output_unique, by = 'r2g')
# check pseudobulk and scenic
pseudobulk_vs_scenic <- merge(x = pseudobulk_output_unique, y = scenic_output_unique, by = 'r2g')
# check pseudobulk and qtls
pseudobulk_vs_qtl <- merge(x = pseudobulk_output_unique, y = qtl_overlap_unique, by = 'r2g')
# and binomial vs scenic
scenic_vs_binomial <- merge(x = scenic_output_unique, y = binomial_output_unique, by = 'r2g')
# scenic vs qtl
scenic_vs_qtl <- merge(x = scenic_output_unique, y = qtl_overlap_unique, by = 'r2g')
# binomial vs qtl
binomial_vs_qtl <- merge(x = binomial_output_unique, y = qtl_overlap_unique, by = 'r2g')


# get the concordances
nrow(pseudobulk_vs_binomial[sign(pseudobulk_vs_binomial[['zscore']]) == sign(pseudobulk_vs_binomial[['meta_z']]), ]) / nrow(pseudobulk_vs_binomial)
# [1] 0.9458333
nrow(pseudobulk_vs_scenic[sign(pseudobulk_vs_scenic[['zscore']]) == sign(pseudobulk_vs_scenic[['rho_R2G']]), ]) / nrow(pseudobulk_vs_scenic)
# [1] 0.8541311
nrow(pseudobulk_vs_qtl[sign(pseudobulk_vs_qtl[['zscore']]) == sign(pseudobulk_vs_qtl[['sign']]), ]) / nrow(pseudobulk_vs_qtl)
# [1] 0.9856528
nrow(scenic_vs_binomial[sign(scenic_vs_binomial[['meta_z']]) == sign(scenic_vs_binomial[['rho_R2G']]), ]) / nrow(scenic_vs_binomial)
# [1] 0.8156997
nrow(scenic_vs_qtl[sign(scenic_vs_qtl[['sign']]) == sign(scenic_vs_qtl[['rho_R2G']]), ]) / nrow(scenic_vs_qtl)
# [1] 0.7747748
nrow(binomial_vs_qtl[sign(binomial_vs_qtl[['meta_z']]) == sign(binomial_vs_qtl[['sign']]), ]) / nrow(binomial_vs_qtl)
# [1] 0.7746358

nrow(pseudobulk_vs_binomial)
# [1] 960
nrow(pseudobulk_vs_scenic)
# [1] 1755
nrow(pseudobulk_vs_qtl)
# [1] 1394
nrow(scenic_vs_binomial)
# [1] 586
nrow(scenic_vs_qtl)
# [1] 666
nrow(binomial_vs_qtl)
# [1] 3501

nrow(binomial_output_unique)
# [1] 12710
nrow(scenic_output_unique)
# [1] 80440
nrow(pseudobulk_output_unique)
# [1] 121935
nrow(qtl_overlap_unique)
# [1] 7677

# plot them as well
plot_grid(
  plot_concondance(pseudobulk_vs_binomial, 'zscore', 'meta_z') + ggtitle('Effects of pseudobulk vs binomial\nCRE detection') + xlab('Pseudobulk Z-score') + ylab('Binomial model Z-score'), 
  plot_concondance(pseudobulk_vs_scenic, 'zscore', 'rho_R2G') + ggtitle('Effects of pseudobulk vs SCENIC+ CRE\ndetection') + xlab('Pseudobulk Z-score') + ylab('SCENIC+ R2G Rho'), 
  plot_concondance(scenic_vs_binomial, 'rho_R2G', 'meta_z') + ggtitle('Effects of SCENIC+ vs binomial\nCRE detection') + xlab('SCENIC+ R2G Rho') + ylab('Binomial model Z-score')
)

