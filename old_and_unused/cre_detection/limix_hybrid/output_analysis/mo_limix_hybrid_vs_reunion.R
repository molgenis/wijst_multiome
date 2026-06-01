#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_limix_hybrid_vs_reunion.R
# Function: compare LIMIX hybrid output to REUNION output
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(ggplot2)
library(ggvenn)
library(cowplot)
library(stringr)
library(IRanges)


####################
# Functions        #
####################


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
  return(output_per_celltype)
}



###################
# Settings        #
###################


####################
# Main Code        #
####################

# location of annotated reunion output
reunion_output_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/reunion_yang2024_cres/df_pbmc_region_tf_gene_link.2_2_cpeaksmatched.tsv.gz'
# location of the hybrid method
hybrid_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/'

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])

# read the location of the genes
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno <- fread(gene_anno_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')

# read hybrid method
# hybrid_output_list <- read_pseudobulk_cre_output_per_celltype(hybrid_output_loc, add_mtc = F, filter_alpha = F, add_global_nominal_threshold = F, add_local_nominal_threshold = F, filename_output = 'qtl_results_all_frac01.txt.gz', alpha_min = .8, alpha_max = 1.2, filter_significance = F)
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
# add annotation of peak location for peaks
hybrid_output <- cbind(hybrid_output, cpeaks_anno[match(hybrid_output[['snp_id']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# get distances between peak and gene
hybrid_distances <- get_closest_flanks(hybrid_output, 'start_hg38', 'end_hg38', 'feature_start', 'feature_end')
# add this to the table
hybrid_output[['distance']] <- hybrid_distances[['min_dist']]
# remove regions with larger p values than our cutoff
# hybrid_output <- hybrid_output[hybrid_output[['p_value']] < hybrid_output[['pval_nominal_threshold_global']], ]
hybrid_output <- hybrid_output[hybrid_output[['global_significance']] == T, ]

# filter based on distance
hybrid_output_unfiltered <- hybrid_output
hybrid_output <- hybrid_output[hybrid_output[['distance']] > 0 & hybrid_output[['distance']] <= 150000, ]

# load reunion
reunion_output <- fread(reunion_output_loc, header = T, sep = '\t')
# make the overlapping feature into the signac format
reunion_output[['overlapping_feature']] <- gsub(':', '-', reunion_output[['overlapping_feature']])

# same filtering for reunion
reunion_output_unfiltered <- reunion_output
reunion_output <- reunion_output[reunion_output[['distance']] > 0 & reunion_output[['distance']] <= 150000, ]

# now look at the overlap
ggvenn::ggvenn(
  data = list('LMM' = unique(hybrid_output[['snp_id']]), 'REUNION' = unique(reunion_output[['overlapping_feature']]))
) + ggtitle('Overlap of CREs between LMM and REUNION')
ggvenn::ggvenn(
  data = list('LMM' = unique(hybrid_output[['feature_id']]), 'REUNION' = unique(reunion_output[['gene_id']]))
) + ggtitle('Overlap of CRE genes between LMM and REUNION')
ggvenn::ggvenn(
  data = list('LMM' = unique(paste(hybrid_output[['snp_id']], hybrid_output[['feature_id']])), 'REUNION' = unique(paste(reunion_output[['overlapping_feature']], reunion_output[['gene_id']])))
) + ggtitle('Overlap of region-gene links between LMM and REUNION')

# sort the hybrid output by significance
hybrid_output <- hybrid_output[order(hybrid_output[['p_value']]), ]

# add the r2g column
hybrid_output[['r2g']] <- paste(hybrid_output[['snp_id']], hybrid_output[['feature_id']])
reunion_output[['r2g']] <- paste(reunion_output[['overlapping_feature']], reunion_output[['gene_id']])
# and keep top results
hybrid_output_top <- hybrid_output[!duplicated(hybrid_output[['r2g']]), ]
reunion_output_top <- reunion_output[!duplicated(reunion_output[['r2g']]), ]

# merge these
matched_output <- merge(reunion_output_top[, c('r2g', 'peak_gene_corr')], hybrid_output_top[, c('r2g', 'r', 'z_score')], by = 'r2g')
# rename columns
colnames(matched_output) <- c('r2g', 'cor_reunion', 'cor_hybrid', 'z_hybrid')
# and make unique
matched_output <- unique(matched_output)

# calculate the concordance
mo_reunion_rho_concordance <- sum(sign(matched_output[['z_hybrid']]) == sign(matched_output[['cor_reunion']])) / nrow(matched_output)
# and the replication rate
mo_reunion_replication <- sum(sign(matched_output[['z_hybrid']]) == sign(matched_output[['cor_reunion']])) / nrow(hybrid_output_top)
# get the minimal and maximum correlations
min_sig_rho_10x <- min(abs(matched_output[['cor_reunion']]))
min_sig_rho_mo <- min(abs(matched_output[['z_hybrid']]))
max_sig_rho_10x <- max(abs(matched_output[['cor_reunion']]))
max_sig_rho_mo <- max(abs(matched_output[['z_hybrid']]))
# plot the correlations
ggplot(data = matched_output, mapping = aes(x = z_hybrid, y = cor_reunion)) + 
  geom_point() +
  xlab('R2G Z in LIMIX LMM') + 
  ylab('R2G Rho in reunion') + 
  ggtitle('R2G correlations in LIMIX LMM\n vs reunion') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # left to right block of non-significant effects
  geom_rect(aes(xmin = -1 * max_sig_rho_mo, xmax = max_sig_rho_mo, ymin = -1 * min_sig_rho_10x, ymax = min_sig_rho_10x), 
            fill = "white", alpha = 0.005) +
  # bottom to top block of non-significant effects
  geom_rect(aes(xmin = -1 * min_sig_rho_mo, xmax = min_sig_rho_mo, ymin = -1 * max_sig_rho_10x, ymax = max_sig_rho_10x), 
            fill = "white", alpha = 0.005) +
  # bottom left block
  geom_rect(aes(xmin = -1 * max_sig_rho_mo, xmax = -1 *min_sig_rho_mo, ymin = -1 * max_sig_rho_10x, ymax = -1 * min_sig_rho_10x), 
            fill = "#0072B2", alpha = 0.01) +
  # bottom right block
  geom_rect(aes(xmin = min_sig_rho_mo, xmax = max_sig_rho_mo, ymin = -1 * max_sig_rho_10x, ymax = -1 * min_sig_rho_10x), 
            fill = "#D55E00", alpha = 0.01) +
  # top left block
  geom_rect(aes(xmin = -1 * max_sig_rho_mo, xmax = -1 *min_sig_rho_mo, ymin = min_sig_rho_10x, ymax = max_sig_rho_10x), 
            fill = "#D55E00", alpha = 0.01) + 
  # top right block
  geom_rect(aes(xmin = min_sig_rho_mo, xmax = max_sig_rho_mo, ymin = max_sig_rho_10x, ymax = min_sig_rho_10x), 
            fill = "#0072B2", alpha = 0.01) + 
  # add the concordance
  annotate("label", x = max_sig_rho_mo * 0.75 , y = max_sig_rho_10x * -0.75, label = paste('concordance', round(mo_reunion_rho_concordance, digits = 2), sep = ':\n')) +
  # add the replication
  annotate("label", x = max_sig_rho_mo * -0.75 , y = max_sig_rho_10x * -0.75, label = paste('replication LMM in Reunion', round(mo_reunion_replication, digits = 2), sep = ':\n')) +
  # add the names of the concordant and non-concordant blocks
  annotate("text", x = max_sig_rho_mo * -0.70 , y = max_sig_rho_10x * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold') +
  # add the names of the concordant and non-concordant blocks
  annotate("text", x = max_sig_rho_mo * 0.70 , y = max_sig_rho_10x * 0.75, label = 'concordant', colour = '#0072B2', fontface = 'bold')

