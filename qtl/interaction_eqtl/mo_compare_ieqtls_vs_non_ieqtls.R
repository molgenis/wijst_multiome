#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_compare_ieqtls_vs_non_ieqtls.R
# Function: compare the characteristics of the eQTLs that are also i-eQTLs versus those that are not
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
iqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
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

# plot the distances per group
p_snp_gene_interaction_distances <- ggplot(data = qtl_output_all_sig_wi, mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene') + 
  ylab('Density') + 
  ggtitle('Distance between variant and gene\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances
# the same, but only for the top variant per gene
p_snp_gene_interaction_distances_top <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and gene\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_top
# now without variant in genes
p_snp_gene_interaction_distances_distnozero <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['distance']] > 0, ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene\n(excluding variants in genes)') + 
  ylab('Density') + 
  ggtitle('Distance between variant and gene\nper interaction category (dist>0)') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_distnozero
# and top variants outside of genes
p_snp_gene_interaction_distances_distnozero_top <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['distance']] > 0 & qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene\n(excluding variants in genes)') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and gene\nper interaction category (dist>0)') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_distnozero_top
# show all together
p_snp_gene_interaction_all <- plot_grid(
  p_snp_gene_interaction_distances, 
  p_snp_gene_interaction_distances_top, 
  p_snp_gene_interaction_distances_distnozero, 
  p_snp_gene_interaction_distances_distnozero_top
)
p_snp_gene_interaction_all
# and top variants outside of genes, but also mono
p_snp_gene_interaction_distances_distnozero_top_mono <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['distance']] > 0 & qtl_output_all_sig_wi[['is_top_variant']] & qtl_output_all_sig_wi[['cell_type']] == 'monocyte', ], mapping = aes(x = distance, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and gene for\n(excluding variants in genes)') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and gene\nper interaction category (dist>0)\nfor monocytes only') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_gene_interaction_distances_distnozero_top_mono

# get tss info
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/genome_annotation.tsv'
gene_anno <- fread(gene_anno_loc, header = T, sep = '\t')
# rename columns to be the same as in limix
colnames(gene_anno) <-c('chrom', 'start', 'end', 'strand', 'gs','Transcription_Start_Site','Transcript_type')
# add TSS info to the QTL output as well
qtl_output_all_sig_wi[['Transcription_Start_Site']] <- gene_anno[match(qtl_output_all_sig_wi[['feature_id']], gene_anno[['gs']]), ][['Transcription_Start_Site']]
# get extra annotations for the pseudobulk output
strand_information_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/LimixExpAnnotationFile.incStrand.txt'
strand_information <- fread(strand_information_loc, header = T, sep = '\t')
# add to the pseudobulk
qtl_output_all_sig_wi[['strand']] <- strand_information[match(qtl_output_all_sig_wi[['feature_id']], strand_information[['feature_id']]), ][['strand']]
# set the TSS
qtl_output_all_sig_wi[['TSS']] <- qtl_output_all_sig_wi[['feature_start']]
# to the end if the strand was negative
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['TSS']] <- qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['feature_end']]
# set to NA if strand info was NA
qtl_output_all_sig_wi[is.na(qtl_output_all_sig_wi[['strand']]), ][['TSS']] <- NA
# calculate distance to tss
qtl_output_all_sig_wi[['tss_dist']] <- qtl_output_all_sig_wi[['TSS']] - qtl_output_all_sig_wi[['snp_position']]
# where the it was on the negative strand, the distance is in the other direction
qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['tss_dist']] <- -1 * qtl_output_all_sig_wi[!is.na(qtl_output_all_sig_wi[['strand']]) & qtl_output_all_sig_wi[['strand']] == -1, ][['tss_dist']]

# plot the distances per group
p_snp_tss_interaction_distances <- ggplot(data = qtl_output_all_sig_wi, mapping = aes(x = tss_dist, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between variant and TSS') + 
  ylab('Density') + 
  ggtitle('Distance between variant and TSS\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_tss_interaction_distances
# plot the distances per group for top effects
p_snp_tss_interaction_distances_top <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = tss_dist, fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Distance between top variant and TSS') + 
  ylab('Density') + 
  ggtitle('Distance between top variant and TSS\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_tss_interaction_distances_top
# plot the distances per group for top effects using absolute distances
p_snp_tss_interaction_distances_top_abs <- ggplot(data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']], ], mapping = aes(x = abs(tss_dist), fill = interaction_direction)) +
  geom_density(alpha = 0.5) +
  xlab('Absolute istance between top variant and TSS') + 
  ylab('Density') + 
  ggtitle('Absolute distance between top variant and TSS\nper interaction category') + 
  labs(fill = 'interaction\ndirection') +
  scale_fill_manual(values = list('none' = 'gray', 'positive' = '#386CB0', 'negative' = '#F0027F')) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_snp_tss_interaction_distances_top_abs
# check if there is a difference
kruskal.test(tss_dist ~ interaction_direction, data = qtl_output_all_sig_wi[qtl_output_all_sig_wi[['is_top_variant']] & !is.na(qtl_output_all_sig_wi[['tss_dist']]), ])
