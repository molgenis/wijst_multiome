#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_hybrid_cre_interaction_overlaps.R
# Function: plot interaction-eQTL at single-cell level with TF or ATAC as interaction terms overlaps across methods
# Example: 
# 
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)
# to plot
library(ggplot2)
library(roycols)
library(cowplot)


####################
# Functions        #
####################


get_tf_ieqtl_outputs <- function(output_loc, cell_types=NULL, output_file='merged/results_fdr.tsv.gz', add_fdr=F) {
  # we'll store the results in a list first
  res_per_ct <- list()
  # we'll list the directories in the output loc
  output_dirs <- list.dirs(output_loc, recursive = F)
  # get the base dirs
  base_dirs <- basename(output_dirs)
  # we'll intersect with the cell types if we have that information
  if (!is.null(cell_types)) {
    output_dirs <- output_dirs[base_dirs %in% cell_types]
    base_dirs <- base_dirs[base_dirs %in% cell_types]
  }
  # we'll check each cell type
  for (path_i in 1 : length(output_dirs)) {
    # get base name
    base_folder <- base_dirs[path_i]
    output_folder <- output_dirs[path_i]
    # get the full path
    full_path <- paste(output_folder, output_file, sep = '/')
    # load the file
    output <- fread(full_path, header = T, sep = '\t')
    # add fdr if needed
    if (add_fdr) {
      # extract each column that has a p value
      p_columns <- colnames(output)[grep('_p$', colnames(output))]
      # do MTC for each of these columns
      for (p_column in p_columns) {
        # check which rows have a value for this p
        p_valid_i <- !is.na(output[[p_column]])
        # add MTC column
        bh_column <- gsub('_p$', '_bh', p_column)
        output[[bh_column]] <- NA
        # extract valid p values, and calculate B&H
        bhs <- p.adjust(output[[p_column]][p_valid_i], method = 'BH')
        # then place those BHs
        output[[bh_column]][p_valid_i] <- bhs
        # add bf column as well
        bf_column <- gsub('_p$', '_bh', p_column)
        # adjust the valid p values
        bfs <- p.adjust(output[[p_column]][p_valid_i], method = 'bonferroni')
        # add the bonferroni column
        output[[bf_column]] <- NA
        # and place the bonferroni adjusted p values
        output[[bf_column]][p_valid_i] <- bfs
      }
    }
    # add the cell type
    output <- cbind(data.table('cell_type' = rep(base_folder, times = nrow(output))), output)
    # put in the list
    res_per_ct[[base_folder]] <- output
  }
  # merge all
  res_all <- rbindlist(res_per_ct, fill = T)
  return(res_all)
}


plot_concondance <- function(dataset_to_compare, d1_effect_column='d1_zscore', d2_effect_column='d2_zscore', concordance_number_font_size=8, concordance_number_label_size=8) {
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
    annotate("label", x = max_sig_z_d1 * 0.75 , y = max_sig_z_d2 * -0.75, label = paste('concordance', concordance, sep = ':\n'), size = concordance_number_font_size) +
    # add the names of the concordant and non-concordant blocks
    annotate("text", x = max_sig_z_d1 * -0.70 , y = max_sig_z_d2 * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold', size = concordance_number_label_size) +
    # add the names of the concordant and non-concordant blocks
    annotate("text", x = max_sig_z_d1 * 0.70 , y = max_sig_z_d2 * 0.75, label = 'concordant', colour = '#0072B2', fontface = 'bold', size = concordance_number_label_size) +
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

# the location of SCENIC output
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# order by the extended
scenic_output <- scenic_output[order(scenic_output[['is_extended']]), ]
# remove what SCENIC thinks is less likely
scenic_output <- scenic_output[scenic_output[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# and keep only what is non-extended if it was both extended and non-extended
scenic_output <- scenic_output[!duplicated(paste(scenic_output[['Region']], scenic_output[['Gene']], scenic_output[['Gene_signature_direction']])), ]
# make the regions cpeak style
scenic_output[['region_cpeaks']] <- gsub(':', '-', scenic_output[['Region']])

# results of the pseudobulk method
ps_tf_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/pseudobulked/'
# result of the single-cell interaction method
sc_tf_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/'
# results of the co-eqtl method
co_tf_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/coeqtl/'

# paths to the files
ps_tf_ieqtl_file <- '/result.tsv.gz'
sc_tf_ieqtl_file <- '/lane_donor_countrna_inclcaqtls/merged/results_fdr.tsv.gz'
co_tf_ieqtl_file <- '/lane_donor_inclcaqtls/merged/results_fdr.tsv.gz'

# get all the results
ps_tf_ieqtls <- get_tf_ieqtl_outputs(ps_tf_ieqtl_loc, output_file = ps_tf_ieqtl_file, add_fdr = T)
sc_tf_ieqtls <- get_tf_ieqtl_outputs(sc_tf_ieqtl_loc, output_file = sc_tf_ieqtl_file, cell_types = c('all'))
co_tf_ieqtls <- get_tf_ieqtl_outputs(co_tf_ieqtl_loc, output_file = co_tf_ieqtl_file, cell_types = c('all'))

# get what is significant
sc_tf_ieqtls_sig <- sc_tf_ieqtls[
    sc_tf_ieqtls[['region:genotype_bh']] < 0.05 & 
    sc_tf_ieqtls[['anova_bh']] < 0.05 & 
    sc_tf_ieqtls[['region_bh']] < 0.05 & 
    sc_tf_ieqtls[['genotype_bh']] < 0.05, 
]
co_tf_ieqtls <- co_tf_ieqtls[
    co_tf_ieqtls[['genotype_bh']] < 0.05, 
]
# merge what is significant in both
sc_and_co_tf_ieqtls_sig <- merge(sc_tf_ieqtls_sig, co_tf_ieqtls, by.x = c('cell_type', 'variant', 'region', 'gene'), by.y = c('cell_type', 'variant', 'region', 'gene'))
# check concordance
sc_and_co_tf_ieqtls_sig_cond <- sum(sign(sc_and_co_tf_ieqtls_sig[['region:genotype_beta']]) == sign(sc_and_co_tf_ieqtls_sig[['genotype_beta.y']])) / nrow(sc_tf_ieqtls_sig)

# merge the sc tf-i-eQTLs with the scenic output
scenic_sc_tf_ieqtls <- merge(scenic_output, sc_tf_ieqtls, by.x = c('Gene_signature_name', 'Gene'), by.y = c('region', 'gene'), allow.cartesian = T)
# check the ieQTL direction and the signs
scenic_sc_tf_ieqtls[['ieqtl_sign']] <- sign(scenic_sc_tf_ieqtls[['region:genotype_beta']])
# and the TF-gene direction
scenic_sc_tf_ieqtls[['tf2g_sign']] <- sign(scenic_sc_tf_ieqtls[['rho_TF2G']])
# add the Z-score for the accessiblity i-eQTLs as well
scenic_sc_tf_ieqtls[['region_z']] <- scenic_sc_tf_ieqtls[['region_beta']] / scenic_sc_tf_ieqtls[['region_se']]
# get what is significant
scenic_sc_tf_ieqtls_sig <- scenic_sc_tf_ieqtls[
  scenic_sc_tf_ieqtls[['region:genotype_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls[['anova_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls[['region_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls[['genotype_bh']] < 0.05, 
]
# make into table
scenic_sc_tf_ieqtls_occ <- data.frame(table(scenic_sc_tf_ieqtls_sig[, c('tf2g_sign', 'ieqtl_sign')]))
# plot this 
scenic_sc_tf_ieqtls_occ[['ieqtl_sign_string']] <- ifelse(scenic_sc_tf_ieqtls_occ[['ieqtl_sign']] == -1, 'TF weakens QTL', 'TF strengthens QTL')
scenic_sc_tf_ieqtls_occ[['tf2g_sign_string']] <- ifelse(scenic_sc_tf_ieqtls_occ[['tf2g_sign']] == -1, 'TF represses gene', 'TF increases gene')
# make a plot
p_scenic_sc_tf_ieqtls_occ <- ggplot(scenic_sc_tf_ieqtls_occ, aes(x = ieqtl_sign_string, y = Freq, fill = tf2g_sign_string)) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  xlab('TF QTL directions') + 
  ylab('Number of significant TF-i-eQTLs') + 
  ggtitle('Direction of significant i-eQTLs per SCENIC+ gene signature direction') + 
  theme(
    legend.title = element_text(size=16), 
    legend.text = element_text(size=14), 
    axis.title.x = element_text(size=16), 
    axis.title.y = element_text(size=16), 
    axis.text.y = element_text(size=14), 
    axis.text.x = element_text(size=14), 
    strip.text.x = element_text(size=14)
  ) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  scale_fill_manual(values=list('TF increases gene' = 'darkred', 'TF represses gene' = 'darkblue'), name = 'SCENIC+ TF direction')

# also try to check the concordance of the TF-gene relations
p_scenic_vs_tf_ieqtls <- plot_concondance(scenic_sc_tf_ieqtls, 'rho_TF2G', 'region_beta') + 
  # with correct labels
  xlab(paste('SCENIC+ TF2G\n(Rho)')) + 
  ylab(paste('Expression ~ TF\n(beta)')) + 
  # and the title
  ggtitle(paste('Effect size of SCENIC+ vs TF-i-eQTLs\n(including non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# use just one region to gene direction
scenic_sc_tf_ieqtls_top_tf2g <- scenic_sc_tf_ieqtls[order(abs(scenic_sc_tf_ieqtls[['rho_TF2G']])), ]
# by keeping the top effect
scenic_sc_tf_ieqtls_top_tf2g <- scenic_sc_tf_ieqtls_top_tf2g[!duplicated(paste(scenic_sc_tf_ieqtls_top_tf2g[['Gene_signature_name']], scenic_sc_tf_ieqtls_top_tf2g[['Gene']])), ]
# also try to check the concordance of the TF-gene relations
p_scenic_vs_tf_ieqtls_top <- plot_concondance(scenic_sc_tf_ieqtls_top_tf2g, 'rho_TF2G', 'region_beta') + 
  # with correct labels
  xlab(paste('SCENIC+ TF2G\n(Rho)')) + 
  ylab(paste('Expression ~ TF\n(beta)')) + 
  # and the title
  ggtitle(paste('Effect size of top SCENIC+ vs TF-i-eQTLs\n(including non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# and with Z
p_scenic_vs_tf_ieqtls_top_z <- plot_concondance(scenic_sc_tf_ieqtls_top_tf2g, 'rho_TF2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ TF2G\n(Rho)')) + 
  ylab(paste('Expression ~ TF\n(Z-score)')) + 
  # and the title
  ggtitle(paste('Effect size of top SCENIC+ vs TF-i-eQTLs\n(including non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# get what is significant for the top
scenic_sc_tf_ieqtls_top_tf2g_sig <- scenic_sc_tf_ieqtls_top_tf2g[
  scenic_sc_tf_ieqtls_top_tf2g[['region:genotype_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls_top_tf2g[['anova_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls_top_tf2g[['region_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls_top_tf2g[['genotype_bh']] < 0.05, 
]
p_scenic_vs_tf_ieqtls_top_z_sig <- plot_concondance(scenic_sc_tf_ieqtls_top_tf2g_sig, 'rho_TF2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ TF2G\n(Rho)')) + 
  ylab(paste('Expression ~ TF\n(Z-score)')) + 
  # and the title
  ggtitle(paste('Effect size of top SCENIC+ vs TF-i-eQTLs\n(only significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# get what is not significant for the top
scenic_sc_tf_ieqtls_top_tf2g_nonsig <- scenic_sc_tf_ieqtls_top_tf2g[!(
  scenic_sc_tf_ieqtls_top_tf2g[['region:genotype_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls_top_tf2g[['anova_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls_top_tf2g[['region_bh']] < 0.05 & 
    scenic_sc_tf_ieqtls_top_tf2g[['genotype_bh']] < 0.05), 
]
p_scenic_vs_tf_ieqtls_top_z_nonsig <- plot_concondance(scenic_sc_tf_ieqtls_top_tf2g_nonsig, 'rho_TF2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ TF2G\n(Rho)')) + 
  ylab(paste('Expression ~ TF\n(Z-score)')) + 
  # and the title
  ggtitle(paste('Effect size of top SCENIC+ vs TF-i-eQTLs\n(only non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# plot the TF plots
p_scenic_vs_tf_ieqtls_top_z_grid <- plot_grid(
  p_scenic_vs_tf_ieqtls_top_z, 
  p_scenic_vs_tf_ieqtls_top_z_sig, 
  p_scenic_vs_tf_ieqtls_top_z_nonsig, 
  nrow = 2, 
  ncol = 2
)
# get what is significant for just region to TF
scenic_sc_tf_ieqtls_top_tf2g_onlytfsig <- scenic_sc_tf_ieqtls_top_tf2g[
  scenic_sc_tf_ieqtls_top_tf2g[['region_bh']] < 0.05,  
]
p_scenic_vs_tf_ieqtls_top_z_onlytfsig <- plot_concondance(scenic_sc_tf_ieqtls_top_tf2g_onlytfsig, 'rho_TF2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ TF2G\n(Rho)')) + 
  ylab(paste('Expression ~ TF\n(Z-score)')) + 
  # and the title
  ggtitle(paste('Effect size of top SCENIC+ vs TF-i-eQTLs\n(only TF significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# save the plot
ggsave('~/multiome/plots/mo_scenic_hybrid_tf_ieqtl_concordance.pdf', plot = p_scenic_vs_tf_ieqtls_top_z_grid, width = 15, height = 15)

# merge the co tf-i-eQTLs with the scenic output
scenic_co_tf_ieqtls <- merge(scenic_output, co_tf_ieqtls, by.x = c('Gene_signature_name', 'Gene'), by.y = c('region', 'gene'), allow.cartesian = T)
# check the region-signature direction and the signs
scenic_co_tf_ieqtls[['ieqtl_sign']] <- sign(scenic_co_tf_ieqtls[['genotype_beta']])
scenic_co_tf_ieqtls_occ <- data.frame(table(scenic_co_tf_ieqtls[, c('Gene_signature_direction', 'ieqtl_sign')]))

# read the region interaction as well
sc_region_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/region_interaction/sc/'
sc_region_ieqtl_file <- '/lane_donor_countrna_inclcaqtls/merged/results_fdr.tsv.gz'
sc_region_ieqtls <- get_tf_ieqtl_outputs(sc_region_ieqtl_loc, output_file = sc_region_ieqtl_file, cell_types = c('all'), add_fdr = F)
# get what is significant
sc_region_ieqtls_sig <- sc_region_ieqtls[
  sc_region_ieqtls[['region:genotype_bh']] < 0.05 & 
    sc_region_ieqtls[['anova_bh']] < 0.05, 
]
# merge the sc region-i-eQTLs with the scenic output
scenic_sc_region_ieqtls <- merge(scenic_output, sc_region_ieqtls, by.x = c('region_cpeaks', 'Gene'), by.y = c('region', 'gene'), allow.cartesian = T)
# check concordance
scenic_sc_region_ieqtls_cond <- sum(sign(scenic_sc_region_ieqtls[['rho_R2G']]) == sign(scenic_sc_region_ieqtls[['region_beta']])) / nrow(scenic_sc_region_ieqtls)
# plot both of them
p_scenic_vs_region_ieqtls <- plot_concondance(scenic_sc_region_ieqtls, 'rho_R2G', 'region_beta') + 
  # with correct labels
  xlab(paste('SCENIC+ R2G Rho')) + 
  ylab(paste('Expression ~ Accessibility beta')) + 
  # and the title
  ggtitle(paste('Effect size of SCENIC+ vs accessibility-i-eQTLs\n(including non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# add the Z-score for the accessiblity i-eQTLs as well
scenic_sc_region_ieqtls[['region_z']] <- scenic_sc_region_ieqtls[['region_beta']] / scenic_sc_region_ieqtls[['region_se']]
# plot both of them
p_scenic_vs_region_ieqtls_z <- plot_concondance(scenic_sc_region_ieqtls, 'rho_R2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ R2G Rho')) + 
  ylab(paste('Expression ~ Accessibility Z-score')) + 
  # and the title
  ggtitle(paste('Effect size of SCENIC+ vs accessibility-i-eQTLs')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# use just one region to gene direction
scenic_sc_region_ieqtls_top_r2g <- scenic_sc_region_ieqtls[order(abs(scenic_sc_region_ieqtls[['rho_R2G']])), ]
# by keeping the top effect
scenic_sc_region_ieqtls_top_r2g <- scenic_sc_region_ieqtls_top_r2g[!duplicated(paste(scenic_sc_region_ieqtls_top_r2g[['Region']], scenic_sc_region_ieqtls_top_r2g[['Gene']])), ]
# plot both of them
p_scenic_vs_region_ieqtls_top <- plot_concondance(scenic_sc_region_ieqtls_top_r2g, 'rho_R2G', 'region_beta') + 
  # with correct labels
  xlab(paste('SCENIC+ R2G Rho')) + 
  ylab(paste('Expression ~ Accessibility beta')) + 
  # and the title
  ggtitle(paste('Effect size of SCENIC+ vs accessibility-i-eQTLs\n(including non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# and Z
p_scenic_vs_region_ieqtls_top_z <- plot_concondance(scenic_sc_region_ieqtls_top_r2g, 'rho_R2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ R2G Rho')) + 
  ylab(paste('Expression ~ Accessibility Z-score')) + 
  # and the title
  ggtitle(paste('Effect size of SCENIC+ vs accessibility-i-eQTLs\n(including non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# get what is not significant for the top
scenic_sc_region_ieqtls_top_r2g_nonsig <- scenic_sc_region_ieqtls_top_r2g[!(
  scenic_sc_region_ieqtls_top_r2g[['region:genotype_bh']] < 0.05 & 
    scenic_sc_region_ieqtls_top_r2g[['anova_bh']] < 0.05), 
]
p_scenic_vs_region_ieqtls_top_z_nonsig <- plot_concondance(scenic_sc_region_ieqtls_top_r2g_nonsig, 'rho_R2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ R2G Rho')) + 
  ylab(paste('Expression ~ Accessibility Z-score')) + 
  # and the title
  ggtitle(paste('Effect size of SCENIC+ vs accessibility-i-eQTLs\n(only non-significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# get what is significant for the top
scenic_sc_region_ieqtls_top_r2g_sig <- scenic_sc_region_ieqtls_top_r2g[
  scenic_sc_region_ieqtls_top_r2g[['region:genotype_bh']] < 0.05 & 
    scenic_sc_region_ieqtls_top_r2g[['anova_bh']] < 0.05, 
]
p_scenic_vs_region_ieqtls_top_z_sig <- plot_concondance(scenic_sc_region_ieqtls_top_r2g_sig, 'rho_R2G', 'region_z') + 
  # with correct labels
  xlab(paste('SCENIC+ R2G Rho')) + 
  ylab(paste('Expression ~ Accessibility Z-score')) + 
  # and the title
  ggtitle(paste('Effect size of SCENIC+ vs accessibility-i-eQTLs\n(only significant)')) + 
  # with bigger labels
  theme(plot.title = element_text(size=18), 
        legend.title = element_text(size=16), 
        legend.text = element_text(size=14),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14),
        strip.text.x = element_text(size=14)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
# plot the accessibility plots
p_scenic_vs_region_ieqtls_top_z_grid <- plot_grid(
  p_scenic_vs_region_ieqtls_top_z, 
  p_scenic_vs_region_ieqtls_top_z_sig, 
  # p_scenic_vs_region_ieqtls_top_z_nonsig, 
  nrow = 2, 
  ncol = 2
)
# save the plot
ggsave('~/multiome/plots/mo_scenic_hybrid_region_ieqtl_concordance.pdf', plot = p_scenic_vs_region_ieqtls_top_z_grid, width = 15, height = 15)

