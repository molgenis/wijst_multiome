#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_compare_tfieqtls_vs_tfaieqtls.R
# Function: compare TF-i-eQTLs coming from TF expression versus ones coming from TF activity
############################################################################################################################


####################
# libraries        #
####################

# for plots
library(ggplot2)
# data table
library(data.table)


####################
# Functions        #
####################


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
    annotate("text", x = max_sig_z_d1 * -0.70 , y = max_sig_z_d2 * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold') +
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
# whether we are in debug mode
debug <- T


###################
# interested data #
###################



# set location of both tables
tfa_ieqtls_annotated_loc <- paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/', 'results_fdr_with_replication.tsv.gz', sep = '/')
tfe_ieqtls_loc <- paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion_tfexpressions/merged/', 'results_fdr.tsv.gz', sep = '/')
# read these
tfa_ieqtls <- fread(tfa_ieqtls_annotated_loc, header = T, sep = '\t')
tfe_ieqtls <- fread(tfe_ieqtls_loc, header = T, sep = '\t')
# rename columns
colnames(tfa_ieqtls) <- gsub('tf', 'tfa', colnames(tfa_ieqtls))
colnames(tfe_ieqtls) <- gsub('region', 'tfe', colnames(tfe_ieqtls))
# add explicit tf column
tfa_ieqtls[['tf']] <- gsub("_(extended|direct).*", "", tfa_ieqtls[['tfa']])
# which will just be the same for the tfe
tfe_ieqtls[['tf']] <- tfe_ieqtls[['tfe']]
# and add a directed tfa beta
tfa_ieqtls[['tfa_beta_directed']] <- tfa_ieqtls[['tfa_beta']]
tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa_beta_directed']] <- -1 * tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa_beta']]
tfa_ieqtls[['tfa:genotype_beta_directed']] <- tfa_ieqtls[['tfa:genotype_beta']]
tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta_directed']] <- -1 * tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta']]
# add z
tfa_ieqtls[['tfa:genotype_z']] <- tfa_ieqtls[['tfa:genotype_beta']] / tfa_ieqtls[['tfa:genotype_se']]
tfa_ieqtls[['tfa:genotype_rep_z']] <- tfa_ieqtls[['tfa:genotype_beta_rep']] / tfa_ieqtls[['tfa:genotype_se_rep']]
tfa_ieqtls[['tfa:genotype_z_directed']] <- tfa_ieqtls[['tfa:genotype_beta_directed']] / tfa_ieqtls[['tfa:genotype_se']]
tfe_ieqtls[['tfe:genotype_z']] <- tfe_ieqtls[['tfe:genotype_beta']] / tfe_ieqtls[['tfe:genotype_se']]

# merge the tf and region qtls
tfa_tfe_ieqtls_overlap <- merge(
  tfa_ieqtls, 
  tfe_ieqtls, 
  by.x = c('variant', 'tf', 'gene'), 
  by.y = c('variant', 'tf', 'gene')
)
# get what is significant in tfa
tfa_tfe_ieqtls_overlap_tfasig <- tfa_tfe_ieqtls_overlap[tfa_tfe_ieqtls_overlap[['significant']], ]
# then redo mtc for replicating
tfa_tfe_ieqtls_overlap_tfasig[['tfe:genotype_bh']] <- p.adjust(tfa_tfe_ieqtls_overlap_tfasig[['tfe:genotype_p']], method = 'BH')
tfa_tfe_ieqtls_overlap_tfasig[['anova_bh.y']] <- p.adjust(tfa_tfe_ieqtls_overlap_tfasig[['anova_p.y']], method = 'BH')
# show concordance
p_tfa_vs_tfe_ieqtls_all <- plot_concondance(tfa_tfe_ieqtls_overlap_tfasig, d1_effect_column = 'tfa:genotype_rep_z', d2_effect_column = 'tfe:genotype_z') + 
  xlab('TFa-i-eQTL interaction Z-score') +
  ylab('TFe-i-eQTL interaction Z-score')
# show the plot
p_tfa_vs_tfe_ieqtls_all
# save plot
ggsave(file = '~/multiome/plots/mo_tfa_vs_tfe_ieqtls_all.pdf', plot = p_tfa_vs_tfe_ieqtls_all, width = 5, height = 5)

# show concordance for nominal significance in the TFe-i-eQTLs
p_tfa_vs_tfe_ieqtls_nominal <- plot_concondance(tfa_tfe_ieqtls_overlap_tfasig[tfa_tfe_ieqtls_overlap_tfasig[['anova_p.y']] < 0.05 &
                                                 tfa_tfe_ieqtls_overlap_tfasig[['tfe:genotype_p']] < 0.05,
                                               ], d1_effect_column = 'tfa:genotype_rep_z', d2_effect_column = 'tfe:genotype_z') + 
  xlab('TFa-i-eQTL interaction Z-score') +
  ylab('TFe-i-eQTL interaction Z-score')
# show the plot
p_tfa_vs_tfe_ieqtls_nominal
# save plot
ggsave(file = '~/multiome/plots/mo_tfa_vs_tfe_ieqtls_nominaltfe.pdf', plot = p_tfa_vs_tfe_ieqtls_nominal, width = 5, height = 5)
