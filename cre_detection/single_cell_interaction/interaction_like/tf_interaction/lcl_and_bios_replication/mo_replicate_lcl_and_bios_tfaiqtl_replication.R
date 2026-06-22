#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_replicate_lcl_and_bios_tfaiqtl_replication.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(ggplot2)


####################
# functions        #
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
# Settings        #
####################

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# read LCL
lcl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/replication/tf_interaction/sc/lcl/mo_var_tfexpression_gene_confinement_inclcaqtls_varinregion_replicated.txt.gz'
lcl <- fread(lcl_loc, header = T, sep = '\t')
# read bios
bios_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/replication/tf_interaction/sc/bios/mo_var_tfexpression_gene_confinement_inclcaqtls_varinregion_replicated.txt.gz'
bios <- fread(bios_loc, header = T, sep = '\t')
# read TFa
tfa_ieqtls_annotated_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr_with_replication.tsv.gz'
tfa_ieqtls <- fread(tfa_ieqtls_annotated_loc, header = T, sep = '\t')
# rename columns
colnames(tfa_ieqtls) <- gsub('tf', 'tfa', colnames(tfa_ieqtls))

# add explicit tf column
tfa_ieqtls[['tf']] <- gsub("_(extended|direct).*", "", tfa_ieqtls[['tfa']])
# and add a directed tfa beta
tfa_ieqtls[['tfa_beta_directed']] <- tfa_ieqtls[['tfa_beta']]
tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa_beta_directed']] <- -1 * tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa_beta']]
tfa_ieqtls[['tfa:genotype_beta_directed']] <- tfa_ieqtls[['tfa:genotype_beta']]
tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta_directed']] <- -1 * tfa_ieqtls[tfa_ieqtls[['tfa_gene_direction']] == '-/+', ][['tfa:genotype_beta']]
# add z
tfa_ieqtls[['tfa:genotype_z']] <- tfa_ieqtls[['tfa:genotype_beta']] / tfa_ieqtls[['tfa:genotype_se']]
tfa_ieqtls[['tfa:genotype_rep_z']] <- tfa_ieqtls[['tfa:genotype_beta_rep']] / tfa_ieqtls[['tfa:genotype_se_rep']]
tfa_ieqtls[['tfa:genotype_z_directed']] <- tfa_ieqtls[['tfa:genotype_beta_directed']] / tfa_ieqtls[['tfa:genotype_se']]
tfa_ieqtls[['genotype_z']] <- tfa_ieqtls[['genotype_beta']] / tfa_ieqtls[['genotype_se']]
# backup original
tfa_ieqtls_full <- tfa_ieqtls
# filter by BH
tfa_ieqtls <- tfa_ieqtls[tfa_ieqtls[['tfa:genotype_bh']] < 0.05, ]

# all
tfa_lcl_all <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa_gene_direction', 'tfa:genotype_z', 'genotype_z', 'tfa:genotype_z_directed')], lcl[, c('variant', 'eregulon', 'feature', 'interaction_z_score_LCL', 'interaction_p_value_LCL')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
tfa_lcl_all <- tfa_lcl_all[complete.cases(tfa_lcl_all), ]
# redo MTC
tfa_lcl_all[['interaction_bh_LCL']] <- p.adjust(tfa_lcl_all[['interaction_p_value_LCL']], method = 'BH')
# filter
tfa_lcl_bh <- tfa_lcl_all[tfa_lcl_all[['interaction_bh_LCL']] < 0.05, ]
# show the concordance
sum(sign(tfa_lcl_bh[['tfa:genotype_z']]) != sign(tfa_lcl_bh[['interaction_z_score_LCL']])) / nrow(tfa_lcl_bh)
# [1] 0.85

# filter LCL
lcl <- lcl[lcl[['interaction_p_value_LCL']] < 0.05, ]
# merge the two
tfa_lcl_nominal <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa_gene_direction', 'tfa:genotype_z', 'genotype_z', 'tfa:genotype_z_directed')], lcl[, c('variant', 'eregulon', 'feature', 'interaction_z_score_LCL')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
# show the concordance
sum(sign(tfa_lcl_nominal[['tfa:genotype_z']]) != sign(tfa_lcl_nominal[['interaction_z_score_LCL']])) / nrow(tfa_lcl_nominal)
# [1] 0.5496599

# filter with qvalue
lcl <- lcl[lcl[['interaction_q_value_LCL']] < 0.05, ]
# merge again
tfa_lcl_emperical <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa_gene_direction', 'tfa:genotype_z', 'genotype_z', 'tfa:genotype_z_directed')], lcl[, c('variant', 'eregulon', 'feature', 'interaction_z_score_LCL')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
# show the concordance
sum(sign(tfa_lcl_emperical[['tfa:genotype_z']]) != sign(tfa_lcl_emperical[['interaction_z_score_LCL']])) / nrow(tfa_lcl_emperical)
# [1] 0.8461538

# swap alleles
tfa_lcl_bh[['tfa:genotype_z_sameallele']] <- -1 * tfa_lcl_bh[['tfa:genotype_z']]
# plot the replication
plot_concondance(tfa_lcl_bh, d1_effect_column = 'tfa:genotype_z_sameallele', d2_effect_column = 'interaction_z_score_LCL') + xlab('TFa-i-eQTL interaction-Z') + ylab('LCL interaction-Z') + ggtitle('Concordance of TFa-i-eQTL interactions with LCL TF-i-eQTL')

# try the same as
tfa_lcl_bh[['tfa:genotype_z_sameallele']] <- -1 * tfa_lcl_bh[['tfa:genotype_z_directed']]
tfa_lcl_bh[['interaction_z_score_LCL_sameallele']] <- tfa_lcl_bh[['interaction_z_score_LCL']]
tfa_lcl_bh[tfa_lcl_bh[['tfa_gene_direction']] == '-/+', ][['interaction_z_score_LCL_sameallele']] <- -1 * tfa_lcl_bh[tfa_lcl_bh[['tfa_gene_direction']] == '-/+', ][['interaction_z_score_LCL']]
tfa_lcl_bh[['interaction_z_score_LCL_sameallele']] <- -1 * tfa_lcl_bh[['interaction_z_score_LCL_sameallele']]
plot_concondance(tfa_lcl_bh, d1_effect_column = 'tfa:genotype_z_directed', d2_effect_column = 'interaction_z_score_LCL_sameallele') + xlab('TFa-i-eQTL interaction-Z') + ylab('LCL interaction-Z') + ggtitle('Concordance of TFa-i-eQTL interactions with LCL TF-i-eQTL')


# all
tfa_bios_all <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa_gene_direction', 'tfa:genotype_z', 'genotype_z', 'tfa:genotype_z_directed')], bios[, c('variant', 'assessed_allele', 'eregulon', 'feature', 'interaction_z_score_BIOS', 'interaction_p_value_BIOS')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
tfa_bios_all <- tfa_bios_all[complete.cases(tfa_bios_all), ]
# redo MTC
tfa_bios_all[['interaction_bh_BIOS']] <- p.adjust(tfa_bios_all[['interaction_p_value_BIOS']], method = 'BH')
# filter
tfa_bios_bh <- tfa_bios_all[tfa_bios_all[['interaction_bh_BIOS']] < 0.05, ]
# show the concordance
sum(sign(tfa_bios_bh[['tfa:genotype_z_directed']]) != sign(tfa_bios_bh[['interaction_z_score_BIOS']])) / nrow(tfa_bios_bh)
# [1] 0.7192488

# filter with qvalue
bios <- bios[bios[['interaction_q_value_BIOS']] < 0.05, ]
# merge again
tfa_bios_emperical <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa_gene_direction', 'tfa:genotype_z', 'genotype_z', 'tfa:genotype_z_directed')], bios[, c('variant', 'assessed_allele', 'eregulon', 'feature', 'interaction_z_score_BIOS')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
# show the concordance
sum(sign(tfa_bios_emperical[['tfa:genotype_z_directed']]) != sign(tfa_bios_emperical[['interaction_z_score_BIOS']])) / nrow(tfa_bios_emperical)
# [1] 0.7225077

# swap alleles
tfa_bios_bh[['tfa:genotype_z_sameallele']] <- -1 * tfa_bios_bh[['tfa:genotype_z_directed']]
tfa_bios_bh[['interaction_z_score_BIOS_z_sameallele']] <- -1 * tfa_bios_bh[['interaction_z_score_BIOS']]
# plot the replication
plot_concondance(tfa_bios_bh, d1_effect_column = 'tfa:genotype_z_sameallele', d2_effect_column = 'interaction_z_score_BIOS') + xlab('TFa-i-eQTL interaction-Z') + ylab('BIOS interaction-Z') + ggtitle('Concordance of TFa-i-eQTL interactions with BIOS TF-i-eQTL')
plot_concondance(tfa_bios_bh, d1_effect_column = 'tfa:genotype_z_directed', d2_effect_column = 'interaction_z_score_BIOS_z_sameallele') + xlab('TFa-i-eQTL interaction-Z') + ylab('BIOS interaction-Z') + ggtitle('Concordance of TFa-i-eQTL interactions with BIOS TF-i-eQTL')


# get top only
tfa_bios_bh_top <- tfa_bios_bh[order(abs(tfa_bios_bh[['interaction_z_score_BIOS']])), ]
tfa_bios_bh_top <- tfa_bios_bh_top[!duplicated(tfa_bios_bh_top[, c('tf', 'gene')]), ]
# plot the replication
plot_concondance(tfa_bios_bh_top, d1_effect_column = 'tfa:genotype_z_sameallele', d2_effect_column = 'interaction_z_score_BIOS') + xlab('TFa-i-eQTL interaction-Z') + ylab('BIOS interaction-Z') + ggtitle('Concordance of TFa-i-eQTL interactions with BIOS TF-i-eQTL (top variants)')

# get top gene only
tfa_bios_bh_top_gene <- tfa_bios_bh_top[!duplicated(tfa_bios_bh_top[, c('gene')]), ]
# plot the replication
plot_concondance(tfa_bios_bh_top_gene, d1_effect_column = 'tfa:genotype_z_sameallele', d2_effect_column = 'interaction_z_score_BIOS') + xlab('TFa-i-eQTL interaction-Z') + ylab('BIOS interaction-Z') + ggtitle('Concordance of TFa-i-eQTL interactions with BIOS TF-i-eQTL (top genes)')

# get the bios overlapping set
tfa_bios_bh_formerge <- tfa_bios_bh[, c('variant', 'tf', 'gene', 'interaction_z_score_BIOS_z_sameallele', 'interaction_p_value_BIOS', 'interaction_bh_BIOS')]
# rename those
colnames(tfa_bios_bh_formerge) <- c('variant', 'tf', 'gene', 'i_z_bios', 'i_p_value_bios', 'i_bh_bios')
# merge with the tfa ieqtls
tfa_ieqtls_full <- merge(tfa_ieqtls_full, tfa_bios_bh_formerge, by = c('variant', 'tf', 'gene'), all.x = T)

# get the lcl overlapping set
tfa_lcl_bh_formerge <- tfa_lcl_bh[, c('variant', 'tf', 'gene', 'interaction_z_score_LCL_sameallele', 'interaction_p_value_LCL', 'interaction_bh_LCL')]
# rename those
colnames(tfa_lcl_bh_formerge) <- c('variant', 'tf', 'gene', 'i_z_lcl', 'i_p_value_lcl', 'i_bh_lcl')
# merge with the tfa ieqtls
tfa_ieqtls_full <- merge(tfa_ieqtls_full, tfa_lcl_bh_formerge, by = c('variant', 'tf', 'gene'), all.x = T)

# save result
tfa_with_bulk_rep_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/replication/lcl/mo_var_tfexpression_gene_confinement_inclcaqtls_varinregion_lclbios.tsv.gz'
write.table(tfa_ieqtls_full, gzfile(tfa_with_bulk_rep_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_sha256_for_file(tfa_with_bulk_rep_loc)

