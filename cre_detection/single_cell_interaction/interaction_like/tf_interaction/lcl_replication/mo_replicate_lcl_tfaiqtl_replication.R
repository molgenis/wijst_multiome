#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_replicate_lcl_tfaiqtl_replication.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)


####################
# Settings        #
####################

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# read LCL
lcl_loc <- 'mo_var_tfexpression_gene_confinement_inclcaqtls_varinregion_replicated.txt.gz'
lcl <- fread(lcl_loc, header = T, sep = '\t')
# read TFa
tfa_ieqtls_annotated_loc <-  paste('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_
donor_countrna_inclcaqtls_varinregion/merged/', 'results_fdr_with_replication.tsv.gz', sep = '/')
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
tfa_ieqtls <- tfa_ieqtls[tfa_ieqtls[['tfa:genotype_bh']] < 0.05, ]

# all
tfa_lcl_all <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa:genotype_z', 'tfa_gene_direction')], lcl[, c('variant', 'eregulon', 'feature', 'interaction_z_score_LCL', 'interaction_p_value_LCL')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
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
tfa_lcl_nominal <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa:genotype_z', 'tfa_gene_direction')], lcl[, c('variant', 'eregulon', 'feature', 'interaction_z_score_LCL')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
# show the concordance
sum(sign(tfa_lcl_nominal[['tfa:genotype_z']]) != sign(tfa_lcl_nominal[['interaction_z_score_LCL']])) / nrow(tfa_lcl_nominal)
# [1] 0.5496599

# filter with qvalue
lcl <- lcl[lcl[['interaction_q_value_LCL']] < 0.05, ]
# merge again
tfa_lcl_emperical <- merge(tfa_ieqtls[, c('variant', 'tf', 'gene', 'tfa:genotype_z', 'tfa_gene_direction')], lcl[, c('variant', 'eregulon', 'feature', 'interaction_z_score_LCL')], by.x = c('variant', 'tf', 'gene'), by.y = c('variant', 'eregulon', 'feature'))
# show the concordance
sum(sign(tfa_lcl_emperical[['tfa:genotype_z']]) != sign(tfa_lcl_emperical[['interaction_z_score_LCL']])) / nrow(tfa_lcl_emperical)
# [1] 0.8461538