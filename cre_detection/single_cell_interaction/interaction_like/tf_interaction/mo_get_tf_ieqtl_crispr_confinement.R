#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_tf_ieqtl_crispr_confinement.R
# Function: get significant TF-i-eQTLs and check if there is circularity where the eGene actually affects TF activity
############################################################################################################################

####################
# libraries        #
####################

library(data.table)


###################
# Functions        #
####################



####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- F


#######################################
# Read TF-i-eQTL original/replication #
#######################################

# location of the tf-i-ieqtls
tf_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr_with_replication.tsv.gz'
# read the inputs
tf_ieqtls <- fread(tf_ieqtl_loc, header = T, sep = '\t')
# # get what is significant
# tf_ieqtls_sig <- tf_ieqtls[
#   !is.na(tf_ieqtls[['tf:genotype_bh']]) &
#   tf_ieqtls[['tf:genotype_bh']] < 0.05 & 
#     tf_ieqtls[['anova_bh']] < 0.05 & 
#     # tf_ieqtls[['region_bh']] < 0.05 & 
#     # tf_ieqtls[['genotype_bh']] < 0.05,
#     tf_ieqtls[['tf_bh']] < 0.05, 
# ]
tf_ieqtls_sig <- tf_ieqtls[tf_ieqtls[['significant']], ]

# location of the eregulon from scenics
scenic_eregs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_signatures.tsv.gz'
# read that
scenic_eregs <- fread(scenic_eregs_loc, header = T, sep = '\t')

# we'll put this in a list
ereg_to_genes_l <- list()
# check each tf
for (eregulon in unique(tf_ieqtls_sig[['tf']])) {
  # get the associated signature genes
  eregulon_sig_genes <- scenic_eregs[
    scenic_eregs[['modality']] == 'Gene_based' & 
    scenic_eregs[['signature_name']] == eregulon, 
  ][['gene_or_region']]
  # subset the table to the i-egenes associated to this tf
  tf_ieqtls_sig_ereg <- tf_ieqtls_sig[tf_ieqtls_sig[['tf']] == eregulon, ]
  # check each of these genes
  for (gene in unique(tf_ieqtls_sig_ereg[['gene']])) {
    # get the gene list excluding this gene
    eregulon_sig_genes_no_iegene <- setdiff(eregulon_sig_genes, gene)
    # make into a dataframe
    eregulon_sig_genes_no_iegene_df <- data.table(
      'eregulon' = rep(eregulon, times = length(eregulon_sig_genes_no_iegene)),
      'egene' = rep(gene, times = length(eregulon_sig_genes_no_iegene)),
      'other_gene' = eregulon_sig_genes_no_iegene
    )
    # put in the list
    ereg_to_genes_l[[paste(eregulon, gene)]] <- eregulon_sig_genes_no_iegene_df
  }
}
# put all together
ereg_to_genes_dt <- rbindlist(ereg_to_genes_l)

# store this result somewhere
ereg_to_genes_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/replication/tf_interaction/sc/mo_ereg_no_iegene.tsv.gz'
write.table(ereg_to_genes_dt, gzfile(ereg_to_genes_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# make checksum
mdfiver::create_sha256_for_file(ereg_to_genes_loc)
