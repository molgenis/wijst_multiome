#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_sccre_coeqtl_tf_interaction_confinement.R
# Function: create variant-TF-gene confinements for TF interaction-eQTL mapping trying to replicate co-eQTLs
# 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(stringr)


####################
# Functions        #
####################



####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- T


####################
# Main code        #
####################

# location of the coeqtl rds
coeqtl_rds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/coeqtl_replication/tf_subset/5DS_Meta_Analysis_Sign_TF_subset.RDS'
#
coeqtl <- readRDS(coeqtl_rds_loc)
# extract the variant-TF-gene
coeqtl_confinement <- unique(coeqtl[, c('snp_id', 'coeGene', 'eGene')])
# set column names
colnames(coeqtl_confinement) <- c('variant', 'TF', 'gene')
# get the eregulons
ereg_names_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_eregnames.txt.gz'
ereg_names <- fread(ereg_names_loc, header = T)[[1]]
# split eregs up by underscore
ereg_annotations <- str_split_fixed(ereg_names, '_', 4)
# set column names
colnames(ereg_annotations) <- c('TF', 'source', 'direction', 'set')
# add the original ereg as well
ereg_annotations <- cbind(data.frame('eregulon' = ereg_names), ereg_annotations)
# keep only the eregulons with the correct direction
ereg_annotations <- ereg_annotations[ereg_annotations[['direction']] %in% c('-/+', '+/+'), ]
# then keep direct if there is also extended
ereg_annotations <- ereg_annotations[order(ereg_annotations[['TF']], ereg_annotations[['source']], ereg_annotations[['direction']]), ]
ereg_annotations <- ereg_annotations[!duplicated(paste(ereg_annotations[['TF']], ereg_annotations[['direction']])), ]
# merge the eregulon onto the TF
coeqtl_eregs <- merge(coeqtl_confinement, ereg_annotations[, c('TF', 'eregulon')], by = 'TF')
# order columns and remove TF column
coeqtl_eregs <- coeqtl_eregs[, c('variant', 'eregulon', 'gene')]
# order to make it easier to search
coeqtl_eregs <- coeqtl_eregs[order(coeqtl_eregs[['variant']], coeqtl_eregs[['eregulon']], coeqtl_eregs[['gene']]), ]
# set where to store the result
coeqtl_confinement_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/coeqtl_replication/tf_subset/coeqtl_replication_set.tsv.gz'
# and write the result
write.table(coeqtl_eregs, gzfile(coeqtl_confinement_output_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(coeqtl_confinement_output_loc)
