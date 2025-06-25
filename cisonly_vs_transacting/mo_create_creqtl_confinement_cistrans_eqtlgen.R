#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_creqtl_confinement_cistrans_eqtlgen.R
# Function: create binomial CRE detection confinement file for cis-only vs trans-acting eQTLs
#
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(mdfiver)


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

# the overlap with the open chromatin
open_chromatin_overlap_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/eqtlgen_replications/cisonly_vs_transacting/qtl_tables/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping_credarnowindow.txt.gz'
# the overlap with the caQTLs
caqtl_overlap_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/eqtlgen_replications/cisonly_vs_transacting/qtl_tables/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping_caqtls.txt.gz'
# where we'll save the result
confinement_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/cisonly_vs_transacting/combined/confinements/cisonly_vs_transacting_full_confinement.tsv.gz'
confinement_r2g_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/cisonly_vs_transacting/combined/confinements/cisonly_vs_transacting_r2g.tsv.gz'

# these are the cell types we will check
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
# this is our openness cutoff
minimal_opennness <- 0.001
# and significance cutoff
significance_cutoff <- 0.05
# read the overlap
open_chromatin_overlap <- fread(open_chromatin_overlap_loc, header = T, sep = '\t')
# create a list with the overlaps
overlaps <- list()
# check each of the cell types
for (cell_type in cell_types) {
  # get the columns for this cell type
  ct_column <- paste0('openness_region_', cell_type)
  openness_column <- paste0('openness_openness_', cell_type)
  # get the rows that have sufficient openness
  entries_ct_open <- open_chromatin_overlap[!is.na(open_chromatin_overlap[[ct_column]]) & open_chromatin_overlap[[openness_column]] >= minimal_opennness, ]
  # get the columns we care about
  entries_ct_open <- entries_ct_open[, c('variant', ..ct_column, 'gene_name', 'type')]
  # rename the columns
  colnames(entries_ct_open) <- c('variant', 'region', 'gene', 'type')
  # add the cell type
  entries_ct_open[['cell_type']] <- cell_type
  # and finally what this was based on
  entries_ct_open[['inclusion']] <- 'openness'
  # add to the list
  overlaps[[paste0('openness', cell_type)]] <- entries_ct_open
}
# clear up memory
rm(open_chromatin_overlap)

# do the same for the caqtl overlap
caqtl_overlap <- fread(caqtl_overlap_loc, header = T, sep = '\t')
# check each of the cell types
for (cell_type in cell_types) {
  # get the columns for this cell type
  significance_column <- paste0('caqtl_p_value_', cell_type)
  # get the rows that have sufficient openness
  entries_ct_qtl <- caqtl_overlap[!is.na(caqtl_overlap[['caqtl_feature_id']]) & !is.na(caqtl_overlap[[significance_column]]) & caqtl_overlap[[significance_column]] < significance_cutoff, ]
  # get the columns we care about
  entries_ct_qtl <- entries_ct_qtl[, c('variant', 'caqtl_feature_id', 'gene_name', 'type')]
  # rename the columns
  colnames(entries_ct_qtl) <- c('variant', 'region', 'gene', 'type')
  # add the cell type
  entries_ct_qtl[['cell_type']] <- cell_type
  # and finally what this was based on
  entries_ct_qtl[['inclusion']] <- 'caqtl'
  # add to the list
  overlaps[[paste0('caqtl', cell_type)]] <- entries_ct_qtl
}

# merge all the tables
overlap_all <- do.call('rbind', overlaps)

# write the table
write.table(overlap_all, gzfile(confinement_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# make checksum
mdfiver::create_md5_for_file(confinement_loc)

# slim down the table
overlap_r2g <- unique(overlap_all[, c('region', 'gene')])
# and write that one
write.table(overlap_r2g, gzfile(confinement_r2g_loc), sep = '\t', row.names = F, col.names = F, quote = F)
# with a checksum again
mdfiver::create_md5_for_file(confinement_r2g_loc)
