#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_correlate_onek1k_pseudobulks.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(progress)


####################
# Functions        #
####################


####################
# Settings        #
####################

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# get the paths to the files
output_loc_new <- '/groups/umcg-franke-scrna/tmp04/external_datasets/onek1k/qtl/interaction_eqtl/sc-eqtlgen/input/L1/combined/'
output_loc_old <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/input/L1/'
# the cell type to check
cell_types <- c('B', 'CD4_T', 'CD8_T', 'DC', 'Mono', 'NK')
# put correlation results in the list
cor_results <- list()
# set progress bar
pb <- progress_bar$new(total = length(cell_types))
# check each cell type
for (cell_type in cell_types) {
  # update progress bar
  pb$tick()
  # get the location of the old
  old_ct_path <- paste0(output_loc_old, '/', cell_type, '.Exp.txt.gz')
  # and the new
  new_ct_path <- paste0(output_loc_new, '/', cell_type, '.Exp.raw.txt.gz')
  # read in the old and new data
  old_ct <- read.table(old_ct_path, header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
  new_ct <- read.table(new_ct_path, header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
  # get the rownames is both
  genes_both <- intersect(rownames(old_ct), rownames(new_ct))
  # and make sure both are in the same order
  old_ct <- old_ct[genes_both, ]
  new_ct <- new_ct[genes_both, ]
  # calculate correlations
  cor_matrix_ct <- cor(old_ct, new_ct, use = "pairwise.complete.obs", method = "spearman")
  # put the correlation results in the list
  cor_results[[cell_type]] <- cor_matrix_ct
}
