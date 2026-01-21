#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_frac_exp_files.R
# Function: create files that show the fraction of cells expressing a gene
#
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(mdfiver)


####################
# Main code        #
####################

# set the location of the ncell with >0 expression per donor
ncells_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/ncells_by_genes/L1/combined/'
# the prepend of the file
ncell_prepend <- ''
# the append
ncell_append <- '.tsv.gz'
# set the location of the smf file
smfs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/'
# the prepend of this file
smfs_prepend <- ''
# the append of this file
smfs_append <- 'smf.tsv.gz'
# set the location of the output file
fracs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/frac_exp/'
# the prepend of this file
fracs_prepend <- ''
# the append of this file
fracs_append <- '.tsv.gz'
# these are the cell types to use
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

# check each cell type
for (cell_type in cell_types) {
  # the path of the ncells file
  ncells_ct_loc <- paste0(ncells_loc, '/', ncell_prepend, cell_type, ncell_append)
  # read the file
  ncells_ct <- fread(ncells_ct_loc, header = T, sep = '\t')
  # path of the smf file
  smf_ct_file <- paste0(smfs_loc, '/', smfs_prepend, cell_type, '/', smfs_append)
  # read that file as well
  smf_ct <- fread(smf_ct_file, header = T, sep = '\t')
  # calculate the non-zero fraction
  ct_nonzero <- data.frame('feature' = ncells_ct[['feature']], 'frac_exp' = rowSums(ncells_ct[, -c('feature')]) / nrow(smf_ct))
  # set the output location
  ct_nonzero_loc <- paste0(fracs_loc, '/', fracs_prepend, cell_type, fracs_append)
  # write the table
  write.table(ct_nonzero, gzfile(ct_nonzero_loc), row.names = F, col.names = T, quote = F, sep = '\t')
  # and make a checksum
  mdfiver::create_md5_for_file(ct_nonzero_loc)
}
