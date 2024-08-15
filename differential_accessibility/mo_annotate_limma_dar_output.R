#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_annotate_limma_dar_output.R
# Function: annotate the DAR output to have the closest gene added
############################################################################################################################

####################
# libraries        #
####################

library(Signac)

####################
# Functions        #
####################


####################
# Main Code        #
####################

# location of the Signac objects
signac_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/'
signac_prepend <- 'mo_cpeaks_filtered_'
signac_append <- '_wstatus_1_80_20240709.rds'

# location of the dar output
dar_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/stimulation/pct01/'
dar_output_prepend <- ''
dar_output_append <- '_condition_final.wpermfdr.tsv.gz'

# the celltypes we'll use
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

# do each cell type
for (cell_type in cell_types) {
  # load the Signac object
  signac_object_loc_full <- paste(signac_object_loc, signac_prepend, tolower(cell_type), signac_append, sep = '')
  signac_object <- readRDS(signac_object_loc_full)
  # read the dar output
  dar_output_loc_full <- paste(dar_output_loc, dar_output_prepend, cell_type, dar_output_append, sep = '')
  dar_output <- read.table(dar_output_loc_full, header = T, sep = '\t')
  # get the genomic ranges from the cell type
  genomic_ranges_celltype <- granges(signac_object)
  # get the closes feature
  genomic_range_to_gene <- ClosestFeature(signac_object, genomic_ranges_celltype)
  # turn into dataframe
  genomic_range_to_gene <- data.frame(genomic_range_to_gene)
  # now merge the granges info onto the dataframe
  dar_output <- cbind(dar_output, data.frame(genomic_range_to_gene[match(dar_output[['feature']], genomic_range_to_gene[['query_region']]), c('tx_id', 'gene_name', 'gene_id', 'gene_biotype', 'type', 'distance')]))
  # write the result
  dar_output_result_loc <- paste(dar_output_loc, dar_output_prepend, cell_type, '_condition_final.wpermfdr.wgene.tsv.gz', sep = '')
  write.table(dar_output, gzfile(dar_output_result_loc), sep = '\t', row.names = F, col.names = T, quote = F)
}
