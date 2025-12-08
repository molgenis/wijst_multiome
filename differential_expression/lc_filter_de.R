#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lc_filter_de.R
# Function: perform differential gene expression analysis using Limma Dream
############################################################################################################################


####################
# libraries        #
####################

library(mdfiver)


####################
# Functions        #
####################


####################
# Main Code        #
####################

# location of the output to replicate
replication_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/replication/covid_replication.tsv'
# read the replication stats
replication <- read.table(replication_loc, header = T, sep = '\t')
# get the list of cell types we need to filter
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK', 'plasmablast')
# the location of the original output
de_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/LONG_COVID/'
# append of DE files
de_output_append <- '_LONG_COVID_final.tsv.gz'
# prepend of DE files
de_output_prepend <- ''
# where to put the subset of data
de_output_subset_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/LONG_COVID_replication/'
# append of DE files
de_output_subset_append <- '_LONG_COVID_final.tsv.gz'
# prepend of DE files
de_output_subset_prepend <- ''
# read the location of the genes
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno <- read.table(gene_anno_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')
# check each of the cell types
for (cell_type in cell_types) {
    # paste full path together
    de_output_loc_full <- paste0(de_output_loc, '/', de_output_prepend, cell_type, de_output_append)
    # read the file
    de_output_ct <- read.table(de_output_loc_full, header = T, sep = '\t')
    # add the ensemble ID
    de_output_ct[['ensemble']] <- gene_anno[match(de_output_ct[['feature']], gene_anno[['gs']]), ][['ens']]
    # filter by what was significant in replication
    de_output_ct <- de_output_ct[de_output_ct[['ensemble']] %in% replication[['ensemble']], ]
    # redo the B&H
    de_output_ct[['adj.P.Val']] <- p.adjust(de_output_ct[['P.Value']], method = 'BH')
    # redo bonferroni
    de_output_ct[['p.bonferroni']] <- p.adjust(de_output_ct[['P.Value']], method = 'bonferroni')
    # paste full output together
    de_output_subset_loc_full <- paste0(de_output_subset_loc, '/', de_output_subset_prepend, cell_type, de_output_subset_append)
    # write the result
    write.table(de_output_ct, gzfile(de_output_subset_loc_full), row.names = F, col.names = T, sep = '\t')
    # make checksum
    mdfiver::create_sha256_for_file(de_output_subset_loc_full)
}
