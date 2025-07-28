#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_scenic_vs_string_comparison.R
# Function: compare the cre outputs of the different methods
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(mdfiver)
library(ggplot2)
library(ggvenn)


####################
# Functions        #
####################


####################
# Settings         #
####################

# whether we are in debug mode
debug <- F


####################
# Main code        #
####################

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# location of the string output
string_output_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/string_database/StringPairsEnsemblGenes.txt.gz'


# read the tables
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
string_output <- fread(string_output_loc, header = T, sep = '\t')

# location of ensemble ID to gene symbol mapping
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno <- fread(gene_anno_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')

# get the ensemble IDs for scenic
scenic_output[['TF_ens']] <- gene_anno[match(scenic_output[['TF']], gene_anno[['gs']]), ][['ens']]
# and for the gene
scenic_output[['Gene_ens']] <- gene_anno[match(scenic_output[['Gene']], gene_anno[['gs']]), ][['ens']]

# add an ordered combination of the two genes
string_output[['g2g']] <- apply(string_output, 1, function(x) {
  # get those genes
  genes <- c(x[['Gene1']], x[['Gene2']])
  # order them
  genes <- genes[order(genes)]
  # paste together
  genes_string <- paste(genes, collapse='_')
  return(genes_string)
})
# do the same for the scenic output
scenic_output[['g2g']] <- apply(scenic_output, 1, function(x) {
  # get those genes
  genes <- c(x[['TF_ens']], x[['Gene_ens']])
  # order them
  genes <- genes[order(genes)]
  # paste together
  genes_string <- paste(genes, collapse='_')
  return(genes_string)
})

# check how many TFs we have, and how many are also in the STRING database
tf_overlap_df <- data.frame(
  'in_string' = c('yes', 'no'), 
  'n' = c(
    length(unique(intersect(scenic_output[['TF_ens']], c(string_output[['Gene1']], string_output[['Gene2']])))),
    length(unique(scenic_output[['TF_ens']])) - length(unique(intersect(scenic_output[['TF_ens']], c(string_output[['Gene1']], string_output[['Gene2']]))))
  )
)
# same for the target genes
gene_overlap_df <- data.frame(
  'in_string' = c('yes', 'no'), 
  'n' = c(
    length(unique(intersect(scenic_output[['Gene_ens']], c(string_output[['Gene1']], string_output[['Gene2']])))),
    length(unique(scenic_output[['Gene_ens']])) - length(unique(intersect(scenic_output[['Gene_ens']], c(string_output[['Gene1']], string_output[['Gene2']]))))
  )
)
# make those into plots
ggplot(data = tf_overlap_df, mapping = aes(x = in_string, y = n, fill = in_string)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = list('yes' = 'darkblue', 'no' = 'darkred')) + 
  xlab('In STRING database') + 
  ylab('Number of TFs') + 
  ggtitle('Overlapping TFs of SCENIC+ in STRING') +
  theme(legend.position = 'none') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
ggplot(data = gene_overlap_df, mapping = aes(x = in_string, y = n, fill = in_string)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = list('yes' = 'darkblue', 'no' = 'darkred')) + 
  xlab('In STRING database') + 
  ylab('Number of target genes') + 
  ggtitle('Overlapping target genes of SCENIC+ in STRING') +
  theme(legend.position = 'none') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))


# check the overlap
ggvenn::ggvenn(
  data = list('SCENIC+' = unique(scenic_output[['g2g']]), 'STRING' = unique(string_output[['g2g']]))
) + ggtitle('Overlap of TF-gene links in SCENIC+ vs STRING')
