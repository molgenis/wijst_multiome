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

random_sample_combinations <- function(true_table, sample_column1='TF_ens', sample_column2='Gene_ens', n_samplings=20) {
  # get the TFs and their occurences
  tf_occurences <- data.frame(table(unique(true_table[, c(..sample_column1, ..sample_column2)])[[sample_column1]]))
  # set column names more descriptive
  colnames(tf_occurences) <- c('sample1', 'occ')
  # get the unique genes
  sample2_values <- unique(true_table[[sample_column2]])
  
  # which we'll store in a list
  samplings <- list()
  # let's do each sampling
  for (sampling_i in 1:n_samplings) {
    # create a list to turn into a tf-gene table
    sampling_tbl_list <- list()
    # check each of the TFs
    for (tf_i in 1 : nrow(tf_occurences)) {
      # grab the tf
      tf <- tf_occurences[tf_i, 'sample1']
      # and the occurences
      occ <- tf_occurences[tf_i, 'occ']
      # now randomly get this many genes
      random_genes <- sample(sample2_values, size = occ)
      # make that into a df
      sampled_df_tf <- data.frame('sample1' = rep(tf, times = occ), 'sample2' = random_genes)
      colnames(sampled_df_tf) <- c(sample_column1, sample_column2)
      # put in the list for this sampling
      sampling_tbl_list[[tf]] <- sampled_df_tf
    }
    # merge the TFs of this sampling
    sampling_tbl <- do.call('rbind', sampling_tbl_list)
    # add the g2g again
    sampling_tbl[['g2g']] <- apply(sampling_tbl, 1, function(x) {
      # get those genes
      genes <- c(x[[sample_column1]], x[[sample_column2]])
      # order them
      genes <- genes[order(genes)]
      # paste together
      genes_string <- paste(genes, collapse='_')
      return(genes_string)
    })
    # put in the list
    samplings[[sampling_i]] <- sampling_tbl
  }
  return(samplings)
}


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


# get some stats for the true set
scenic_tf_gene_in_string <- length(unique(intersect(scenic_output[['g2g']], string_output[['g2g']])))
scenic_tf_gene_not_in_string <- length(unique(scenic_output[['g2g']])) - scenic_tf_gene_in_string

# get random samplings
samplings <- random_sample_combinations(true_table=scenic_output, sample_column1='TF_ens', sample_column2='Gene_ens', n_samplings=100)

# check each of the samplings
sampling_stats <- list()
for (sampling_i in 1 : length(samplings)) {
  # extract the random sampling
  sampling_tbl <- samplings[[sampling_i]]
  # check which of the random samplings are in string
  sampling_tf_gene_in_string <- length(unique(intersect(sampling_tbl[['g2g']], string_output[['g2g']])))
  sampling_tf_gene_not_in_string <- length(unique(sampling_tbl[['g2g']])) - sampling_tf_gene_in_string
  # make contingency table
  contingency_table <- matrix(
    c(scenic_tf_gene_in_string, scenic_tf_gene_not_in_string,
      sampling_tf_gene_in_string, sampling_tf_gene_not_in_string),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      set = c("scenic", "random"),
      string = c("in_string", "no_string")
    )
  )
  # do fisher-exact
  fexact <- fisher.test(contingency_table)
  # put in the list
  sampling_stats[[sampling_i]] <- fexact
}

# get the odds ratios of each comparison
ors <- rep(NA, times = length(sampling_stats))
# and the p-values
ps <- rep(NA, times = length(sampling_stats))
for (fexact_i in 1:length(sampling_stats)) {
  # grab the odds ratio
  ors[fexact_i] <- sampling_stats[[fexact_i]]$estimate
  ps[fexact_i] <- sampling_stats[[fexact_i]]$p.value
}
# plot the odds ratios
ggplot(data = data.frame(x = ors), mapping = aes(x = x)) + 
  geom_density(fill = 'darkred', alpha=.5) +
  xlab('Odds ratios') +
  ylab('Density') +
  ggtitle('Odds ratios of STRING enrichment in\nSCENIC vs random samplings from SCENIC') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
grid.text(paste('p =', formatC(max(ps), format = "e", digits = 2)), x = 0.8, y = 0.85, gp = gpar(fontsize = 16))

# plot the expected vs observed
max_p_index <- which(ps == max(ps))
# get that specific sampling
max_p_sampling <- samplings[[max_p_index]]
# and get the expected overlap
max_p_overlap <- length(unique(intersect(max_p_sampling[['g2g']], string_output[['g2g']])))
# make into a plot
p_overlap_expected_observed <- ggplot(data = data.frame('group' = c('Observed', 'Expected'), 'n' = c(scenic_tf_gene_in_string, max_p_overlap)), mapping = aes(x = group, y = n, fill= group)) + 
  geom_bar(stat = 'identity') + 
  ylim(c(0, scenic_tf_gene_in_string*1.25)) +
  scale_fill_manual(values = list('Observed' = 'darkgreen', 'Expected' = 'gray')) +
  ylab('Number of TF-target gene combinations\noverlapping with STRING') +
  xlab('') +
  ggtitle('Odds ratios of STRING enrichment in\nSCENIC+ vs random samplings from SCENIC+') +
  theme(legend.position = 'none') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  # horizontal comparison
  geom_segment(data = NULL, aes(x = 'Expected', xend = 'Observed', y = scenic_tf_gene_in_string*1.1, yend = scenic_tf_gene_in_string*1.1), colour = "darkgray", size = 1) +
  # vertical comparison on Expected side
  geom_segment(data = NULL, aes(x = 'Expected', xend = 'Expected', y = scenic_tf_gene_in_string*1.1, yend = max_p_overlap*1.05), colour = "darkgray", size = 1) +
  # vertical comparison on Observed side
  geom_segment(data = NULL, aes(x = 'Observed', xend = 'Observed', y = scenic_tf_gene_in_string*1.1, yend = scenic_tf_gene_in_string*1.05), colour = "darkgray", size = 1) +
  # add OR text
  annotate('text', x = 1.5, y = scenic_tf_gene_in_string*1.2, label = paste('OR =', 1.6), size = 5) +
  # add p text
  annotate('text', x = 1.5, y = scenic_tf_gene_in_string*1.15, label = paste('p =', formatC(max(ps), format = "e", digits = 2)), size = 5) +
  # make ticks bigger
  theme(axis.text.x = element_text(size = 14)) +
  # make y label bigger
  theme(axis.title.y = element_text(size = 14))

# show the plot
p_overlap_expected_observed
