#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_scenic_eregulon_vs_iegenes.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(ggplot2)


####################
# Functions        #
####################


rss_table_to_long <- function(rss_table) {
  # get the categories
  categories <- rss_table[[1]]
  # get the regulons
  regulons <- colnames(rss_table)[-1]
  # create the new table
  long_table <- data.frame('regulon' = regulons)
  # next go through the categories
  for (category in categories) {
    # add that to the table
    long_table[[category]] <- as.vector(unlist(rss_table[rss_table[[1]] == category, c(-1)]))
  }
  return(long_table)
}


#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @returns a list with the egenes per cell type
#' 
get_qtls_per_celltype_limix <- function(qtl_output_loc, output_file='/inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', significance_cutoffs=list('feature_q_value' = 0.05, 'total_bf_eigen' = 0.05), verbose=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  qtls_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(qtl_output_loc, '/', cell_type, '/', output_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    # filter the results on significance
    if (!is.null(significance_cutoffs)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
      }
      # check each cutoff
      for(significance_column in names(significance_cutoffs)) {
        # get the cutoff
        significance_cutoff <- significance_cutoffs[[significance_column]]
        # filter
        cell_type_output <- cell_type_output[
          !is.na(cell_type_output[[significance_column]]) &
            cell_type_output[[significance_column]] < significance_cutoff, 
        ]
      }
      if (verbose) {
        print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
      }
    }
    # add to the list
    qtls_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(qtls_per_celltype)
}


get_de_outputs <- function(mast_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', lfc_cutoff=NULL, to_ens=F, symbols.to.ensg.mapping='genes.tsv', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w'), de_method='mast_seurat', condition_name='timepoint'){
  # set up per cell type
  de_per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # set up per stim combination
    de_per_condition <- list()
    # check each stim
    for(stim in stims){
      for(stim2 in stims){
        try({
          if(stim != stim2){
            print(paste(cell_type, stim, stim2, sep = ' '))
            # paste the filepath together
            filepath <- NULL
            if ("limma" == de_method) {
              #filepath <- paste(mast_output_loc, cell_type, '_', condition_name,  stim, '_', stim2, '.tsv.gz', sep = '')
              filepath <- paste(mast_output_loc, cell_type, '_', condition_name, '.tsv.gz', sep = '')
            }
            else {
              filepath <- paste(mast_output_loc, cell_type, stim,stim2, '.tsv.gz', sep = '')
            }
            if (file.exists(filepath)) {
              # read the mast output
              mast <- read.table(filepath, header=T)
              # filter to only include the significant results
              mast <- mast[mast[[pval_column]] <= 0.05, ]
              # filter for only the positive lfc if required
              if(only_positive){
                mast <- mast[mast[[lfc_column]] < 0, ]
              }
              # filter for only the positive lfc if required
              if(only_negative){
                mast <- mast[mast[[lfc_column]] > 0, ]
              }
              # filter only ones with strong enough effect if required
              if (!is.null(lfc_cutoff)) {
                mast <- mast[abs(mast[[lfc_column]]) >= lfc_cutoff, ]
              }
              # confine in some way if reporting a max number of genes
              if(!is.null(max)){
                # by p if required
                if(max_by_pval){
                  mast <- mast[order(mast[[p_val_column]]), ]
                }
                # by lfc otherwise
                else{
                  mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
                }
                # subset to the number we requested if max was set
                mast <- mast[1:max,]
              }
              # grab the genes from the column names
              genes <- rownames(mast)
              if ("limma" == de_method) {
                genes <- mast[['feature']]
              }
              # convert the symbols to ensemble IDs
              if (to_ens) {
                mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
                mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
                genes <- mapping[match(genes, mapping$V2),"V1"]
              }
              # otherwise change the Seurat replacement back
              else{
                #genes <- gsub("-", "_", genes)
              }
              # change gene names
              if ("limma" == de_method) {
                mast[['feature']] <- genes
              }
              else {
                rownames(mast) <- genes
              }
              de_per_condition[[paste(stim, stim2, sep = '')]] <- mast
            }
          }
        })
      }
    }
    de_per_ct[[cell_type]] <- de_per_condition
  }
  return(de_per_ct)
}


get_color_coding_dict <- function() {
  # medhigh
  color_coding_dict <- list()
  color_coding_dict[["B"]] <- "#71BC4B"
  #color_coding_dict[['CD4_T_cells']] <- '#7FC97F'
  color_coding_dict[['CD4_T_cells']] <- '#153057'
  color_coding_dict[['CD4T']] <- '#153057'
  #color_coding_dict[['CD8_T_cells']] <- '#BEAED4'
  color_coding_dict[['CD8_T_cells']] <- '#009DDB'
  color_coding_dict[['CD8T']] <- '#009DDB'
  #color_coding_dict[['Dendritic_cells']] <- '#FDC086'
  color_coding_dict[['Dendritic_cells']] <- '#965EC8'
  color_coding_dict[['DC']] <- '#965EC8'
  color_coding_dict[['Endothelial_cells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia_cells']] <- '#F0027F'
  color_coding_dict[['Mast_cells']] <- '#BF5B17'
  color_coding_dict[['Mature_absorptive_enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature_secretory_enterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory_B']] <- '#D95F02'
  color_coding_dict[['Microfold_cell']] <- '#BEAED4'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Monocyte']] <- '#EDBA1B'
  color_coding_dict[['Naive_B_cells']] <- '#FDC086'
  color_coding_dict[['NK']] <- '#E64B50'
  #color_coding_dict[['Plasma_cells']] <- '#E7298A'
  color_coding_dict[['Plasma_cells']] <- '#DB8E00'
  color_coding_dict[['Stem_cells']] <- '#66A61E'
  color_coding_dict[['Stromal_cells']] <- '#8DD3C7'
  #color_coding_dict[['T_others']] <- '#A6761D'
  color_coding_dict[['T_others']] <- '#FF63B6'
  color_coding_dict[['Transit_amplifying_cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  #color_coding_dict[['CD4+ T cells']] <- '#7FC97F'
  color_coding_dict[['CD4+ T cells']] <- '#153057'
  color_coding_dict[['CD4+ T']] <- '#153057'
  #color_coding_dict[['CD8+ T cells']] <- '#BEAED4'
  color_coding_dict[['CD8+ T cells']] <- '#009DDB'
  color_coding_dict[['CD8+ T']] <- '#009DDB'
  #color_coding_dict[['Dendritic cells']] <- '#FDC086'
  color_coding_dict[['Dendritic cells']] <- '#965EC8'
  color_coding_dict[['Endothelial cells']] <- '#FFFFB3'
  color_coding_dict[['Endothelial\ncells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia cells']] <- '#F0027F'
  color_coding_dict[['MAST cells']] <- '#BF5B17'
  color_coding_dict[['Mature absorptive enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature\nabsorptive\nenterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature secretory enterocytes']] <- '#1B9E77'
  color_coding_dict[['Mature secretory\nenterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory B cells']] <- '#D95F02'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Microfold cells']] <- '#BEAED4'
  color_coding_dict[['Monocytes']] <- '#EDBA1B'
  color_coding_dict[['Naive B cells']] <- '#FDC086'
  #color_coding_dict[['Plasma cells']] <- '#E7298A'
  color_coding_dict[['Plasma cells']] <- '#DB8E00'
  color_coding_dict[['Stem cells']] <- '#66A61E'
  color_coding_dict[['Stromal cells']] <- '#8DD3C7'
  #color_coding_dict[['other T cells']] <- '#A6761D'
  color_coding_dict[['other T cells']] <- '#FF63B6'
  color_coding_dict[['Transit amplifying cells']] <- '#FF7F00'
  color_coding_dict[['Transit\namplifying cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  # up and down regulation will be added to, we need a whitening percentage
  pct_whitening <- 40
  # then we will check each cell type
  for (cell_type in names(color_coding_dict)) {
    # the up color is the same as the regular one
    color_coding_dict[[paste(cell_type, 'up')]] <- color_coding_dict[[cell_type]]
    # but the down one will have a more faded colour
    color_coding_dict[[paste(cell_type, 'down')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    # we'll do something similiar when we have multiple conditions
    color_coding_dict[[paste(cell_type, 'combined')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'UT')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, '24hCA')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "black"))(100)[pct_whitening]
  }
  # general
  color_coding_dict[['AI']] <- 'darkblue'
  color_coding_dict[['NI']] <- 'darkred'
  color_coding_dict[['Actively Inflamed']] <- 'darkblue'
  color_coding_dict[['Non-Inflamed']] <- 'darkred'
  color_coding_dict[['UT']] <- 'lightgrey'
  color_coding_dict[['24hCA']] <- 'forestgreen'
  return(color_coding_dict)
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
# get the enrichment statistics
rss_condition_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_condition.tsv.gz'
rss_cell_type_major_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_cell_type_lowerres.tsv.gz'
rss_cell_type_minor_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_cell_type_highres.tsv.gz'

# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# and enrichments
rss_condition <- fread(rss_condition_loc, header = T, sep = '\t')

# convert rss table to a long table
rss_condition_long <- rss_table_to_long(rss_condition)

# also make a plottble table
rss_condition_plottable <- rbind(
  data.frame('rss' = rss_condition_long[['UT']], 'condition' = rep('UT', times = nrow(rss_condition_long))), 
  data.frame('rss' = rss_condition_long[['24hCA']], 'condition' = rep('24hCA', times = nrow(rss_condition_long)))
)

# and plot that
ggplot(data = rss_condition_plottable, mapping = aes(x = rss, fill = condition)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = get_color_coding_dict()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('UT and 24hCA enrichment scores') +
  xlab('RSS eRegulon enrichment score') + 
  ylab('Density')

# calculate the difference between the two conditions
rss_condition_long[['delta']] <- abs(rss_condition_long[['24hCA']] - rss_condition_long[['UT']])

# plot the delta as wel
ggplot(data = rss_condition_long, mapping = aes(x = delta)) +
  geom_density(alpha = 0.5, fill = 'lightgray') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('UT and 24hCA enrichment score deltas') +
  xlab('RSS eRegulon enrichment score') + 
  ylab('Density')

# get the regulons more than .6
# condition_24hca_regulons <- rss_condition_long[rss_condition_long[['24hCA']] >= .5, ][['regulon']]
# sort by 24hCA score
rss_condition_long <- rss_condition_long[order(rss_condition_long[['24hCA']], decreasing = T), ]
# take the top 25 %
condition_24hca_regulons <- rss_condition_long[1 : round(.33 *nrow(rss_condition_long)), ][['regulon']]

# location of eQTL
ieqtl_output_combined_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# get all the ieqtl data
ieqtl_output <- get_qtls_per_celltype_limix(ieqtl_output_combined_loc)

# get the genes where the effect is stronger after stim
stim_iegenes <- list()
# and all the genes
all_iegenes <- list()
# check each cell type
for (ct in names(ieqtl_output)) {
  # grab the set
  ct_ieqtl_set <- ieqtl_output[[ct]]
  # extract the genes
  stim_iegenes[[ct]] <- ct_ieqtl_set[ct_ieqtl_set[['beta']] > 0, ][['feature']]
  all_iegenes[[ct]] <- ct_ieqtl_set[['feature']]
}
# add the cells into one list  
stim_iegenes <- unique(do.call('c', stim_iegenes))
all_iegenes <- unique(do.call('c', all_iegenes))

# get all genes
all_scenic_genes <- unique(scenic_output[['Gene']])
# and the stim genes
stim_scenic_genes <- unique(scenic_output[scenic_output[['Gene_signature_name']] %in% condition_24hca_regulons, ][['Gene']])

# get the non stim scenic genes
unstim_scenic_genes <- setdiff(all_scenic_genes, stim_scenic_genes)
# get the unstim iegenes
unstim_iegenes <- setdiff(all_iegenes, stim_iegenes)

# get the overlaps
stim_gene_in_scenic_stim_n <- length(intersect(stim_iegenes, stim_scenic_genes))
stim_gene_not_in_scenic_stim_n <- length(stim_iegenes) - stim_gene_in_scenic_stim_n
nonstim_gene_in_scenic_stim_n <- length(intersect(unstim_iegenes, stim_scenic_genes))
nonstim_gene_not_in_scenic_stim_n <- length(unstim_iegenes) - nonstim_gene_in_scenic_stim_n

# make contingency table
contingency_table <- matrix(
  c(stim_gene_in_scenic_stim_n, stim_gene_not_in_scenic_stim_n,
    nonstim_gene_in_scenic_stim_n, nonstim_gene_not_in_scenic_stim_n),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    ieqtl = c("stim", "unstim"),
    scenic = c("stimscenic", "notstimscenic")
  )
)
# show the table
contingency_table
#                   scenic
# ieqtl    stimscenic notstimscenic
# stim           66           100
# unstim         24            41
# do fisher-exact
fexact <- fisher.test(contingency_table)
# show fexact result
fexact
# Fisher's Exact Test for Count Data
# 
# data:  contingency_table
# p-value = 0.7648
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   0.6009864 2.1429326
# sample estimates:
#   odds ratio 
# 1.126917 
# 

# location of the DE output
all_de_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/'
# just get the DE genes
all_de_output <- get_de_outputs(all_de_output_loc, de_method = 'limma', pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, stims = c('24hCA', 'UT'), condition_name = 'condition_final')
# now collect all the DE genes, and the stim specific genes
all_de_genes <- list()
stim_de_genes <- list()
# check each cell type
for (ct in names(all_de_output)) {
  # extract the table
  de_output_ct <- all_de_output[[ct]][['UT24hCA']]
  # get all the genes
  all_de_genes[[ct]] <- de_output_ct[['feature']]
  # and the stim specific ones
  stim_de_genes[[ct]] <- de_output_ct[de_output_ct[['logFC']] > 0, ][['feature']]
}
# merge the gene lists
all_de_genes <- unique(do.call('c', all_de_genes))
stim_de_genes <- unique(do.call('c', stim_de_genes))
# get the setdiff
unstim_de_genes <- setdiff(all_de_genes, stim_de_genes)

# get the overlaps
de_stim_gene_in_scenic_stim_n <- length(intersect(stim_de_genes, stim_scenic_genes))
de_stim_gene_not_in_scenic_stim_n <- length(stim_de_genes) - de_stim_gene_in_scenic_stim_n
nonstim_de_gene_in_scenic_stim_n <- length(intersect(unstim_de_genes, stim_scenic_genes))
nonstim_de_gene_not_in_scenic_stim_n <- length(unstim_de_genes) - nonstim_de_gene_in_scenic_stim_n
# make contingency table
de_contingency_table <- matrix(
  c(de_stim_gene_in_scenic_stim_n, de_stim_gene_not_in_scenic_stim_n,
    nonstim_de_gene_in_scenic_stim_n, nonstim_de_gene_not_in_scenic_stim_n),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    de = c("stim", "unstim"),
    scenic = c("stimscenic", "notstimscenic")
  )
)
# show the table
de_contingency_table
#                 scenic
# de       stimscenic notstimscenic
# stim         4675          8988
# unstim       1076          6905
# do fisher-exact
de_fexact <- fisher.test(de_contingency_table)
# show fexact result
de_fexact
# 
# Fisher's Exact Test for Count Data
# 
# data:  de_contingency_table
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   3.100262 3.595132
# sample estimates:
#   odds ratio 
# 3.337702 

# sort by delta score
rss_condition_long <- rss_condition_long[order(rss_condition_long[['delta']], decreasing = T), ]
# take the top 25 %
condition_different_regulons <- rss_condition_long[1 : round(.33 *nrow(rss_condition_long)), ][['regulon']]
# get the genes that were specific to a condition
condition_scenic_genes <- unique(scenic_output[scenic_output[['Gene_signature_name']] %in% condition_different_regulons, ][['Gene']])
# get the condition scenic genes
non_condition_scenic_genes <- setdiff(all_scenic_genes, condition_scenic_genes)

# get the overlaps
scenic_condition_in_iegene_n <- length(intersect(all_iegenes, condition_scenic_genes))
scenic_condition_not_in_iegene_n <- length(condition_scenic_genes) - scenic_condition_in_iegene_n
scenic_noncondition_in_iegene_n <- length(intersect(non_condition_scenic_genes, all_iegenes))
scenic_noncondition_not_in_iegene_n <- length(non_condition_scenic_genes) - scenic_noncondition_in_iegene_n
# make contingency table
condition_ie_contingency_table <- matrix(
  c(scenic_condition_in_iegene_n, scenic_condition_not_in_iegene_n,
    scenic_noncondition_in_iegene_n, scenic_noncondition_not_in_iegene_n),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    scenic = c("condvar", "noncondvar"), 
    iegene = c("iegene", "notiegene")
  )
)
# show contingency table
condition_ie_contingency_table
#                   iegene
# scenic   stimscenic notstimscenic
# stim           76          4446
# unstim         27          2565
# do fisher-exact
condition_ie_fexact <- fisher.test(condition_ie_contingency_table)
# show fexact result
condition_ie_fexact
# 
# Fisher's Exact Test for Count Data
# 
# data:  condition_ie_contingency_table
# p-value = 0.03034
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.031263 2.628559
# sample estimates:
# odds ratio 
#   1.623828

# get the overlaps with DE again
scenic_condition_in_degene_n <- length(intersect(all_de_genes, condition_scenic_genes))
scenic_condition_not_in_degene_n <- length(condition_scenic_genes) - scenic_condition_in_degene_n
scenic_noncondition_in_degene_n <- length(intersect(non_condition_scenic_genes, all_de_genes))
scenic_noncondition_not_in_degene_n <- length(non_condition_scenic_genes) - scenic_noncondition_in_degene_n
# make contingency table
condition_de_contingency_table <- matrix(
  c(scenic_condition_in_degene_n, scenic_condition_not_in_degene_n,
    scenic_noncondition_in_degene_n, scenic_noncondition_not_in_degene_n),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    scenic = c("condvar", "noncondvar"), 
    degene = c("degene", "notdegene")
  )
)
# show contingency table
condition_de_contingency_table
#             degene
# scenic       degene notdegene
# condvar      4492        30
# noncondvar   2559        33
# do fisher-exact
condition_de_fexact <- fisher.test(condition_de_contingency_table)
# show fexact result
condition_de_fexact
# Fisher's Exact Test for Count Data
# 
# data:  condition_de_contingency_table
# p-value = 0.01183
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.138599 3.285574
# sample estimates:
# odds ratio 
#   1.930777 
