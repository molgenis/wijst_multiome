#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_scenic_vs_hic_comparison.R
# Function: compare the cre outputs against the hiC data
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(mdfiver)
library(ggplot2)
library(ggvenn)
library(doParallel)


####################
# Functions        #
####################

get_closest_flanks <- function(position_table, left_flank_column1, right_flank_column1, left_flank_column2, right_flank_column2) {
  # get the distance between left flanks
  dist_left_flank1_to_left_flank2 <- position_table[[left_flank_column1]] - position_table[[left_flank_column2]]
  # distance between the right flanks
  dist_right_flank1_to_right_flank2 <- position_table[[right_flank_column1]] - position_table[[right_flank_column2]]
  # distance between left flank 1 and right flank 2
  dist_left_flank1_to_right_flank2 <- position_table[[left_flank_column1]] - position_table[[right_flank_column2]]
  # distance between right flank 1 and left flank 2
  dist_right_flank1_to_left_flank2 <- position_table[[right_flank_column1]] - position_table[[left_flank_column2]]
  # put in a table for convenience sake
  distances_tbl <- data.table(
    'lf1_to_lf2' = dist_left_flank1_to_left_flank2, 
    'rf1_to_rf2' = dist_right_flank1_to_right_flank2, 
    'lf1_to_rf2' = dist_left_flank1_to_right_flank2, 
    'rf1_to_lf2' = dist_right_flank1_to_left_flank2
  )
  # add the minimum absolute distance
  distances_tbl[['min_dist']] <- apply(distances_tbl, 1, function(x) {
    return(min(abs(x)))
  })
  # but set this to zero if any of the flanks end in the bodies
  #              -----
  #                 ++++
  distances_tbl[(distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['lf1_to_rf2']] > 0) |
                  #                   ----
                #                 ++++
                (distances_tbl[['rf1_to_lf2']] > 0 & distances_tbl[['lf1_to_rf2']] < 0) |
                  #                   ----
                #                 +++++++++
                (distances_tbl[['lf1_to_lf2']] > 0 & distances_tbl[['rf1_to_rf2']] < 0) |
                  #                 ---------
                #                   ++++
                (distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['rf1_to_rf2']] > 0)
                , 'min_dist'] <- 0
  return(distances_tbl)
}


randomly_sample_regions_per_gene <- function(true_table, table_to_sample_from, region_column_true='region', gene_column_true='gene', region_column_sampling='region', gene_column_sampling='gene', distance_column_true=NULL, distance_column_sampling=NULL, distance_overshoot=.25) {
  # get all unique genes
  true_genes <- unique(true_table[[gene_column_true]])
  # we'll store the samplings per gene first
  samplings_per_gene <- list()
  # check each gene
  #for (true_gene in true_genes) {
  samplings_per_gene <- foreach (i = 1: length(true_genes)) %dopar% {
    # get the gene
    true_gene <- true_genes[i]
    # get the sample table for that gene
    sample_from_gene <- table_to_sample_from[table_to_sample_from[[gene_column_sampling]] == true_gene, ]
    # check if we can sample
    if (nrow(sample_from_gene) > 0) {
      # initialize the sampling result
      sampled_for_gene <- NULL
      # if we don't care about distance, this is pretty simple
      if (is.null(distance_column_true) & is.null(distance_column_sampling)) {
        # we just randomly sample the number of regions for that gene from the original table
        n_random_to_sample <- nrow(unique(true_table[true_table[[gene_column_true]] == true_gene, c(..gene_column_true, ..region_column_true)]))
        random_region_indices <- sample(x = 1 : nrow(sample_from_gene), size = n_random_to_sample)
        # grab those random regions using the indices from the table to sample from, and put that in a table
        sampled_for_gene <- data.frame('gene' = rep(true_gene, times = n_random_to_sample), 'region' = sample_from_gene[[region_column_sampling]][random_region_indices])
      } else if ((is.null(distance_column_true) & !is.null(distance_column_sampling)) | (!is.null(distance_column_true) & is.null(distance_column_sampling))) {
        stop('distance_column_true and distance_column_sampling both need to be either true or false')
      } else {
        stop('TODO')
      }
      # put the sampling in the list
      #samplings_per_gene[[true_gene]] <- sampled_for_gene
      return(sampled_for_gene)
    }
    else {
      warning(paste0('could not sample for ', true_gene, ' due to it not being in the sample table, this one will be skipped!'))
    }
  }
  # merge all
  sampled_all <- do.call('rbind', samplings_per_gene)
  return(sampled_all)
}

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
# number of parallel cores
registerDoParallel(cores=4)


####################
# Main code        #
####################

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# location of the hiC output
hic_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cpeaks_to_screenv4_blood_hic.tsv.gz'
# location of 150k window file
window_pairs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/featureVariantFile.w150k.filtered0.0001_cts.txt'

# read the tables
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
hic <- fread(hic_loc, header = T, sep = '\t')
window_pairs <- fread(window_pairs_loc, header = T, sep = '\t')

# rename the columns for the window pairs
colnames(window_pairs) <- c('gene', 'region')

# location of ensemble ID to gene symbol mapping
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/input/LimixAnnotationFile.txt'
gene_anno <- fread(gene_anno_loc, header = T, sep = '\t')

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])

# get the ensemble IDs for scenic
scenic_output[['TF_ens']] <- gene_anno[match(scenic_output[['TF']], gene_anno[['gs']]), ][['ens']]
# and for the gene
scenic_output[['Gene_ens']] <- gene_anno[match(scenic_output[['Gene']], gene_anno[['gs']]), ][['ens']]

# add the chromosome and positions of the region
window_pairs[, c('region_chr', 'region_start', 'region_end')] <- cpeaks_anno[match(window_pairs[['region']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')]
# and the gene information as well
window_pairs[, c('gene_chr', 'gene_start', 'gene_end')] <- gene_anno[match(window_pairs[['gene']], gene_anno[['feature_id']]), c('chromosome', 'start', 'end')]
# making the gene chr the same format as the region format
window_pairs[['gene_chr']] <- paste0('chr', window_pairs[['gene_chr']])
# get the distances
window_pairs_distances <- get_closest_flanks(window_pairs, 'region_start', 'region_end', 'gene_start', 'gene_end')
# add those distances
window_pairs[['distance']] <- window_pairs_distances[['min_dist']]
# make sure to make the distance NA if the chromosomes were different (should never be the case, but just in case that might be the input)
window_pairs[window_pairs[['region_chr']] != window_pairs[['gene_chr']], 'distance'] <- NA

# we'll do the same thing for the hi-C data, to see if the distances are comparable
hic[, c('region_chr', 'region_start', 'region_end')] <- cpeaks_anno[match(hic[['region']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')]
# and the gene information as well
hic[, c('gene_chr', 'gene_start', 'gene_end')] <- gene_anno[match(hic[['gene']], gene_anno[['feature_id']]), c('chromosome', 'start', 'end')]
# making the gene chr the same format as the region format
hic[['gene_chr']] <- paste0('chr', hic[['gene_chr']])
# get the distances
hic_distances <- get_closest_flanks(hic, 'region_start', 'region_end', 'gene_start', 'gene_end')
# add those distances
hic[['distance']] <- hic_distances[['min_dist']]
# make sure to make the distance NA if the chromosomes were different (should never be the case, but just in case that might be the input)
hic[hic[['region_chr']] != hic[['gene_chr']], 'distance'] <- NA
# add region to gene column
hic[['r2g']] <- paste(hic[['region']], hic[['gene']])

# and finally for scenic as well
scenic_output[, c('region_chr', 'region_start', 'region_end', 'signac_region_name')] <- cpeaks_anno[match(scenic_output[['Region']], cpeaks_anno[['scenic_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38', 'signac_hg38')]
# and the gene information as well
scenic_output[, c('gene_chr', 'gene_start', 'gene_end')] <- gene_anno[match(scenic_output[['Gene']], gene_anno[['feature_id']]), c('chromosome', 'start', 'end')]
# making the gene chr the same format as the region format
scenic_output[['gene_chr']] <- paste0('chr', scenic_output[['gene_chr']])
# get the distances
scenic_output_distances <- get_closest_flanks(scenic_output, 'region_start', 'region_end', 'gene_start', 'gene_end')
# add those distances
scenic_output[['distance']] <- scenic_output_distances[['min_dist']]
# make sure to make the distance NA if the chromosomes were different (should never be the case, but just in case that might be the input)
scenic_output[scenic_output[['region_chr']] != scenic_output[['gene_chr']], 'distance'] <- NA
# add a region to gene column
scenic_output[['r2g']] <- paste(scenic_output[['signac_region_name']], scenic_output[['Gene']])

# subset to autosomal for scenic
scenic_output_autosomal <- scenic_output[scenic_output[['region_chr']] %in% paste0('chr', 1:22), ]

# randomly sample regions to genes without taking the region into consideration
random_scenic_samplings_noregion <- list()
for (i in 1:20) {
  random_scenic_samplings_noregion[[i]] <- randomly_sample_regions_per_gene(scenic_output_autosomal, window_pairs, region_column_true = 'signac_region_name', gene_column_true = 'Gene')
}

# get true overlap
scenic_r_gene_in_hic <- length(unique(intersect(scenic_output_autosomal[['r2g']], hic[['r2g']])))
scenic_r_gene_not_in_hic <- length(unique(scenic_output_autosomal[['r2g']])) - scenic_r_gene_in_string

# check each of the samplings
sampling_stats <- list()
for (sampling_i in 1 : length(random_scenic_samplings_noregion)) {
  # extract the random sampling
  sampling_tbl <- random_scenic_samplings_noregion[[sampling_i]]
  # add region to gene
  sampling_tbl[['r2g']] <- paste(sampling_tbl[['region']], sampling_tbl[['gene']])
  # check which of the random samplings are in string
  sampling_r_gene_in_hic <- length(unique(intersect(sampling_tbl[['r2g']], hic[['r2g']])))
  sampling_r_gene_not_in_hic <- length(unique(sampling_tbl[['r2g']])) - sampling_r_gene_in_hic
  # make contingency table
  contingency_table <- matrix(
    c(scenic_r_gene_in_hic, scenic_r_gene_not_in_hic,
      sampling_r_gene_in_hic, sampling_r_gene_not_in_hic),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      set = c("scenic", "random"),
      string = c("in_hic", "no_hic")
    )
  )
  # do fisher-exact
  fexact <- fisher.test(contingency_table)
  # put in the list
  sampling_stats[[sampling_i]] <- fexact
}

