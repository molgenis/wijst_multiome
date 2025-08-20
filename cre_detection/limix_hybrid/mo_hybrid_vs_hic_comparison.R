#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_hybrid_vs_hic_comparison.R
# Function: compare the limix hybrid CRE output against the hiC data
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(mdfiver)
library(ggplot2)
library(ggvenn)
library(doParallel)
library(svMisc)


####################
# Functions        #
####################

#' Calculate Minimum Distances Between Flanking Positions
#'
#' This function computes the pairwise distances between left and right flanking positions
#' from two sets of genomic coordinates provided in a data table. It returns a table of
#' distances and the minimum absolute distance for each row, with special handling for
#' overlapping flanks.
#'
#' @param position_table A `data.table` or `data.frame` containing genomic positions.
#' @param left_flank_column1 A string specifying the column name for the left flank of the first set.
#' @param right_flank_column1 A string specifying the column name for the right flank of the first set.
#' @param left_flank_column2 A string specifying the column name for the left flank of the second set.
#' @param right_flank_column2 A string specifying the column name for the right flank of the second set.
#'
#' @return A `data.table` with the following columns:
#' \describe{
#'   \item{lf1_to_lf2}{Distance from left flank 1 to left flank 2}
#'   \item{rf1_to_rf2}{Distance from right flank 1 to right flank 2}
#'   \item{lf1_to_rf2}{Distance from left flank 1 to right flank 2}
#'   \item{rf1_to_lf2}{Distance from right flank 1 to left flank 2}
#'   \item{min_dist}{Minimum absolute distance among the above, or 0 if flanks overlap}
#' }
#'
#' @examples
#' library(data.table)
#' dt <- data.table(
#'   lf1 = c(100, 200),
#'   rf1 = c(150, 250),
#'   lf2 = c(130, 220),
#'   rf2 = c(180, 270)
#' )
#' get_closest_flanks(dt, "lf1", "rf1", "lf2", "rf2")
#'
#' @export
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


#' Randomly Sample Regions Per Gene with Optional Distance Constraints
#'
#' This function performs gene-wise sampling of regions from a reference table, matching the number of regions in a "true" table. 
#' Optionally, it can constrain sampling based on distance comparisons, retrying several times if sampling fails.
#'
#' @param true_table A `data.frame` or `data.table` containing the original regions and genes.
#' @param table_to_sample_from A `data.frame` or `data.table` with regions to sample from.
#' @param region_column_true A string giving the column name for regions in `true_table`. Default is `"region"`.
#' @param gene_column_true A string giving the column name for genes in `true_table`. Default is `"gene"`.
#' @param region_column_sampling A string specifying the region column in `table_to_sample_from`. Default is `"region"`.
#' @param gene_column_sampling A string specifying the gene column in `table_to_sample_from`. Default is `"gene"`.
#' @param distance_column_true Optional string for the column holding distance values in `true_table`. If `NULL`, distance is ignored.
#' @param distance_column_sampling Optional string for the column holding distance values in `table_to_sample_from`.
#' @param distance_overshoot A numeric value defining the tolerance window around distances. Can be a proportion or fixed value. Default is `0.25`.
#' @param n_retries Integer number of times to retry sampling if suitable regions aren't found. Default is `10`.
#' @param filter_trues Logical. If `TRUE`, removes true regions from `table_to_sample_from` before sampling. Default is `FALSE`.
#'
#' @return A `data.frame` containing sampled regions and their corresponding genes.
#'
#' @examples
#' # Sample regions per gene ignoring distance
#' sampled <- randomly_sample_regions_per_gene(
#'   true_table = true_df,
#'   table_to_sample_from = background_df,
#'   region_column_true = "region",
#'   gene_column_true = "gene",
#'   filter_trues = TRUE
#' )
#'
#' @export
randomly_sample_regions_per_gene <- function(true_table, table_to_sample_from, region_column_true='region', gene_column_true='gene', region_column_sampling='region', gene_column_sampling='gene', distance_column_true=NULL, distance_column_sampling=NULL, distance_overshoot=.25, n_retries=10, filter_trues=F) {
  # get all unique genes
  true_genes <- unique(true_table[[gene_column_true]])
  # we'll store the samplings per gene first
  samplings_per_gene <- list()
  # check each gene
  samplings_per_gene <- foreach (i = 1: length(true_genes)) %dopar% {
    #for (i in 1:length(true_genes)) {
    # get the gene
    true_gene <- true_genes[i]
    # get the sample table for that gene
    sample_from_gene <- table_to_sample_from[table_to_sample_from[[gene_column_sampling]] == true_gene, ]
    # remove the true effects from the sampling if requested
    if (filter_trues) {
      sample_from_gene <- sample_from_gene[!(sample_from_gene[[region_column_sampling]] %in% true_table[true_table[[gene_column_true]] == true_gene, ][[region_column_true]]), ]
    }
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
        # set the number of retries
        n_tried <- 0
        # and whether we succeeded
        succeeded <- F
        while(n_tried <= n_retries & succeeded == F) {
          # and sort by distance
          sample_from_gene <- sample_from_gene[order(sample_from_gene[[distance_column_sampling]]), ]
          # get the actual data for this gene
          true_for_gene <- unique(true_table[true_table[[gene_column_true]] == true_gene, c(..gene_column_true, ..region_column_true, ..distance_column_true)])
          # sort by distance
          true_for_gene <- true_for_gene[order(true_for_gene[[distance_column_true]]), ]
          # get the regions for this gene
          true_regions <- true_for_gene[[region_column_true]]
          # make a vector of the samples genes
          sampled_for_gene_vector <- rep(NA, times = length(true_regions))
          # check each true region
          for (i_region in 1: length(true_regions)) {
            # get the region based on index
            true_region <- true_regions[i_region]
            # get the distance for this region
            region_distance <- true_table[true_table[[region_column_true]] == true_region, ][[distance_column_true]][1]
            # calculate the flanks
            flank_left <- NULL
            flank_right <- NULL
            # we'll do either as fraction of actual size
            if (distance_overshoot < 0) {
              stop('distance overshoot needs to be either a fraction or a positive number larger than 1')
            }
            else if (distance_overshoot > 1) {
              flank_left <- max(0, region_distance - distance_overshoot)
              flank_right <- region_distance + distance_overshoot
            }
            else if (distance_overshoot >= 0) {
              flank_left <- region_distance * (1 - distance_overshoot)
              flank_right <- region_distance * (1 + distance_overshoot)
            }
            # subset to regions that fall into this distance
            sampled_regions_comparable_distance <- sample_from_gene[
              #sample_from_gene[[region_column_sampling]] != true_region & 
              sample_from_gene[[distance_column_sampling]] >= flank_left & 
                sample_from_gene[[distance_column_sampling]] <= flank_right, ]
            # if there is anything to select, we can continue on
            if (nrow(sampled_regions_comparable_distance) > 0) {
              # then randomly select a region
              random_region_i <- sample(1 : nrow(sampled_regions_comparable_distance), 1)
              # and get that region
              region_randomly_selected <- sampled_regions_comparable_distance[random_region_i, ][[region_column_sampling]]
              # and put that into the vector
              sampled_for_gene_vector[i_region] <- region_randomly_selected
              # now remove that random region from out pool, so we don't select it again
              sample_from_gene <- sample_from_gene[sample_from_gene[[region_column_sampling]] != region_randomly_selected, ]
            }
            else {
              # otherwise we might have to try again
              n_tried <- n_tried + 1
              # reset 
              sampled_for_gene <- NULL
              # and break the loop
              break
            }
          }
          # if we got to the end, we succeeded
          succeeded <- T
        }
        # if we didn't succeed, let the user know
        if (!succeeded) {
          warning(paste0('could not sample for ', true_gene, ' after ', n_retries, ' retries, will skip this gene!'))
          sampled_for_gene <- data.frame('gene' = c(), 'region' = c())
        }
        else {
          sampled_for_gene <- data.frame('gene' = rep(true_gene, times = length(sampled_for_gene_vector)), 'region' = sampled_for_gene_vector)
        }
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


read_pseudobulk_cre_output_per_celltype <- function(pseudobulk_output_folder, cell_types=NULL, filename_output='qtl_results_all.txt.gz', significance_column='empirical_feature_p_value', significance_cutoff=0.05, add_mtc=T, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value', add_global_nominal_threshold=F, add_local_nominal_threshold=F, global_nominal_threshold_column_to_add='pval_nominal_threshold_global', local_nominal_threshold_column_to_add='pval_nominal_threshold_local', alpha_column='alpha_param', beta_column='beta_param', nominal_p_column='p_value', filter_alpha=T, alpha_min=.2, alpha_max=5, filter_significance=T) {
  # list all the files in the directory
  cell_type_folders <- list.dirs(pseudobulk_output_folder, full.names = F, recursive = F)
  # intersect the cell type folders with the cell types we are interested in
  if (!is.null(cell_types)) {
    cell_type_folders <- intersect(cell_type_folders, cell_types)
  }
  # save the results in a list
  output_per_celltype <- list()
  # now check each cell type
  for (cell_type in cell_type_folders) {
    # paste together the full file path
    cell_type_output_loc <- paste(pseudobulk_output_folder, cell_type, filename_output, sep = '/')
    # check if the file exists
    if (file.exists(cell_type_output_loc)) {
      # read this file
      cell_type_output <- fread(cell_type_output_loc, header = T, sep = '\t')
      # make sure there are no duplicates
      cell_type_output <- unique(cell_type_output)
      # filter on alpha if requested
      if (filter_alpha) {
        cell_type_output <- cell_type_output[!(cell_type_output[[alpha_column]] > alpha_max | cell_type_output[[alpha_column]] < alpha_min), ]
      }
      
      # get the features and the emperical p value
      if (add_mtc) {
        # subset to what we need
        cell_type_output_features <- NULL
        # which is a bit if we care about the nominal threshold
        if (add_global_nominal_threshold) {
          cell_type_output_features <- cell_type_output[, c(..feature_mtc_column, ..mtc_column, ..nominal_p_column, ..alpha_column, ..beta_column), with = F]
        }
        # even less if we don't try to get the nominal threshold as well
        else {
          cell_type_output_features <- cell_type_output[, c(..feature_mtc_column, ..mtc_column), with = F]
        }
        # remove the wherever we dont have our significance
        cell_type_output_features <- cell_type_output_features[!is.na(cell_type_output_features[[significance_column]]) & cell_type_output_features[[significance_column]] >= 0, ]
        # order by significance
        cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
        # keep only the first entry
        cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
        # set the values that are larger than 1, to be 1, problem with precision
        cell_type_output_features[cell_type_output_features[[mtc_column]] > 1, mtc_column] <- 1
        # add multiple testing correction
        cell_type_output_features[['qvalue']] <- qvalue(cell_type_output_features[[mtc_column]])$qvalues
        # now add back to the original table
        cell_type_output[[mtc_column_to_add]] <- cell_type_output_features[match(cell_type_output[[feature_mtc_column]], cell_type_output_features[[feature_mtc_column]]), 'qvalue'][['qvalue']]
        # based on this MTC column, we can now also add a cuttoff
        if (add_local_nominal_threshold) {
          cell_type_output_local_threshold <- calculate_nominal_thresholds(cell_type_output_features, fdr=significance_cutoff, pval_col=nominal_p_column, nominal_threshold_column='nomthres', cutoff_column = 'qvalue', alpha_column = alpha_column, beta_column = beta_column)
          # now add the nominal threshold to the full table
          cell_type_output[[local_nominal_threshold_column_to_add]] <- cell_type_output_local_threshold[match(cell_type_output[[feature_mtc_column]], cell_type_output_local_threshold[[feature_mtc_column]]), 'nomthres'][['nomthres']]
        }
        if(add_global_nominal_threshold) {
          # filter the output to significant MTC hits
          cell_type_output_features_significant <- cell_type_output_features[cell_type_output_features[['qvalue']] < significance_cutoff, ]
          # and get the maximum significant nominal value
          global_p_cutoff <- max(cell_type_output_features_significant[[nominal_p_column]])
          # add that to the table
          cell_type_output[[global_nominal_threshold_column_to_add]] <- global_p_cutoff
        }
      }
      # filter the file if requested
      if (filter_significance) {
        cell_type_output <- cell_type_output[
          cell_type_output[[significance_column]] < significance_cutoff, 
        ]
      }
      # add the cell type
      cell_type_output[['cell_type']] <- cell_type
      # put in the list
      output_per_celltype[[cell_type]] <- cell_type_output
    }
    else {
      warning(paste('folder exists at', cell_type_output_loc, 'but no file is there'))
    }
  }
  return(output_per_celltype)
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
# location of CRE output
cre_output_combined_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/'

# read the tables
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

# read hybrid method
hybrid_output_list <- read_pseudobulk_cre_output_per_celltype(cre_output_combined_loc, add_mtc = F, filter_alpha = F, add_global_nominal_threshold = F, add_local_nominal_threshold = F, filename_output = 'qtl_results_all_frac01.txt.gz', alpha_min = .2, alpha_max = 5, cell_types = c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
# merge them
hybrid_output <- do.call('rbind', hybrid_output_list)
# add z score
hybrid_output[['zscore']] <- hybrid_output[['beta']] / hybrid_output[['beta_se']]
# add a correlation based on the Z score, taking sample size and removing 2 + 10 PCs to get the degrees of freedom
hybrid_output[['r']] <- hybrid_output[['zscore']] / sqrt(hybrid_output[['zscore']]^2 + (hybrid_output[['n_samples']][1] - 12))
# add a p based z
hybrid_output[['z_from_p']] <- qnorm(hybrid_output[['p_value']] / 2) * -1 * sign(hybrid_output[['beta']])
# and add the distances
hybrid_output <- cbind(hybrid_output, cpeaks_anno[match(hybrid_output[['snp_id']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
hybrid_distances <- get_closest_flanks(hybrid_output, 'start_hg38', 'end_hg38', 'feature_start', 'feature_end')
hybrid_output[['distance']] <- hybrid_distances[['min_dist']]
hybrid_output[['category']] <- 'hybrid'
hybrid_output[['screen']] <- cpeaks_anno[match(hybrid_output[['snp_id']], cpeaks_anno[['signac_hg38']]), ][['screen_all']]
# split the output back into the list
for (cell_type in unique(hybrid_output[['cell_type']])) {
  # and put that in the list
  hybrid_output_list[[cell_type]] <- hybrid_output[!is.na(hybrid_output[['cell_type']]) & hybrid_output[['cell_type']] == cell_type, ]
}

# get significant region-gene pairs in SCENIC+
hybrid_output_significant <- hybrid_output[!is.na(hybrid_output[['feature_q_value']]) &
                                             hybrid_output[['feature_q_value']] < 0.05 &
                                             hybrid_output[['pval_nominal_threshold_global']] <= hybrid_output[['p_value']] &
                                             hybrid_output[['distance']] > 0 & 
                                             hybrid_output[['distance']] <= 150000, ]
# and non significant
hybrid_output_nonsignificant <- hybrid_output[!is.na(hybrid_output[['feature_q_value']]) &
                                             !(hybrid_output[['feature_q_value']] < 0.05) |
                                             !(hybrid_output[['pval_nominal_threshold_global']] <= hybrid_output[['p_value']]) &
                                             hybrid_output[['distance']] > 0 & 
                                             hybrid_output[['distance']] <= 150000, ]
# get the unique limix-sc region-gene pairs
hybrid_region_gene_pairs_sig <- unique(paste(hybrid_output_significant[['snp_id']], hybrid_output_significant[['feature_id']]))
hybrid_region_gene_pairs_nonsig <- unique(paste(hybrid_output_nonsignificant[['snp_id']], hybrid_output_nonsignificant[['feature_id']]))
# get the unique region-gene pairs in hiC
hic_region_gene_pairs <- hic[['r2g']]
# get significant limix-sc region-gene pairs in hic
hybrid_region_gene_pairs_sig_in_hic <- intersect(hybrid_region_gene_pairs_sig, hic_region_gene_pairs)
# get non-significant limix-sc region-gene pairs in hic
hybrid_region_gene_pairs_nonsig_in_hic <- intersect(hybrid_region_gene_pairs_nonsig, hic_region_gene_pairs)
# and the not in hic
hybrid_region_gene_pairs_sig_not_in_hic <- setdiff(hybrid_region_gene_pairs_sig, hic_region_gene_pairs)
hybrid_region_gene_pairs_nonsig_not_in_hic <- setdiff(hybrid_region_gene_pairs_nonsig, hic_region_gene_pairs)
# and the lengths
hybrid_region_gene_pairs_sig_in_hic_n <- length(hybrid_region_gene_pairs_sig_in_hic)
hybrid_region_gene_pairs_nonsig_in_hic_n <- length(hybrid_region_gene_pairs_nonsig_in_hic)
hybrid_region_gene_pairs_sig_not_in_hic_n <- length(hybrid_region_gene_pairs_sig_not_in_hic)
hybrid_region_gene_pairs_nonsig_not_in_hic_n <- length(hybrid_region_gene_pairs_nonsig_not_in_hic)

# make contingency table
contingency_table <- matrix(
  c(hybrid_region_gene_pairs_sig_in_hic_n, hybrid_region_gene_pairs_sig_not_in_hic_n,
    hybrid_region_gene_pairs_nonsig_in_hic_n, hybrid_region_gene_pairs_nonsig_not_in_hic_n),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    set = c("sig", "nonsig"),
    string = c("in_hic", "no_hic")
  )
)
# show the table
contingency_table
# set      in_hic no_hic
# sig       848   8248
# nonsig   7055  91712

# do fisher-exact
fexact <- fisher.test(contingency_table)
# show fexact result
fexact
# Fisher's Exact Test for Count Data
# 
# data:  contingency_table
# p-value = 1.698e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.238845 1.440579
# sample estimates:
# odds ratio 
#   1.336516 

# determine number of samplings
n_samplings <- 5
# now do a sampling per cell type
samplings_cell_types <- list()
# calculate the total samplings
n_samplings_tot <- n_samplings * length(names(hybrid_output_list))
# and initialize that number
sampling_i_tot <- 0
# check each cell type
for (cell_type in names(hybrid_output_list)) {
  # extract that table
  cres_hybrid_ct <- hybrid_output_list[[cell_type]]
  # remove anything that is not qvalue significant
  cres_hybrid_ct <- cres_hybrid_ct[cres_hybrid_ct[['feature_q_value']] < 0.05, ]
  # and where our cutoff defines it as not significant
  cres_hybrid_ct <- cres_hybrid_ct[cres_hybrid_ct[['p_value']] <= cres_hybrid_ct[['pval_nominal_threshold_global']], ]
  # and where the distance is within frame of what we expect
  cres_hybrid_ct <- cres_hybrid_ct[cres_hybrid_ct[['distance']] > 0 & cres_hybrid_ct[['distance']] <= 150000, ]
  # do x amount of samplings
  random_scenic_samplings_in_vs_out <- list()
  for (i in 1:n_samplings) {
    # increase sampling number
    sampling_i_tot <- sampling_i_tot + 1
    # show which iteration we are doing
    progress(sampling_i_tot, n_samplings_tot, progress.bar = TRUE)
    # do the sampling
    random_scenic_samplings_in_vs_out[[i]] <- randomly_sample_regions_per_gene(cres_hybrid_ct, window_pairs, region_column_true = 'snp_id', gene_column_true = 'feature_id', distance_column_true = 'distance', distance_column_sampling = 'distance', distance_overshoot = 10000, filter_trues = T)
  }
  # put that in the list
  samplings_cell_types[[cell_type]] <- random_scenic_samplings_in_vs_out
}
