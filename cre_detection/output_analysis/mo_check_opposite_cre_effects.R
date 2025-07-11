#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_opposite_cre_effects.R
# Function: plot the CREs that seem to have an opposite direction in between the different CRE mapping methods
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(ggplot2)
library(roycols)
library(data.table)
library(qvalue)
library(cowplot)
library(snpStats)


####################
# Functions        #
####################

add_binarized_assay <- function(seurat_object, assay_to_binarize='peaks', layer_to_binarize='counts', assay_to_add='binpeaks') {
  # set the default assay
  DefaultAssay(seurat_object) <- assay_to_binarize
  # extract the assay
  assay_data <- Seurat::GetAssayData(seurat_object, layer = layer_to_binarize, assay = assay_to_binarize)
  # binarize the chromatin data if requested
  assay_data@x[assay_data@x > 0] <- 1
  # add the assay data back
  if (layer_to_binarize == 'counts') {
    seurat_object[[assay_to_add]] <- CreateAssayObject(counts = assay_data, cells = colnames(assay_data))
  }
  else {
    seurat_object[[assay_to_add]] <- CreateAssayObject(data = assay_data, cells = colnames(assay_data))
  }
  return(seurat_object)
}


normalize_mj <- function(seurat_object, assay_to_add='MJ', assay_to_normalize='RNA', layer_to_normalize='counts', slot_to_normalize='counts') {
  # set the assay
  DefaultAssay(seurat_object) <- assay_to_normalize
  # we'll extract the counts
  count_matrix <- NULL
  # with the native version depending on the Seurat object config
  if ('layers' %in% slotNames(seurat_object[[assay_to_normalize]])) {
    count_matrix <- Seurat::GetAssayData(seurat_object, layer = layer_to_normalize, assay = assay_to_normalize)
  }
  else {
    count_matrix <- Seurat::GetAssayData(seurat_object, slot = slot_to_normalize, assay = assay_to_normalize)
  }
  # ignore genes that are never expressed
  count_matrix <-  count_matrix[which(rowSums(count_matrix) != 0), ]
  # create new object to store the counts in
  norm_count_matrix <- count_matrix
  # do mean sample-sum normalization
  sample_sum_info = colSums(norm_count_matrix)
  mean_sample_sum = mean(sample_sum_info)
  sample_scale = sample_sum_info / mean_sample_sum
  # divide each column by sample_scale
  norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
  if ('layers' %in% slotNames(seurat_object[[assay_to_normalize]])) {
    print('using Seurat v5 style \'layer\'')
    seurat_object[[assay_to_add]] <- CreateAssay5Object(data = norm_count_matrix)
    
  } else {
    print('using Seurat v3/4 style \'slot\'')
    seurat_object[[assay_to_add]] <- CreateAssayObject(data = norm_count_matrix)
  }
  return(seurat_object)
}


plot_binarized_vs_other <- function(seurat_object, binarized_feature_name, other_feature_name, binarized_assay='binpeaks', binarized_layer='counts', other_assay='MJ', other_layer='data', pointless=F, legendless=T, ylim=NULL, paper_style=T, angle_labels=F, better_colors=T) {
  # extract the binarized feature
  binarized_features <- Seurat::GetAssayData(seurat_object, layer = binarized_layer, assay = binarized_assay)
  # get the specific featurs
  binarized_feature <- as.vector(unlist(binarized_features[binarized_feature_name, ]))
  # make into a character, as it is binarized
  binarized_feature <- as.character(binarized_feature)
  # extract the other feature
  other_features <- Seurat::GetAssayData(seurat_object, layer = other_layer, assay = other_assay)
  # get the specific feature
  other_feature <- as.vector(unlist(other_features[other_feature_name, ]))
  # merge these
  features_both <- data.frame(x = binarized_feature, y = other_feature)
  # make into a boxplot
  p <- ggplot(data = features_both, mapping = aes(x = x, y = y, fill = x)) +
    geom_boxplot() +
    xlab(binarized_feature_name) +
    ylab(other_feature_name) +
    geom_boxplot(outlier.shape = NA) + 
    # and add jitter
    geom_jitter(size = 0.5, alpha = 0.5)
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  if (better_colors) {
    p <- p + scale_fill_manual(values = roycols::get_color_list(features_both[['x']]))
  }
  return(p)
}

#' Calculate Nominal Thresholds
#'
#' This function calculates nominal thresholds for p-values based on a given false discovery rate (FDR).
#'
#' @param res_df A data frame containing the results with p-values and other relevant columns.
#' @param fdr A numeric value specifying the false discovery rate threshold. Default is 0.05.
#' @param pval_col A character string specifying the name of the column with p-values. Default is 'p_value'.
#' @param nominal_threshold_column A character string specifying the name of the column to store the nominal thresholds. Default is 'pval_nominal_threshold'.
#' @param cutoff_column A character string specifying the name of the column with feature q-values. Default is 'feature_q_value'.
#' @param alpha_column A character string specifying the name of the column with alpha parameters for the beta distribution. Default is 'alpha_param'.
#' @param beta_column A character string specifying the name of the column with beta parameters for the beta distribution. Default is 'beta_param'.
#'
#' @return A data frame with an additional column for nominal thresholds.
#' @export
#'
#' @examples
#' \dontrun{
#'   res_df <- data.frame(
#'     p_value = runif(100),
#'     feature_q_value = runif(100),
#'     alpha_param = rep(1, 100),
#'     beta_param = rep(1, 100)
#'   )
#'   calculate_nominal_thresholds(res_df)
#' }
calculate_nominal_thresholds <- function(res_df, fdr=0.05, pval_col='p_value', nominal_threshold_column='pval_nominal_threshold', cutoff_column='feature_q_value', alpha_column='alpha_param', beta_column='beta_param') {
  # get the lowerbound p values, so the ones that are smaller than the FDR
  indices_lb <- res_df[[cutoff_column]] < fdr
  lb <- as.vector(res_df[indices_lb, ][[pval_col]])
  # put then in ascending order
  lb <- lb[order(lb)]
  # get the upperbound p values, so the ones that are bigger than the FDR
  indices_ub <- res_df[[cutoff_column]] > fdr
  ub <- as.vector(res_df[indices_ub, ][[pval_col]])
  # and order them
  ub <- ub[order(ub)]
  
  # if we have any significant effects, we can get a cutoff
  if (length(lb) > 0) {
    # get the highest (p) significant value
    highest_in_lb <- tail(lb, 1)
    # if there are any non significant effects
    if (length(ub) > 0) {
      # get the lowest (p) non-significant value
      lowest_in_ub <- head(ub, 1)
      # and calculate the threshold
      pthreshold <- (highest_in_lb + lowest_in_ub) / 2
    } else {
      # otherwise the highest effect will just be the cutoff
      pthreshold <- highest_in_lb
    }
    # ge the threshold, based on the shapes of the beta distribution and the significance threshold
    res_df[[nominal_threshold_column]] <- stats::qbeta(pthreshold, as.vector(res_df[[alpha_column]]), as.vector(res_df[[beta_column]]))
  }
  else {
    # otherwise it would have to be zero
    res_df[[nominal_threshold_column]] <- 0
  }
  return(res_df)
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
  return(output_per_celltype)
}


read_binomial_output_per_celltype <- function(binomial_output_folder, cell_types=NULL, filename_output='meta_result.tsv.gz', significance_column='meta_q', significance_cutoff=0.05) {
  # just use the pseudobulk function
  output_per_celltype <- read_pseudobulk_cre_output_per_celltype(binomial_output_folder, cell_types = cell_types, filename_output = filename_output, significance_column = significance_column, significance_cutoff = significance_cutoff, add_mtc = F, add_global_nominal_threshold = F, add_local_nominal_threshold = F, filter_alpha = F)
  return(output_per_celltype)
}


plot_bino_z_scores <- function(zscore_matrix, region, gene, region_column='region', gene_column='gene', pointless=F, legendless=T, ylim=NULL, paper_style=T, angle_labels=F, better_colors=T) {
  # get the columns that are not the region and gene
  zscore_matrix_numbers <- zscore_matrix[, -c(..region_column, ..gene_column)]
  # subset to the region and gene
  zscores_region_gene <- as.vector(unlist(zscore_matrix_numbers[zscore_matrix[[region_column]] == region & zscore_matrix[[gene_column]] == gene, ]))
  # make into a table
  zscores_plot_frame <- data.frame(x = rep('samples', times = length(zscores_region_gene)), y = zscores_region_gene)
  # make into a boxplot
  p <- ggplot(data = zscores_plot_frame, mapping = aes(x = x, y = y, fill = x)) +
    geom_boxplot() +
    xlab('samples') +
    ylab('Z-score') +
    geom_boxplot(outlier.shape = NA) + 
    # and add jitter
    geom_jitter(size = 0.5, alpha = 0.5)
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  if (better_colors) {
    p <- p + scale_fill_manual(values = roycols::get_color_list(zscores_plot_frame[['x']]))
  }
  return(p)
}

plot_bino_z_distribution <- function(zscore_matrix, region, gene, region_column='region', gene_column='gene', pointless=F, legendless=T, ylim=NULL, paper_style=T, angle_labels=F) {
  # get the columns that are not the region and gene
  zscore_matrix_numbers <- zscore_matrix[, -c(..region_column, ..gene_column)]
  # subset to the region and gene
  zscores_region_gene <- as.vector(unlist(zscore_matrix_numbers[zscore_matrix[[region_column]] == region & zscore_matrix[[gene_column]] == gene, ]))
  # make into a table
  zscores_plot_frame <- data.frame(x = zscores_region_gene)
  # make into a boxplot
  p <- ggplot(data = zscores_plot_frame, mapping = aes(x = x)) +
    geom_density() +
    xlab('Z-score') +
    ylab('Density') +
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p)
}


plot_pseudobulk_cre <- function(accessibility_table, expression_table, region, gene, region_column='feature_id', gene_column='feature_id', pointless=F, legendless=T, ylim=NULL, paper_style=T, angle_labels=F, smf=NULL, genotype_object=NULL, variant=NULL, smf_sample_column='phenotype_id', smf_genotype_column='genotype_id', better_colors=T) {
  # extract the region data
  region_data <- accessibility_table[accessibility_table[[region_column]] == region, -c(..region_column)]
  # make into table
  region_plot_df <- data.table('sample' = names(region_data), 'accessibility' = as.vector(unlist(region_data)))
  # extract the expression data
  gene_data <- expression_table[expression_table[[gene_column]] == gene, -c(..gene_column)]
  # make into table
  gene_plot_df <- data.table('sample' = names(gene_data), 'expression' = as.vector(unlist(gene_data)))
  # merge these
  plot_df <- merge(x = region_plot_df, y = gene_plot_df, by = 'sample')
  # we'll initialize the plot
  p <- NULL
  # we'll add the genotype info if we have it
  if (!is.null(smf) & !is.null(genotype_object)) {
    # add the genotype id to the plot
    plot_df[['donor']] <- smf[match(plot_df[['sample']], smf[[smf_sample_column]]), ][[smf_genotype_column]]
    # get the genotypes for the donors
    genotypes_df <- data.frame(genotypes$genotypes[intersect(plot_df[['donor']], rownames(genotypes$genotypes)), variant])
    # add the donor as explicit column
    genotypes_df[['donor']] <- rownames(genotypes_df)
    # rename columns
    colnames(genotypes_df) <- c('genotype', 'donor')
    # and make the genotype a character string
    genotypes_df[['genotype']] <- as.character(genotypes_df[['genotype']])
    # merge that onto the plot df
    plot_df <- merge(plot_df, genotypes_df, by = 'donor')
    # plot these against one another
    p <- ggplot(data = plot_df, mapping = aes(x = accessibility, y = expression, col = genotype))
  }
  else {
    # plot these against one another
    p <- ggplot(data = plot_df, mapping = aes(x = accessibility, y = expression))
  }
  p <- p +
    geom_point() +
    geom_smooth(method='lm', formula = y ~ x) +
    xlab(paste(region, 'accessibility')) +
    ylab(paste(gene, 'expression'))
  if(pointless){
      p <- p + theme(axis.text.x=element_blank(), 
                     axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  if (better_colors & 'genotype' %in% colnames(plot_df)) {
    p <- p + scale_color_manual(values = roycols::get_color_list(genotypes_df[['genotype']]))
  }
  return(p)
}


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')
set.seed(7777)


#####################
# Get genotype data #
#####################

# location of the genotype data
genotype_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/genotype/imputed_hg38_all_anc'
# read the genotypes
genotypes <- read.plink(
  bed = paste(genotype_loc, '.bed', sep = ''),
  bim = paste(genotype_loc, '.bim', sep = ''),
  fam = paste(genotype_loc, '.fam', sep = '')
)
# location of the smfs
smf_eqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/input/L1/combined/smf.txt'
smf_caqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/smf.txt'
# read them both
smf_eqtl <- fread(smf_eqtl_loc, header = T, sep = '\t')
smf_caqtl <- fread(smf_caqtl_loc, header = T, sep = '\t')
# merge them
smf <- rbind(smf_eqtl, smf_caqtl)
# and make unique
smf <- unique(smf)
# turn double semicolon into double dot
smf[['phenotype_id']] <- gsub(';;', '..', smf[['phenotype_id']])


########################
# Add layers to Seurat #
########################

# location of the Seurat objects
mono_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_monocyte_1_80_20240521.rds'
# read the object
mono_object <- readRDS(mono_object_loc)
# add the binarized assay
mono_object <- add_binarized_assay(mono_object)
# add pflogpf
mono_object <- normalize_mj(mono_object)


####################
# Read CRE outputs #
####################

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# location of the pseudobulk CRE mapping
pseudobulk_output_folder <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/eQTA_v2/L1/'
# location of the binomial method
binomial_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/combined/betas_ps/'

# read the pseudobulk outputs
pseudobulk_output_ut_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, 'UT', sep = '/'), add_mtc = T, filter_alpha = T, add_global_nominal_threshold = T, add_local_nominal_threshold = T)
pseudobulk_output_24hca_list <- read_pseudobulk_cre_output_per_celltype(paste(pseudobulk_output_folder, '24hCA', sep = '/'), add_mtc = T, filter_alpha = T, add_global_nominal_threshold = T, add_local_nominal_threshold = T)
# merge the cell types
pseudobulk_output_ut <- do.call('rbind', pseudobulk_output_ut_list)
pseudobulk_output_24hca <- do.call('rbind', pseudobulk_output_24hca_list)
# add the condition
pseudobulk_output_ut[['condition']] <- 'UT'
pseudobulk_output_24hca[['condition']] <- '24hCA'
# add z score
pseudobulk_output_ut[['zscore']] <- pseudobulk_output_ut[['beta']] / pseudobulk_output_ut[['beta_se']]
pseudobulk_output_24hca[['zscore']] <- pseudobulk_output_24hca[['beta']] / pseudobulk_output_24hca[['beta_se']]
# add a correlation based on the Z score, taking sample size and removing 2 + 10 PCs to get the degrees of freedom
pseudobulk_output_ut[['r']] <- pseudobulk_output_ut[['zscore']] / sqrt(pseudobulk_output_ut[['zscore']]^2 + (pseudobulk_output_ut[['n_samples']][1] - 12))
pseudobulk_output_24hca[['r']] <- pseudobulk_output_24hca[['zscore']] / sqrt(pseudobulk_output_24hca[['zscore']]^2 + (pseudobulk_output_24hca[['n_samples']][1] - 12))
# add a p based z
pseudobulk_output_ut[['z_from_p']] <- qnorm(1 - pseudobulk_output_ut[['p_value']] / 2) * sign(pseudobulk_output_ut[['beta']])
pseudobulk_output_24hca[['z_from_p']] <- qnorm(1 - pseudobulk_output_24hca[['p_value']] / 2) * sign(pseudobulk_output_24hca[['beta']])
pseudobulk_output_ut[['z_from_p']] <- qnorm(pseudobulk_output_ut[['p_value']] / 2) * -1 * sign(pseudobulk_output_ut[['beta']])
pseudobulk_output_24hca[['z_from_p']] <- qnorm(pseudobulk_output_24hca[['p_value']] / 2) * -1 * sign(pseudobulk_output_24hca[['beta']])
# merge them
pseudobulk_output <- do.call('rbind', list(pseudobulk_output_ut, pseudobulk_output_24hca))
# filter on theshold
pseudobulk_output_unfiltered <- pseudobulk_output
# check significance threshold
pseudobulk_output <- pseudobulk_output[
  !is.na(pseudobulk_output[['p_value']]) & !is.na(pseudobulk_output[['pval_nominal_threshold_global']]) & 
    pseudobulk_output[['p_value']] < pseudobulk_output[['pval_nominal_threshold_global']] &
    !is.na(pseudobulk_output[['feature_q_value']]) & pseudobulk_output[['feature_q_value']] < 0.05 , ]

# read the binomial results
binomial_output_list <- read_binomial_output_per_celltype(binomial_output_loc)
# merge them
binomial_output <- do.call('rbind', binomial_output_list)
# add a correlation based on the Z score, taking sample size and removing 2 + 10 PCs to get the degrees of freedom
binomial_output[['r']] <- binomial_output[['meta_z']] / sqrt(binomial_output[['meta_z']]^2 + (binomial_output[['n_sample']][1] - 2))

# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')


#############################
# Read binomial model input #
#############################

# the locations of the files
mono_betas_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/combined/betas_ps/monocyte/betas.tsv.gz'
mono_ses_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/combined/betas_ps/monocyte/ses.tsv.gz'
# read the files
mono_betas <- fread(mono_betas_loc, header = T, sep = '\t', na.strings = c('NA', 'NaN', 'nan'))
mono_ses <- fread(mono_ses_loc, header = T, sep = '\t', na.strings = c('NA', 'NaN', 'nan'))
# calculate the z scores
mono_zs <- mono_betas[, -c('region', 'gene')] / mono_ses[, -c('region', 'gene')]
# add region and gene back
mono_zs <- cbind(mono_betas[, c('region', 'gene')], mono_zs)


###############################
# Read pseudobulk model input #
###############################

# the location of the files
mono_ut_eqtl_input_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/L1/UT/monocyte.qtlInput.txt.gz'
mono_24hca_eqtl_input_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/L1/24hCA/monocyte.qtlInput.txt.gz'
mono_ut_caqtl_input_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/monocyte.qtlInput.txt.gz'
mono_24hca_caqtl_input_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/24hCA/monocyte.qtlInput.txt.gz'
# read all of these
mono_ut_eqtl_input <- read.table(mono_ut_eqtl_input_loc, header = T, sep = '\t', row.names = 1)
mono_24hca_eqtl_input <- read.table(mono_24hca_eqtl_input_loc, header = T, sep = '\t', row.names = 1)
mono_ut_caqtl_input <- read.table(mono_ut_caqtl_input_loc, header = T, sep = '\t', row.names = 1)
mono_24hca_caqtl_input <- read.table(mono_24hca_caqtl_input_loc, header = T, sep = '\t', row.names = 1)
# add the features as explicit columns
mono_ut_eqtl_input <- cbind('feature_id' = rownames(mono_ut_eqtl_input), mono_ut_eqtl_input)
mono_24hca_eqtl_input <- cbind('feature_id' = rownames(mono_24hca_eqtl_input), mono_24hca_eqtl_input)
mono_ut_caqtl_input <- cbind('feature_id' = rownames(mono_ut_caqtl_input), mono_ut_caqtl_input)
mono_24hca_caqtl_input <- cbind('feature_id' = rownames(mono_24hca_caqtl_input), mono_24hca_caqtl_input)
# that way they can be safely converted to data.table
mono_ut_eqtl_input <- data.table(mono_ut_eqtl_input)
mono_24hca_eqtl_input <- data.table(mono_24hca_eqtl_input)
mono_ut_caqtl_input <- data.table(mono_ut_caqtl_input)
mono_24hca_caqtl_input <- data.table(mono_24hca_caqtl_input)


#########################
# Plot opposite effects #
#########################

# plot an effect that is opposite in pseudobulk vs binomial as single cell
plot_binarized_vs_other(mono_object, 'chr17-45585724-45587049', 'AC126544.2')
# plot the pseudobulk effect
plot_grid(
  plot_pseudobulk_cre(mono_ut_caqtl_input, mono_ut_eqtl_input, 'chr17-45585724-45587049', 'AC126544.2') + ggtitle('UT'), 
  plot_pseudobulk_cre(mono_24hca_caqtl_input, mono_24hca_eqtl_input, 'chr17-45585724-45587049', 'AC126544.2') + ggtitle('24hCA'), 
  nrow = 1, 
  ncol = 2
)
# or with genotype
plot_grid(
  plot_pseudobulk_cre(mono_ut_caqtl_input, mono_ut_eqtl_input, 'chr17-45585724-45587049', 'AC126544.2', smf = smf, genotype_object = genotypes, variant = '17:46108697:C:T') + ggtitle('UT'), 
  plot_pseudobulk_cre(mono_24hca_caqtl_input, mono_24hca_eqtl_input, 'chr17-45585724-45587049', 'AC126544.2', smf = smf, genotype_object = genotypes, variant = '17:46108697:C:T') + ggtitle('24hCA'), 
  nrow = 1, 
  ncol = 2
)
# plot the z scores from the per-sample binomial model
plot_bino_z_scores(mono_zs, 'chr17-45585724-45587049', 'AC126544.2')
# and their distribution
plot_bino_z_distribution(mono_zs, 'chr17-45585724-45587049', 'AC126544.2')
