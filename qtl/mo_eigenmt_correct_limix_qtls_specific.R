#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_eigenmt_correct_limix_qtls_specific.R
# Function: perform eigenMT MTC on QTL summary stats coming from LIMIX-QTL
############################################################################################################################


####################
# libraries        #
####################

library(snpStats)
library(data.table)
library(r2r)
library(rhdf5)
library(ReigenMT)
library(qvalue)


####################
# Functions        #
####################


#' Get PLINK Genotypes for a Chromosome
#'
#' This function reads PLINK genotype files for a specified chromosome.
#'
#' @param genotype_loc Character. The location of the genotype files.
#' @param chromosome Integer. The chromosome number.
#' @param genotype_prepend Character. The prefix for the genotype files. Default is 'EUR_imputed_hg38_varFiltered_chr'.
#' @param genotype_append Character. The suffix for the genotype files. Default is ''.
#' @return A list containing the genotype data.
#' @export
#' @examples
#' genotypes <- get_plink_genotypes_chromosome("path/to/genotypes", 1)
get_plink_genotypes_chromosome <- function(genotype_loc, chromosome, genotype_prepend='EUR_imputed_hg38_varFiltered_chr', genotype_append='', variants_to_load=NULL) {
  # get full path
  genotypes_loc_full <- paste(genotype_loc, genotype_prepend, chromosome, genotype_append, sep = '')
  # initialize variable
  genotypes <- NULL
  # check if we need to intersect first
  if (!is.null(variants_to_load)){
    # get variants present
    variants_present <- fread(paste(genotypes_loc_full, 'bim', sep = '.'))[[2]]
    # check which are overlapping
    variants_overlapping <- intersect(variants_to_load, variants_present)
    # check if we are missing anything
    if (length(variants_to_load) > length(variants_overlapping)) {
      warning(paste('asking for more variants than present in genotype data!', as.character(length(variants_to_load) - length(variants_overlapping)), 'of', as.character(length(variants_to_load)), 'are missing!'))
    }
    genotypes <- read.plink(
      bed = paste(genotypes_loc_full, '.bed', sep = ''),
      bim = paste(genotypes_loc_full, '.bim', sep = ''),
      fam = paste(genotypes_loc_full, '.fam', sep = ''), 
      select.snps = variants_overlapping
    )
  } else {
    genotypes <- read.plink(
      bed = paste(genotypes_loc_full, '.bed', sep = ''),
      bim = paste(genotypes_loc_full, '.bim', sep = ''),
      fam = paste(genotypes_loc_full, '.fam', sep = '')
    )
  }
  
  return(genotypes)
}

#' Correct All Chunks for a Chromosome
#'
#' This function corrects all chunks for a specified chromosome using eigenMT.
#'
#' @param genotypes_chromosome List. The genotype data for the chromosome.
#' @param var_locations_chromosome Data frame. The variant locations for the chromosome.
#' @param chunks_loc Character. The location of the chunk files.
#' @param chunk_pattern Character. The pattern to match chunk files.
#' @param pvalue_column Character. The name of the p-value column. Default is 'p_value'.
#' @return A data frame containing the corrected chunks.
#' @export
#' @examples
#' corrected_chunks <- correct_all_chunks_chromosome(genotypes, var_locations, "path/to/chunks", "chunk_pattern")
correct_all_chunks_chromosome <- function(genotypes_chromosome, var_locations_chromosome, chunks_loc, chunk_pattern, pvalue_column='p_value') {
  # list all files
  chunk_files <- list.files(chunks_loc)
  # filter by regular expression
  chunk_files <- chunk_files[grepl(chunk_pattern, chunk_files)]
  # save in a list
  all_chunks <- list()
  # check each file
  for (chunk in chunk_files) {
    # load summary statistics
    summary_stats_h5 <- ReigenMT::limix_h5_to_sumstats_format(paste(chunks_loc, chunk, sep = '/'))
    # perform eigenMT correction
    eigen_corrected_chunk <- ReigenMT::eigenmt(summary_stats = summary_stats_h5, genotypes = genotypes_chromosome, genotype_to_position = var_locations_chromosome, variant_column_summary_stats = 'snp_id', feature_column_summary_stats = 'feature', var_explained_threshold = 0.975, pvalue_column = pvalue_column)
    # put in the list
    all_chunks[[chunk]] <- eigen_corrected_chunk
  }
  # merge the corrected chunks
  chunks_merged <- do.call('rbind', all_chunks)
  return(chunks_merged)
}

#' Correct All Chunks for All Chromosomes
#'
#' This function corrects all chunks for all specified chromosomes using eigenMT.
#'
#' @param genotype_loc Character. The location of the genotype files.
#' @param chunks_loc Character. The location of the chunk files.
#' @param chromosomes Integer vector. The chromosomes to process. Default is 1:22.
#' @param genotype_prepend Character. The prefix for the genotype files. Default is 'EUR_imputed_hg38_varFiltered_chr'.
#' @param genotype_append Character. The suffix for the genotype files. Default is ''.
#' @param pvalue_column Character. The name of the p-value column. Default is 'p_value'.
#' @param qtl_results_prepend Character. The prefix for the QTL results files. Default is 'qtl_results_'.
#' @return A data frame containing the corrected chunks for all chromosomes.
#' @export
#' @examples
#' corrected_all <- correct_all_chunks_all_chromosomes("path/to/genotypes", "path/to/chunks")
correct_all_chunks_all_chromosomes <- function(genotype_loc, chunks_loc, chromosomes=1:22, genotype_prepend='EUR_imputed_hg38_varFiltered_chr', genotype_append='', pvalue_column='p_value', qtl_results_prepend='qtl_results_', verbose=T) {
  all_chrom_chunks <- list()
  # check each chromosome passed
  for (chrom in chromosomes) {
    if (verbose) {
      message(paste('loading chromosome', chrom, 'genotype data\n'))
    }
    # lead the genotype data in plink format
    genotypes_chromosome <- get_plink_genotypes_chromosome(genotype_loc, chrom, genotype_prepend = genotype_prepend, genotype_append = genotype_append)
    if (verbose) {
      message(paste('loading chromosome', chrom, 'variant positions\n'))
    }
    # extract the locations of the variants from the genotype data
    var_locations_chromosome <- get_position_hashtab(genotypes_chromosome)
    # create the pattern for this chromosome, to list the h5 files
    h5_pattern_chromosome <- paste(qtl_results_prepend, chrom, '_\\d+_\\d+.h5', sep = '')
    if (verbose) {
      message(paste('processing chromosome', chrom, 'QTL chunks\n'))
    }
    # do eigenMT MTC for this chunk
    chunks_chromosome <- correct_all_chunks_chromosome(
      genotypes_chromosome = genotypes_chromosome, 
      var_locations_chromosome = var_locations_chromosome, 
      chunks_loc = chunks_loc, 
      chunk_pattern = h5_pattern_chromosome, 
      pvalue_column = pvalue_column
    )
    # put in the list
    all_chrom_chunks[[chrom]] <- chunks_chromosome
  }
  # merge the chunks
  all_chrom_chunks_merged <- do.call('rbind', all_chrom_chunks)
  return(all_chrom_chunks_merged)
}


perform_qvalue_correction <- function(cell_type_output, mtc_column='empirical_feature_p_value', feature_mtc_column='feature_id', mtc_column_to_add='feature_q_value') {
  # subset to only the important columns
  cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column), with = F]
  # order by significance
  cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
  # remove where the feature is smaller than zero
  cell_type_output_features <- cell_type_output_features[!(cell_type_output_features[[mtc_column]] < 0), ]
  # keep only the first entry
  cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
  # set the values that are larger than 1, to be 1, problem with precision
  cell_type_output_features[cell_type_output_features[[mtc_column]] > 1, mtc_column] <- 1
  # add multiple testing correction
  cell_type_output_features[['qvalue']] <- qvalue::qvalue(cell_type_output_features[[mtc_column]])$qvalues
  # now add back to the original table
  cell_type_output[[mtc_column_to_add]] <- cell_type_output_features[match(cell_type_output[[feature_mtc_column]], cell_type_output_features[[feature_mtc_column]]), 'qvalue'][['qvalue']]
  return(cell_type_output)
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


perform_nominal_threshold_calculation <- function(cell_type_output, significance_cutoff=0.05, mtc_column='feature_q_value', feature_mtc_column='feature_id', nominal_threshold_column='pval_nominal_threshold_global', alpha_column='alpha_param', beta_column='beta_param', nominal_p_column='p_value') {
  # subset to only the important columns
  cell_type_output_features <- cell_type_output[, c(feature_mtc_column, mtc_column, alpha_column, beta_column, nominal_p_column), with = F]
  # order by significance
  cell_type_output_features <- cell_type_output_features[order(cell_type_output_features[[mtc_column]]), ]
  # remove where the feature is smaller than zero
  cell_type_output_features <- cell_type_output_features[!(cell_type_output_features[[mtc_column]] < 0), ]
  # keep only the first entry
  cell_type_output_features[!duplicated(cell_type_output_features[[feature_mtc_column]]), ]
  # get the nominal p value cutoff based on the p values and the beta distribution
  cell_type_output_global_threshold <- calculate_nominal_thresholds(cell_type_output_features, fdr=significance_cutoff, pval_col=nominal_p_column, nominal_threshold_column='nomthres', cutoff_column='qvalue', alpha_column = alpha_column, beta_column = beta_column)
  # now add the nominal threshold to the full table
  cell_type_output[[nominal_threshold_column]] <- cell_type_output_global_threshold[match(cell_type_output[[feature_mtc_column]], cell_type_output_global_threshold[[feature_mtc_column]]), ][['nomthres']]
  return(cell_type_output)
}

correct_feature <- function(genotypes_chromosome, var_locations_chromosome, summary_stats, pvalue_column='p_value', variant_column_summary_stats = 'snp_id',  feature_column_summary_stats='feature', var_explained_threshold = 0.975) {
  # perform eigenMT correction
  eigen_corrected_chunk <- ReigenMT::eigenmt(summary_stats = summary_stats, genotypes = genotypes_chromosome, genotype_to_position = var_locations_chromosome, variant_column_summary_stats = variant_column_summary_stats, feature_column_summary_stats = feature_column_summary_stats, var_explained_threshold = var_explained_threshold, pvalue_column = pvalue_column)
  return(eigen_corrected_chunk)
}



correct_all_chunks_all_chromosomes_from_single_file <- function(genotype_loc, qtl_data, chromosomes=1:22, genotype_prepend='EUR_imputed_hg38_varFiltered_chr', genotype_append='', pvalue_column='p_value', feature_column='feature_id', chrom_column_feature='feature_chromosome', variant_column='snp_id', verbose=T) {
  # store result per chromosome
  all_chrom_chunks <- list()
  # check each chromosome passed
  for (chrom in chromosomes) {
    # subset to chromosome
    qtl_data_chrom <- qtl_data[
      qtl_data[[chrom_column_feature]] == chrom, 
    ]
    # get the variants tested for this chromosome
    variants_chrom <- unique(qtl_data_chrom[[variant_column]])
    # say what we are doing
    if (verbose) {
      message(paste('loading chromosome', chrom, 'genotype data\n'))
    }
    # lead the genotype data in plink format
    genotypes_chromosome <- get_plink_genotypes_chromosome(genotype_loc, chrom, genotype_prepend = genotype_prepend, genotype_append = genotype_append, variants_to_load = variants_chrom)
    if (verbose) {
      message(paste('loading chromosome', chrom, 'variant positions\n'))
    }
    # extract the locations of the variants from the genotype data
    var_locations_chromosome <- get_position_hashtab(genotypes_chromosome)
    # now check each feature
    for (feature in unique(qtl_data_chrom[[feature_column]])) {
      # subset to that feature
      qtl_data_feature <- qtl_data_chrom[qtl_data_chrom[[feature_column]] == feature, ]
      # do eigenMT MTC for this chunk
      chunks_chromosome <- correct_feature(
        genotypes_chromosome = genotypes_chromosome, 
        var_locations_chromosome = var_locations_chromosome, 
        summary_stats = qtl_data_feature,
        pvalue_column = pvalue_column, 
        variant_column_summary_stats = variant_column,  
        feature_column_summary_stats= feature_column
      )
      # put in the list
      all_chrom_chunks[[feature]] <- chunks_chromosome
    }
  }
  # merge the chunks
  all_chrom_chunks_merged <- do.call('rbind', all_chrom_chunks)
  return(all_chrom_chunks_merged)
}


####################
# Settings        #
####################


####################
# Main Code        #
####################

# genotypes 
genotypes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/'

# location of the eQTL interactions
eqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'

# check each cell type in the eQTLs
eqtls_per_celltype <- list()
for (ct in list.dirs(eqtl_loc, recursive = F, full.names = F)) {
# for (ct in c('CD4T', 'CD8T', 'DC', 'monocyte', 'NK')) {
  print(ct)
  qtl_data_loc <- paste(eqtl_loc, '/', ct, '/qtl_results_all.txt.gz', sep = '')
  qtl_data <- fread(qtl_data_loc, header = T, sep = '\t')
  ct_res <- correct_all_chunks_all_chromosomes_from_single_file(
    genotype_loc = genotypes_loc, 
    qtl_data = qtl_data, 
    chromosomes = c(8), 
    genotype_prepend = 'EUR_imputed_hg38_varFiltered_chr', 
    genotype_append = '', 
    pvalue_column = 'p_value'
  )
  # set output loc
  qtl_data_out_loc <- paste(eqtl_loc, '/', ct, '/qtl_results_eigenmt_chr8.tsv.gz', sep = '')
  # store the result
  write.table(ct_res, gzfile(qtl_data_out_loc), row.names = F, col.names = T, sep = '\t')
  # make checksum
  mdfiver::create_sha256_for_file(qtl_data_out_loc)
}
