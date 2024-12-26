#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_eigenmt_correct_limix_qtls.R
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
get_plink_genotypes_chromosome <- function(genotype_loc, chromosome, genotype_prepend='EUR_imputed_hg38_varFiltered_chr', genotype_append='') {
  genotypes_loc <- paste(genotype_loc, genotype_prepend, chromosome, genotype_append, sep = '')
  genotypes <- read.plink(
    bed = paste(genotypes_loc, '.bed', sep = ''),
    bim = paste(genotypes_loc, '.bim', sep = ''),
    fam = paste(genotypes_loc, '.fam', sep = '')
  )
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
  chunk_files <- list.files(chunks_loc)
  chunk_files <- chunk_files[grepl(chunk_pattern, chunk_files)]
  all_chunks <- list()
  for (chunk in chunk_files) {
    summary_stats_h5 <- ReigenMT::limix_h5_to_sumstats_format(paste(chunks_loc, chunk, sep = '/'))
    eigen_corrected_chunk <- eigenmt(summary_stats = summary_stats_h5, genotypes = genotypes_chromosome, genotype_to_position = var_locations_chromosome, variant_column_summary_stats = 'snp_id', feature_column_summary_stats = 'feature', var_explained_threshold = 0.975, pvalue_column = pvalue_column)
    all_chunks[[chunk]] <- eigen_corrected_chunk
  }
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
correct_all_chunks_all_chromosomes <- function(genotype_loc, chunks_loc, chromosomes=1:22, genotype_prepend='EUR_imputed_hg38_varFiltered_chr', genotype_append='', pvalue_column='p_value', qtl_results_prepend='qtl_results_') {
  all_chrom_chunks <- list()
  for (chrom in chromosomes) {
    genotypes_chromosome <- get_plink_genotypes_chromosome(genotype_loc, chrom, genotype_prepend = genotype_prepend, genotype_append = genotype_append)
    var_locations_chromosome <- ReigenMT::get_position_hashtab(genotypes_chromosome)
    h5_pattern_chromosome <- paste(qtl_results_prepend, chrom, '_\\d+_\\d+.h5', sep = '')
    chunks_chromosome <- correct_all_chunks_chromosome(
      genotypes_chromosome = genotypes_chromosome,
      var_locations_chromosome = var_locations_chromosome,
      chunks_loc = chunks_loc,
      chunk_pattern = h5_pattern_chromosome,
      pvalue_column = pvalue_column
    )
    all_chrom_chunks[[chrom]] <- chunks_chromosome
  }
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
genotypes_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/'
summary_stats_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/'

# location of the eQTL interactions
eqtl_interaction_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/L1/'
# location of the caQTL interactions
caqtl_interaction_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/output/L1/'

# check each cell type in the eQTLs
eqtl_interactions_per_celltype <- list()
for (ct in list.dirs(eqtl_interaction_loc, recursive = F)) {
  ct_res <- correct_all_chunks_all_chromosomes(
    genotype_loc = genotypes_loc,
    chunks_loc = paste(eqtl_interaction_loc, '/', ct, '/inflammation_final/qtl/', sep = ''),
    chromosomes = 1:22,
    genotype_prepend = 'EUR_imputed_hg38_varFiltered_chr',
    genotype_append = '',
    pvalue_column = 'p_value',
    qtl_results_prepend = 'iqtl_results_'
  )
  eqtl_interactions_per_celltype[[ct]] <- ct_res
}
