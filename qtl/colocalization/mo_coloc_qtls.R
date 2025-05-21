#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_coloc_qtls.R
# Function: colocalize the eQTLs with the caQTLs
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(susieR)
library(Rfast) # not in container
library(coloc)

####################
# Functions        #
####################


calculate_distances_per_chromosome <- function(region_annotation1, region_annotation2, annotation1_chrom_column='chromosome', annotation2_chrom_column='chromosome', annotation1_name_column='feature_id', annotation2_name_column='feature_id', annotation1_start_column='start', annotation2_start_column='start', annotation1_end_column='end', annotation2_end_column='end') {
  # get the chromosomes from the first table
  chroms_ann1 <- unique(region_annotation1[[annotation1_chrom_column]])
  # and the ones from the second table
  chroms_ann2 <- unique(region_annotation2[[annotation2_chrom_column]])
  # only keep the chromosomes we have in both tables
  chroms_both <- intersect(chroms_ann1, chroms_ann2)
  # we'll first store the results per chromosome
  distances_per_chrom <- list()
  # and check each chromosome
  for (chrom in chroms_both) {
    # subset both of the tables on that chromosome
    ann1_chrom <- region_annotation1[!is.na(region_annotation1[[annotation1_chrom_column]]) & region_annotation1[[annotation1_chrom_column]] == chrom, ]
    ann2_chrom <- region_annotation2[!is.na(region_annotation2[[annotation2_chrom_column]]) & region_annotation2[[annotation2_chrom_column]] == chrom, ]
    # extract the regions coordinates and metadata
    ann1_chrom <- ann1_chrom[, c(..annotation1_chrom_column, ..annotation1_name_column, ..annotation1_start_column, ..annotation1_end_column)]
    ann2_chrom <- ann2_chrom[, c(..annotation2_name_column, ..annotation2_start_column, ..annotation2_end_column)]
    # rename columns so they are unique
    colnames(ann1_chrom) <- c('chromosome', 'name_region1', 'start_region1', 'end_region1')
    colnames(ann2_chrom) <- c('name_region2', 'start_region2', 'end_region2')
    # make every possible combination
    ann_both_chrom <- data.table(merge(data.frame(ann1_chrom), data.frame(ann2_chrom), by = NULL))
    # now calculate the distances
    ann_both_chrom[['start1_to_start2']] <- ann_both_chrom[['start_region1']] - ann_both_chrom[['start_region2']]
    ann_both_chrom[['start1_to_end2']] <- ann_both_chrom[['start_region1']] - ann_both_chrom[['end_region2']]
    ann_both_chrom[['end1_to_start2']] <- ann_both_chrom[['end_region1']] - ann_both_chrom[['start_region2']]
    ann_both_chrom[['end1_to_end2']] <- ann_both_chrom[['end_region1']] - ann_both_chrom[['end_region2']]
    # and put the result in the list
    distances_per_chrom[[as.character(chrom)]] <- ann_both_chrom
  }
  # merge the results per chromosome
  distances_all <- do.call('rbind', distances_per_chrom)
  return(distances_all)
}


#' Colocalization of Fine-Mapping Datasets
#'
#' This function performs colocalization analysis on two fine-mapping datasets using Bayes factors.
#'
#' @param dataset1_fm_table A data frame containing the fine-mapping results for the first dataset.
#' @param dataset2_fm_table A data frame containing the fine-mapping results for the second dataset.
#' @param dataset1_name A character string specifying the name of the first dataset. Default is 'dataset1'.
#' @param dataset2_name A character string specifying the name of the second dataset. Default is 'dataset2'.
#' @param feature1_column A character string specifying the column name for features in the first dataset. Default is 'feature'.
#' @param feature2_column A character string specifying the column name for features in the second dataset. Default is 'feature'.
#' @param snp1_column A character string specifying the column name for SNPs in the first dataset. Default is 'variant'.
#' @param snp2_column A character string specifying the column name for SNPs in the second dataset. Default is 'variant'.
#' @param cs_columns1 A vector of column indices or names specifying the credible set columns in the first dataset. Default is columns matching '^CS\\d+$'.
#' @param cs_columns2 A vector of column indices or names specifying the credible set columns in the second dataset. Default is columns matching '^CS\\d+$'.
#'
#' @return A list of results for each feature combination that is present in the datasets. Each element of the list contains:
#' \item{feature1}{The feature being analyzed from the first list.}
#' \item{feature2}{The feature being analyzed from the second list.}
#' \item{dataset1}{The name of the first dataset.}
#' \item{dataset2}{The name of the second dataset.}
#' \item{variants1}{The variants from the first dataset.}
#' \item{variants2}{The variants from the second dataset.}
#' \item{result}{The result of the colocalization analysis.}
#' \item{table}{A data frame summarizing the colocalization results.}
#'
#' @examples
#' \dontrun{
#' dataset1 <- data.frame(feature = c('gene1', 'gene1', 'gene2'),
#'                        variant = c('rs1', 'rs2', 'rs3'),
#'                        CS1 = c(0.8, 0.1, 0.6),
#'                        CS2 = c(0.2, 0.9, 0.4))
#' dataset2 <- data.frame(feature = c('gene1', 'gene2', 'gene2'),
#'                        variant = c('rs1', 'rs3', 'rs4'),
#'                        CS1 = c(0.7, 0.5, 0.3),
#'                        CS2 = c(0.3, 0.5, 0.7))
#' result <- coloc_datasets(dataset1, dataset2)
#' }
#'
coloc_datasets <- function(dataset1_fm_table, 
                           dataset2_fm_table, 
                           dataset1_name='dataset1', 
                           dataset2_name='dataset2', 
                           feature1_column='feature', 
                           feature2_column='feature', 
                           snp1_column='variant', 
                           snp2_column='variant', 
                           cs_columns1=grep('^CS\\d+$', colnames(dataset1_fm_table)), 
                           cs_columns2=grep('^CS\\d+$', colnames(dataset2_fm_table))) {
  # extract the features for table 1
  features_table1 <- unique(dataset1_fm_table[[feature1_column]])
  # and table 2
  features_table2 <- unique(dataset2_fm_table[[feature2_column]])
  # we'll store results per feature
  coloc_per_feature_combination <- list()
  # now do each feature in dataset 1
  for (feature1 in features_table1) {
    # against each feature in dataset 2
    for (feature2 in features_table2) {
      # subset each table to those features
      table1_feature <- dataset1_fm_table[!is.na(dataset1_fm_table[[feature1_column]]) &
                                            dataset1_fm_table[[feature1_column]] == feature1, ]
      table2_feature <- dataset2_fm_table[!is.na(dataset2_fm_table[[feature2_column]]) &
                                            dataset2_fm_table[[feature2_column]] == feature2, ]
      # extract the variants
      table1_feature_variants <- table1_feature[[snp1_column]]
      table2_feature_variants <- table2_feature[[snp2_column]]
      # get just the credible set information
      table1_cs <- table1_feature[, cs_columns1]
      table2_cs <- table2_feature[, cs_columns2]
      # remove completely NA columns, these can be there due to padding when the finemapped results were aggretated into a single table
      table1_cs <- table1_cs[, !apply(table1_cs, 2, function(col) all(is.na(col)))]
      table2_cs <- table2_cs[, !apply(table2_cs, 2, function(col) all(is.na(col)))]
      # transpose them
      table1_cs_t <- t(table1_cs)
      table2_cs_t <- t(table2_cs)
      # and set the variants as the column names
      colnames(table1_cs_t) <- table1_feature_variants
      colnames(table2_cs_t) <- table2_feature_variants
      # finally do the actual coloc
      coloc_bfbf <- coloc.bf_bf(table1_cs_t, table2_cs_t, overlap.min = 0.1)
      # extract the result table
      res_table <- as.data.frame(coloc_bfbf[['summary']])
      # add the variants
      res_table[['variant1']] <- table1_feature_variants[res_table[['idx1']]]
      res_table[['variant2']] <- table2_feature_variants[res_table[['idx2']]]
      # now add the features and the two datasets
      res_table <- cbind(
        data.frame(
          'dataset1' = rep(dataset1_name, times = nrow(res_table)), 
          'dataset2' = rep(dataset2_name, times = nrow(res_table)), 
          'trait1' = rep(feature1, times = nrow(res_table)), 
          'trait2' = rep(feature2, times = nrow(res_table))),
        res_table
      )
      # make a list with the results
      result_list <- list(
        'feature1' = feature1,
        'feature2' = feature2,
        'dataset1' = dataset1_name, 
        'dataset2' = dataset2_name, 
        'variants1' = table1_feature_variants, 
        'variants2' = table2_feature_variants, 
        'sets1' = table1_cs_t, 
        'sets2' = table2_cs_t, 
        'result' = coloc_bfbf, 
        'table' = res_table
      )
      # and put into the bigger list
      coloc_per_feature_combination[[paste0(feature1, '_vs_', feature2)]] <- result_list
    }
  }
  return(coloc_per_feature_combination)
}


####################
# Settings        #
####################

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# location of the finemapped eQTLs
finemapped_eqtls_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# location of the finemapped caQTLs
finemapped_caqtls_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_caqtl/sc-eqtlgen/output/combined_significant/L1/'

# location of gene annotations
gene_annotation_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/input/LimixAnnotationFile.txt'
# location of the ATAC regions
region_annotation_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/annotations/pct01/LimixAnnotationFile.tsv.gz'

# where we will save the results
coloc_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/colocalization/eqtl_caqtl/'

# cell types to consider
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

# prepend of the eqtls
eqtl_finemapping_prepend <- ''
# append of the eQTLs
eqtl_finemapping_append <- '_finemapped.tsv.gz'
# prepend of the caQTLs
caqtl_finemapping_prepend <- ''
# append of the caQTLs
caqtl_finemapping_append <- '_finemapped.tsv.gz'
# prepend of coloc
coloc_prepend <- ''
coloc_append <- '_eqtl_caqtl_coloc.rds'
coloc_tbl_append <- '_eqtl_caqtl_coloc.tsv.gz'

# name of the first data
dataset1_name <- 'caQTLs'
dataset2_name <- 'eQTLs'

# read the annotation files
gene_annotation <- fread(gene_annotation_loc, header = T, sep = '\t')
region_annotation <- fread(region_annotation_loc, header = T, sep = '\t')

# # get the distances between the genes and the regions
# region_gene_distances <- calculate_distances_per_chromosome(region_annotation, gene_annotation)

# check each cell type
for (cell_type in cell_types) {
  # paste together the paths
  eqtl_path <- paste0(finemapped_eqtls_loc, '/', eqtl_finemapping_prepend, cell_type, eqtl_finemapping_append)
  caqtl_path <- paste0(finemapped_caqtls_loc, '/', caqtl_finemapping_prepend, cell_type, caqtl_finemapping_append)
  # read these files
  eqtls <- fread(eqtl_path, header = T, sep = '\t')
  caqtls <- fread(caqtl_path, header = T, sep = '\t')
  # store result per chromosome
  fm_per_chrom <- list()
  # also with tables
  fm_tbl_per_chrom <- list()
  # check each of the chromosomes
  # for (chrom in unique(region_gene_distances[['chromosome']])) {
  for (chrom in unique(gene_annotation[['chromosome']])) {
    # subset to the annotations for that chromosome
    # region_gene_distances_chrom <- region_gene_distances[region_gene_distances[['chromosome']] == chrom
    #                                                      !is.na(region_gene_distances[['chromosome']]), ]
    # # get the regions for this chromosome
    # genes_chrom <- unique(region_gene_distances_chrom[['name_region1']])
    # peaks_chrom <- unique(region_gene_distances_chrom[['name_region2']])
    genes_chrom <- gene_annotation[gene_annotation[['chromosome']] == chrom, ][['feature_id']]
    peaks_chrom <- region_annotation[region_annotation[['chromosome']] == chrom, ][['feature_id']]
    # subset the finemapping outputs
    eqtls_chrom <- eqtls[eqtls[['feature']] %in% genes_chrom, ]
    caqtls_chrom <- caqtls[caqtls[['feature']] %in% peaks_chrom, ]
    # perform the coloc
    coloc_chrom <- coloc_datasets(
      data.frame(caqtls_chrom), 
      data.frame(eqtls_chrom), 
      dataset1_name = dataset1_name, 
      dataset2_name = dataset2_name
    )
    # put in the list
    fm_per_chrom[[as.character(chrom)]] <- coloc_chrom
    # put in list of table
    for (coloc_traits in names(coloc_chrom)) {
      # but only if there are multiple rows, and as such a result
      if(nrow(coloc_chrom[[coloc_traits]][['table']] > 0)) {
        fm_tbl_per_chrom[[coloc_traits]] <- coloc_chrom[[coloc_traits]][['table']]
      }
    }
  }
  # paste together the output location
  coloc_ct_output_loc <- paste0(coloc_output_loc, '/', coloc_prepend, cell_type, coloc_append)
  # write the output
  saveRDS(fm_per_chrom, coloc_ct_output_loc)
  # make a checksum
  mdfiver::create_md5_for_file(coloc_ct_output_loc)
  # merge the tables
  fm_tbl_all_chrom <- do.call('rbind', fm_tbl_per_chrom)
  # write the table
  coloc_tbl_ct_output_loc <- paste0(coloc_output_loc, '/', coloc_prepend, cell_type, coloc_tbl_append)
  write.table(fm_tbl_all_chrom, gzfile(coloc_tbl_ct_output_loc), row.names = F, col.names = T, sep = '\t')
  mdfiver::create_md5_for_file(coloc_tbl_ct_output_loc)
}
