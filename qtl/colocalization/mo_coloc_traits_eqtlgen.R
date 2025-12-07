#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_coloc_traits_eqtlgen.R
# Function: perform colocalization on two finemapped traits
# Example: 
# Rscript ~/mo_coloc_traits_eqtlgen.R \
#   --dataset1_in /groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/eqtl/finemapping/output/elmentaite_adult_martin_immune/cell_type_medhigh_inflammationsplit_mincor/unconfined/pcs/CD4_T_cells_merged.tsv.gz \
#   --dataset2_in_directory /scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/eqtl-gwas-susie-coloc/output_gwas_finemap_all_20250919/GWAS_finemap/ \
#   --dataset2_in_prepend White_blood_cell_count__ \
#   --dataset2_in_append ___gwas.txt.gz \
#   --dataset1_name CD4T \
#   --dataset2_name White_blood_cell_count \
#   --output_loc /groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/eqtl/colocalizations/celltypes/elmentaite_adult_martin_immune/cell_type_medhigh_inflammationsplit_mincor/unconfined/pcs/CD4_T_cells_vs_CD8_T_cells.tsv.gz \
#   --binary_rds_loc /groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/eqtl/colocalizations/celltypes/elmentaite_adult_martin_immune/cell_type_medhigh_inflammationsplit_mincor/unconfined/pcs/CD4_T_cells_vs_CD8_T_cells.rds
#
############################################################################################################################


####################
# libraries        #
####################

library(coloc)
library(optparse)
library(mdfiver)


####################
# Functions        #
####################

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
#' @return A list of results for each feature that is present in both datasets. Each element of the list contains:
#' \item{feature}{The feature being analyzed.}
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
  # get the features we have in both
  features_both <- intersect(features_table1, features_table2)
  # we'll store results per feature
  coloc_per_feature <- list()
  # and thus check each features
  for (feature_both in features_both) {
    # subset each table to those features
    table1_feature <- dataset1_fm_table[!is.na(dataset1_fm_table[[feature1_column]]) &
                                        dataset1_fm_table[[feature1_column]] == feature_both, ]
    table2_feature <- dataset2_fm_table[!is.na(dataset2_fm_table[[feature2_column]]) &
                                        dataset2_fm_table[[feature2_column]] == feature_both, ]
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
    coloc_bfbf <- coloc.bf_bf(table1_cs_t, table2_cs_t)
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
        'trait' = rep(feature_both, times = nrow(res_table))), 
      res_table
    )
    # make a list with the results
    result_list <- list(
      'feature' = feature_both, 
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
    coloc_per_feature[[feature_both]] <- result_list
  }
  # return the bigger list
  return(coloc_per_feature)
}


####################
# Settings        #
####################

set.seed(7777)


####################
# Debug            #
####################

#do_debug()


####################
# Main Code        #
####################


# make command line options
option_list <- list(
  make_option(c("-t", "--dataset1_in"), type="character", default=NULL, 
              help="finemapping output file of first dataset", metavar="character"),
  make_option(c("-d", "--dataset2_in_directory"), type="character", default=NULL, 
              help="finemapping output file of second dataset first part of filename", metavar="character"),
  make_option(c("-p", "--dataset2_in_prepend"), type="character", default=NULL, 
              help="finemapping output file of second dataset first part of filename", metavar="character"),
  make_option(c("-s", "--dataset2_in_append"), type="character", default=NULL, 
              help="finemapping output file of second dataset first part of filename", metavar="character"),
  make_option(c("-n", "--dataset1_name"), type="character", default='dataset1', 
              help="name of the first dataset, to put in the output [default]", metavar="character"), 
  make_option(c("-a", "--dataset2_name"), type="character", default='dataset2', 
              help="name of the second dataset, to put in the output  [default]", metavar="character"),
  make_option(c("-o", "--output_loc"), type="character", default=NULL, 
              help="tab separated output location of the colocalization", metavar="character"), 
  make_option(c("-b", "--binary_rds_loc"), type="character", default=NULL, 
              help="optional output location of full binary RDS output", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize some values
dataset1_in <- NULL
dataset2_in <- NULL
output_loc <- NULL
# and these will always be set
dataset1_name <- opt[['dataset1_name']]
dataset2_name <- opt[['dataset2_name']]
# this one is optional
binary_rds_loc <- NULL

# check if we have all the relevant parameters
if (is.null(opt[['dataset1_in']])) {
  stop(paste('-t/--dataset1_in is an obligatory parameter\n'))
} else {
  dataset1_in <- opt[['dataset1_in']]
}

if (is.null(opt[['dataset2_in_directory']])) {
  stop(paste('-d/--dataset2_in_directory is an obligatory parameter\n'))
} else {
  dataset2_in_directory <- opt[['dataset2_in_directory']]
}
if (is.null(opt[['dataset2_in_prepend']])) {
  stop(paste('-p/--dataset2_in_prepend is an obligatory parameter\n'))
} else {
  dataset2_in_prepend <- opt[['dataset2_in_prepend']]
}
if (is.null(opt[['dataset2_in_append']])) {
  stop(paste('-s/--dataset2_in_append is an obligatory parameter\n'))
} else {
  dataset2_in <- opt[['dataset2_in_append']]
}

if (is.null(opt[['output_loc']])) {
  stop(paste('-o/--output_loc is an obligatory parameter\n'))
} else {
  output_loc <- opt[['output_loc']]
}

# check if the directory for the output exists, otherwise we would fail at the very last step
output_dir <- dirname(output_loc)
if (!(dir.exists(output_dir))) {
  stop(paste('directory for output table', output_loc, 'does not exist. This would mean the final step writing to file would fail, please create the directory first\n'))
}
# check if we are not accidentally overwriting the source file
if (dataset1_in == output_loc) {
  stop(paste('output file', output_loc, 'is the same as the dataset1 input file', dataset1_in, '. You should not overwrite your input! Please select a different output file.'))
}
if (dataset2_in == output_loc) {
  stop(paste('output file', output_loc, 'is the same as the dataset2 input file', dataset2_in, '. You should not overwrite your input! Please select a different output file.'))
}
# and if we are not looking at exactly the same file
if (dataset1_in == dataset2_in) {
  stop(paste('dataset1 file', dataset1_in, 'is the same as the dataset2 input file', dataset2_in, '. You need to supply two different files'))
}
# and finally, check if the names of the datasets are different
if (dataset1_name == dataset2_name) {
  stop(paste('names of the datasets that were supplied', dataset1_name, 'are the same, either make these different or use the default \'dataset1\' and \'dataset2\''))
}

# check if we are also writing the binary rds
if (!is.null(opt[['binary_rds_loc']])) {
  binary_rds_loc <- opt[['binary_rds_loc']]
  # also check if that folder exists
  rds_dir <- dirname(binary_rds_loc)
  if (!(dir.exists(rds_dir))) {
    stop(paste('directory for output rds', binary_rds_loc, 'does not exist. This would mean the final step writing to file would fail, please create the directory first\n'))
  }
  # and check for overwriting source files
  if (dataset1_in == binary_rds_loc) {
    stop(paste('rds output file', binary_rds_loc, 'is the same as the dataset1 input file', dataset1_in, '. You should not overwrite your input! Please select a different output file.'))
  }
  if (dataset2_in == binary_rds_loc) {
    stop(paste('rds output file', binary_rds_loc, 'is the same as the dataset2 input file', dataset2_in, '. You should not overwrite your input! Please select a different output file.'))
  }
}

# read the files
dataset1_fm_table <- read.table(dataset1_in, header = T, sep = '\t')
# make the regex for the dataset2 files
d2_files_regex <- paste0('^', dataset2_in_prepend, '(.+)', dataset2_in_append, '$')
# list all files in the directory
all_d2_files <- list.files(dataset2_in_directory, full.names = F)
# filter to those matching our regex
matching_d2_files <- all_d2_files[grepl(d2_files_regex, all_d2_files)]
# create a list to keep results per file
all_results_list <- list()
# loop over the files
for (d2_file in matching_d2_files) {
  # paste the full path
  d2_file_full <- file.path(dataset2_in_directory, d2_file)
  # read the file
  dataset2_fm_table <- read.table(dataset2_in, header = T, sep = '\t')
  # replace some of the column names
  colnames(dataset1_fm_table) <- gsub('lbf_cs_', 'CS', colnames(dataset1_fm_table))
  # get the result
  coloc_bfbf <- coloc_datasets(
    dataset1_fm_table, 
    dataset2_fm_table, 
    dataset1_name, 
    dataset2_name
  )
  # and put into the big list
  all_results_list[[d2_file]] <- coloc_bfbf
}

# write the rds if we can
if (!is.null(binary_rds_loc)) {
  saveRDS(all_results_list, binary_rds_loc)
  # with a checksum
  mdfiver::create_md5_for_file(binary_rds_loc)
}
# next extract all the tables
tbls <- list()
for (file_name in names(all_results_list)) {
    for (trait in names(all_results_list[[file_name]])) {
        tbls[[paste0(file_name, trait)]] <- all_results_list[[file_name]][[trait]][['table']]
        # add the filename to that table
        tbls[[paste0(file_name, trait)]][['filename']] <- file_name
    }
}
# merge all these
results_table_all <- do.call('rbind', tbls)

# make the output table
output_path_full <- output_loc
# gzip if that ends with
if (grepl('.gz$', output_loc)) {
  output_path_full <- gzfile(output_loc)
}
# write the result
write.table(results_table_all, output_path_full, row.names = F, col.names = T, sep = '\t', quote = F)
# and make a checksum
mdfiver::create_md5_for_file(output_loc)