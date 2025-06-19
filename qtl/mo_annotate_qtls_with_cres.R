#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_annotate_qtls_with_cres.R
# Function: annotate QTL outputs with DAR or CRE outputs
# Example: Rscript mo_annotate_qtls_with_cres.R \
# --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/monocyte/qtl_results_all_qval_allchroms_fdr005_significant.txt.gz \
# --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/monocyte/qtl_results_all_qval_allchroms_fdr005_significant_credar.txt.gz \
# --add_chrom chr
#
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(IRanges)
library(optparse)
library(stringr)
library(mdfiver)


####################
# Functions        #
####################


#' Merge QTL Data with Differentially Accessible Regions (DARs)
#'
#' This function identifies overlaps between QTL variants and DARs by chromosome and position,
#' and annotates the QTL data with the overlapping DAR region identifiers.
#'
#' @param qtl_input A data frame containing QTL data. Must include columns for variant ID, chromosome, and position.
#' @param dar_input A data frame containing DAR data. Must include a column with region strings in the format "chr:start-end".
#' @param variant_column_qtls Character string specifying the column name in `qtl_input` for variant IDs. Default is `'snp_id'`.
#' @param chromosome_column_qtls Character string specifying the column name in `qtl_input` for chromosome. Default is `'snp_chromosome'`.
#' @param position_column_qtls Character string specifying the column name in `qtl_input` for variant position. Default is `'snp_position'`.
#' @param region_column_dars Character string specifying the column name in `dar_input` that contains region strings. Default is `'region'`.
#'
#' @return A data frame containing the original QTL data with an additional column `overlapping_dar` indicating the DAR region each variant overlaps with (if any).
#'
#' @importFrom IRanges IRanges findOverlaps pintersect
#' @importFrom data.table data.table
#' @importFrom stringr str_split_fixed
#'
#' @examples
#' \dontrun{
#' merged_data <- qtl_merge_with_dars(qtl_df, dar_df)
#' }
#'
#' @export
qtl_merge_with_dars <- function(qtl_input, dar_input, variant_column_qtls='snp_id', chromosome_column_qtls='snp_chromosome', position_column_qtls='snp_position', region_column_dars='region', variant_window_left=0, variant_window_right=0) {
  # create a table of the DARs
  dar_table <- data.frame('chrom' = rep(NA, times = nrow(dar_input)), 'start' = rep(NA, times = nrow(dar_input)), 'end' = rep(NA, times = nrow(dar_input)))
  # split the dar output into regions
  dar_table[c('chrom', 'start', 'end')] <- str_split_fixed(dar_input[[region_column_dars]], '\\-|\\:', 3)
  # make start and end numeric
  dar_table[['start']] <- as.numeric(dar_table[['start']])
  dar_table[['end']] <- as.numeric(dar_table[['end']])
  # add the name of the region to that
  dar_table[[region_column_dars]] <- dar_input[[region_column_dars]]
  # make unique
  dar_table <- unique(dar_table)
  # get the chromosomes in the qtl data
  qtl_chroms <- unique(qtl_input[[chromosome_column_qtls]])
  # and in the dars
  dars_chroms <- unique(dar_table[['chrom']])
  # only do the ones present in both
  chroms_both <- intersect(qtl_chroms, dars_chroms)
  # we'll save the results in a list
  overlaps_per_chrom <- list()
  # and check each chromosome
  for (chrom in chroms_both) {
    # subset to this chrom
    dar_input_regions_chromosome <- dar_table[!is.na(dar_table[['chrom']]) & dar_table[['chrom']] == chrom, ]
    qtl_regions_chromosome <- qtl_input[!is.na(qtl_input[[chromosome_column_qtls]]) & qtl_input[[chromosome_column_qtls]] == chrom, ]
    # then subset to unique variants for the QTLs
    qtl_regions_chromosome_variants <- unique(qtl_regions_chromosome[, c(..variant_column_qtls, ..chromosome_column_qtls, ..position_column_qtls)])
    # turn into iranges objects
    dar_input_chromosome_iranges <- IRanges(start = dar_input_regions_chromosome[['start']], end = dar_input_regions_chromosome[['end']])
    # with the window supplied, making sure we don't get negative numbers
    qtl_region_chromosome_variants_windowed_starts <- qtl_regions_chromosome_variants[[position_column_qtls]] - variant_window_left
    qtl_region_chromosome_variants_windowed_starts[qtl_region_chromosome_variants_windowed_starts < 1] <- 1
    qtl_chromosome_iranges <- IRanges(start = qtl_region_chromosome_variants_windowed_starts, end = qtl_regions_chromosome_variants[[position_column_qtls]] + variant_window_right)
    # find overlaps
    feature_chromosome_overlaps <- findOverlaps(dar_input_chromosome_iranges, qtl_chromosome_iranges)
    # extract overlapping ranges
    overlapping_ranges <- pintersect(qtl_chromosome_iranges[subjectHits(feature_chromosome_overlaps)], dar_input_chromosome_iranges[queryHits(feature_chromosome_overlaps)])
    # create a  table for the overlaps
    overlaps_table <- data.table(
      'variant_id' = qtl_regions_chromosome_variants[[variant_column_qtls]][subjectHits(feature_chromosome_overlaps)],
      'overlapping_dar' = dar_input_regions_chromosome[[region_column_dars]][queryHits(feature_chromosome_overlaps)]
    )
    # add these positions
    qtl_regions_chromosome <- merge(qtl_regions_chromosome, overlaps_table, by.x = variant_column_qtls, by.y = 'variant_id', all.x = T, allow.cartesian=TRUE)
    # and put in the list
    overlaps_per_chrom[[as.character(chrom)]] <- qtl_regions_chromosome
  }
  # merge it all
  overlaps_all <- do.call('rbind', overlaps_per_chrom)
  return(overlaps_all)
}


#' Merge QTL Data with SCENIC CRE Annotations
#'
#' This function identifies overlaps between QTL variants and cis-regulatory elements (CREs) from SCENIC output,
#' and annotates the QTL data with overlapping CREs, transcription factors (TFs), and target genes.
#'
#' @param qtl_input A data frame containing QTL data. Must include columns for variant ID, chromosome, and position.
#' @param scenic_input A data frame containing SCENIC CRE annotations. Must include a region column in the format "chr:start-end", and columns for associated genes and transcription factors.
#' @param variant_column_qtls Character string specifying the column name in `qtl_input` for variant IDs. Default is `'snp_id'`.
#' @param chromosome_column_qtls Character string specifying the column name in `qtl_input` for chromosome. Default is `'snp_chromosome'`.
#' @param position_column_qtls Character string specifying the column name in `qtl_input` for variant position. Default is `'snp_position'`.
#' @param region_column_scenic Character string specifying the column name in `scenic_input` that contains region strings. Default is `'Region'`.
#' @param gene_column_scenic Character string specifying the column name in `scenic_input` for the target gene. Default is `'Gene'`.
#' @param tf_column_scenic Character string specifying the column name in `scenic_input` for the transcription factor. Default is `'TF'`.
#'
#' @return A data frame containing the original QTL data with additional columns:
#' \itemize{
#' \item `overlapping_cre`: the CRE region that overlaps with the variant (if any),
#' \item `cre_tf`: the transcription factor associated with the CRE,
#' \item `cre_gene`: the gene associated with the CRE.
#' }
#'
#' @importFrom IRanges IRanges findOverlaps pintersect
#' @importFrom data.table data.table
#' @importFrom stringr str_split_fixed
#'
#' @examples
#' \dontrun{
#' merged_data <- qtl_merge_with_scenic(qtl_df, scenic_df)
#' }
#'
#' @export
qtl_merge_with_scenic <- function(qtl_input, scenic_input, variant_column_qtls='snp_id', chromosome_column_qtls='snp_chromosome', position_column_qtls='snp_position', region_column_scenic='Region', gene_column_scenic='Gene', tf_column_scenic='TF', variant_window_left=0, variant_window_right=0) {
  # create a table of the scenics
  cre_table <- data.frame('chrom' = rep(NA, times = nrow(scenic_input)), 'start' = rep(NA, times = nrow(scenic_input)), 'end' = rep(NA, times = nrow(scenic_input)))
  # split the scenic output into regions
  cre_table[c('chrom', 'start', 'end')] <- str_split_fixed(scenic_input[[region_column_scenic]], '\\-|\\:', 3)
  # make start and end numeric
  cre_table[['start']] <- as.numeric(cre_table[['start']])
  cre_table[['end']] <- as.numeric(cre_table[['end']])
  # add the name of the region to that
  cre_table[[region_column_scenic]] <- scenic_input[[region_column_scenic]]
  # make that unique
  cre_table <- unique(cre_table)
  # add the gene and TF with a more specific name in the scenic data
  scenic_input[['cre_tf']] <- scenic_input[[tf_column_scenic]]
  scenic_input[['cre_gene']] <- scenic_input[[gene_column_scenic]]
  # get the chromosomes in the qtl data
  qtl_chroms <- unique(qtl_input[[chromosome_column_qtls]])
  # and in the dars
  scenic_chroms <- unique(cre_table[['chrom']])
  # only do the ones present in both
  chroms_both <- intersect(qtl_chroms, scenic_chroms)
  # we'll save the results in a list
  overlaps_per_chrom <- list()
  # and check each chromosome
  for (chrom in chroms_both) {
    # subset to this chrom
    cre_input_regions_chromosome <- cre_table[!is.na(cre_table[['chrom']]) & cre_table[['chrom']] == chrom, ]
    qtl_regions_chromosome <- qtl_input[!is.na(qtl_input[[chromosome_column_qtls]]) & qtl_input[[chromosome_column_qtls]] == chrom, ]
    # then subset to unique variants for the QTLs
    qtl_regions_chromosome_variants <- unique(qtl_regions_chromosome[, c(..variant_column_qtls, ..chromosome_column_qtls, ..position_column_qtls)])
    # turn into iranges objects
    cre_input_chromosome_iranges <- IRanges(start = cre_input_regions_chromosome[['start']], end = cre_input_regions_chromosome[['end']])
    # using the window, and making sure that we dont get negative numbers
    qtl_region_chromosome_variants_windowed_starts <- qtl_regions_chromosome_variants[[position_column_qtls]] - variant_window_left
    qtl_region_chromosome_variants_windowed_starts[qtl_region_chromosome_variants_windowed_starts < 1] <- 1
    qtl_chromosome_iranges <- IRanges(start = qtl_region_chromosome_variants_windowed_starts, end = qtl_regions_chromosome_variants[[position_column_qtls]] + variant_window_right)
    # find overlaps
    feature_chromosome_overlaps <- findOverlaps(cre_input_chromosome_iranges, qtl_chromosome_iranges)
    # extract overlapping ranges
    overlapping_ranges <- pintersect(qtl_chromosome_iranges[subjectHits(feature_chromosome_overlaps)], cre_input_chromosome_iranges[queryHits(feature_chromosome_overlaps)])
    # create a  table for the overlaps
    overlaps_table <- data.table(
      'feature_id' = qtl_regions_chromosome_variants[[variant_column_qtls]][subjectHits(feature_chromosome_overlaps)],
      'overlapping_cre' = cre_input_regions_chromosome[[region_column_scenic]][queryHits(feature_chromosome_overlaps)]
    )
    # add these positions
    qtl_regions_chromosome <- merge(qtl_regions_chromosome, overlaps_table, by.x = variant_column_qtls, by.y = 'feature_id', all.x = T, allow.cartesian=TRUE)
    # take the ones with a cre
    qtl_regions_chromosome_cre <- qtl_regions_chromosome[!is.na(qtl_regions_chromosome[['overlapping_cre']]), ]
    # now we need to add the eregulon for each of these cres
    overlaps_with_scenic <- merge(qtl_regions_chromosome_cre, scenic_input[, c(..region_column_scenic, 'cre_tf', 'cre_gene')], by.x = 'overlapping_cre', by.y = region_column_scenic, all.x = T, allow.cartesian=TRUE)
    # also take the ones without a cre
    qtl_regions_chromosome_creless <- qtl_regions_chromosome[is.na(qtl_regions_chromosome[['overlapping_cre']]), ]
    # and put those together
    overlaps_with_scenic <- rbind(overlaps_with_scenic, qtl_regions_chromosome_creless, fill = T)
    # order the columns so the CRE stuff is at the back
    column_order <- c(setdiff(colnames(overlaps_with_scenic), c('overlapping_cre', 'cre_tf', 'cre_gene')), 'overlapping_cre', 'cre_tf', 'cre_gene')
    overlaps_with_scenic <- overlaps_with_scenic[, ..column_order]
    # put in list
    overlaps_per_chrom[[as.character(chrom)]] <- overlaps_with_scenic
  }
  overlaps_all <- do.call('rbind', overlaps_per_chrom)
  return(overlaps_all)
}



#' Merge QTL Data with Openness Regions Based on Genomic Overlaps
#'
#' This function merges QTL (Quantitative Trait Loci) data with genomic openness data
#' by identifying overlaps between QTL variant positions and openness regions on the same chromosome.
#' It returns the QTL data annotated with overlapping openness region identifiers and openness values.
#'
#' @param qtl_input A `data.table` or `data.frame` containing QTL data.
#' @param openness_input A `data.table` or `data.frame` containing openness region data.
#' @param variant_column_qtls Name of the column in `qtl_input` that contains variant identifiers. Default is `'snp_id'`.
#' @param chromosome_column_qtls Name of the column in `qtl_input` that contains chromosome identifiers. Default is `'snp_chromosome'`.
#' @param position_column_qtls Name of the column in `qtl_input` that contains variant positions. Default is `'snp_position'`.
#' @param openness_column_region Name of the column in `openness_input` that contains openness region identifiers. Default is `'name'`.
#' @param opennes_column_openness Name of the column in `openness_input` that contains openness values. Default is `'pct_exp'`.
#' @param openness_column_chromosome Name of the column in `openness_input` that contains chromosome identifiers. Default is `'#chrom'`.
#' @param openness_column_start Name of the column in `openness_input` that contains region start positions. Default is `'start'`.
#' @param openness_column_end Name of the column in `openness_input` that contains region end positions. Default is `'end'`.
#' @param overlapping_region_column Name of the output column for overlapping openness region identifiers. Default is `'openness_region'`.
#' @param overlapping_openness_column Name of the output column for overlapping openness values. Default is `'openness_openness'`.
#'
#' @return A `data.table` containing the original QTL data with additional columns for overlapping openness region and openness value.
#'
#' @import data.table
#' @importFrom IRanges IRanges findOverlaps pintersect
#' @export
#'
#' @examples
#' # Example usage:
#' # result <- qtl_merge_with_openness(qtl_data, openness_data)
qtl_merge_with_openness <- function(qtl_input, openness_input, variant_column_qtls='snp_id', chromosome_column_qtls='snp_chromosome', position_column_qtls='snp_position', openness_column_region='name', opennes_column_openness='pct_exp', openness_column_chromosome='#chrom', openness_column_start='start', openness_column_end='end', overlapping_region_column='openness_region', overlapping_openness_column='openness_openness', variant_window_left=0, variant_window_right=0) {
  # get the chromosomes in the qtl data
  qtl_chroms <- unique(qtl_input[[chromosome_column_qtls]])
  # and in the dars
  openness_chroms <- unique(openness_input[[openness_column_chromosome]])
  # only do the ones present in both
  chroms_both <- intersect(qtl_chroms, openness_chroms)
  # we'll save the results in a list
  overlaps_per_chrom <- list()
  # and check each chromosome
  for (chrom in chroms_both) {
    # subset to this chrom
    openness_input_regions_chromosome <- openness_input[!is.na(openness_input[[openness_column_chromosome]]) & openness_input[[openness_column_chromosome]] == chrom, ]
    qtl_regions_chromosome <- qtl_input[!is.na(qtl_input[[chromosome_column_qtls]]) & qtl_input[[chromosome_column_qtls]] == chrom, ]
    # then subset to unique variants for the QTLs
    qtl_regions_chromosome_variants <- unique(qtl_regions_chromosome[, c(..variant_column_qtls, ..chromosome_column_qtls, ..position_column_qtls)])
    # turn into iranges objects
    openness_input_chromosome_iranges <- IRanges(start = openness_input_regions_chromosome[[openness_column_start]], end = openness_input_regions_chromosome[[openness_column_end]])
    # with the window supplied (making sure that we don't get a negative number)
    qtl_region_chromosome_variants_windowed_starts <- qtl_regions_chromosome_variants[[position_column_qtls]] - variant_window_left
    qtl_region_chromosome_variants_windowed_starts[qtl_region_chromosome_variants_windowed_starts < 1] <- 1
    qtl_chromosome_iranges <- IRanges(start = qtl_region_chromosome_variants_windowed_starts, end = qtl_regions_chromosome_variants[[position_column_qtls]] + variant_window_right)
    # find overlaps
    feature_chromosome_overlaps <- findOverlaps(openness_input_chromosome_iranges, qtl_chromosome_iranges)
    # extract overlapping ranges
    overlapping_ranges <- pintersect(qtl_chromosome_iranges[subjectHits(feature_chromosome_overlaps)], openness_input_chromosome_iranges[queryHits(feature_chromosome_overlaps)])
    # create a  table for the overlaps
    overlaps_table <- data.table(
      'variant_id' = qtl_regions_chromosome_variants[[variant_column_qtls]][subjectHits(feature_chromosome_overlaps)],
      'openness_region' = openness_input_regions_chromosome[[openness_column_region]][queryHits(feature_chromosome_overlaps)],
      'openness_openness' = openness_input_regions_chromosome[[opennes_column_openness]][queryHits(feature_chromosome_overlaps)]
    )
    # name them as we wanted
    colnames(overlaps_table) <- c('variant_id', overlapping_region_column, overlapping_openness_column)
    # add these positions
    qtl_regions_chromosome <- merge(qtl_regions_chromosome, overlaps_table, by.x = variant_column_qtls, by.y = 'variant_id', all.x = T, allow.cartesian=TRUE)
    # put in list
    overlaps_per_chrom[[as.character(chrom)]] <- qtl_regions_chromosome
  }
  overlaps_all <- do.call('rbind', overlaps_per_chrom)
  return(overlaps_all)
}

####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- F


####################
# Main code        #
####################

# location of the CREs
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'

# location of the DARs
dars_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/dar_detection/merged_major_and_minor_celltypes_120topics/wilcoxon/merged_major_and_minor_celltypes_120topics_dars.tsv.gz'

# location of openness files
openness_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/'
# prepend and append
openness_prepend <- 'mo_peaks_lane1to80_'
openness_append <- '.bed'
# and the openness cell types
openness_cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')


# make command line options
option_list <- list(
  make_option(c("-i", "--in"), type="character", default=NULL, 
              help="input QTL file to add information to", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output QTL file to save", metavar="character"), 
  make_option(c("-v", "--variant_column"), type="character", default='snp_id', 
              help="column denoting the variant [default: %default]", metavar="character"),
  make_option(c("-f", "--feature_column"), type="character", default='feature_id', 
              help="column denoting the feature [default: %default]", metavar="character"), 
  make_option(c("-c", "--chromosome_column"), type="character", default='snp_chromosome', 
              help="column chromosome of the variant [default: %default]", metavar="character"), 
  make_option(c("-p", "--position_column"), type="character", default='snp_position', 
              help="column denoting the position of the variant on the chromosome [default: %default]", metavar="character"), 
  make_option(c("-a", "--add_chrom"), type="character", default=NULL, 
              help="prepend to add to the chromsome column of the variant, something like 'chr' is common, leave parameter out for no prepend", metavar="character"), 
  make_option(c("-s", "--gene_chunk_size"), type="numeric", default=100, 
              help="the number of genes to process in a chunk [default: %default]", metavar="numeric"), 
  make_option(c("-w", "--variant_window_size"), type="numeric", default=0, 
              help="the window around the variant to use for overlaps (default is to only look at the variant position itself) [default: %default]", metavar="numeric"), 
  make_option(c("-r", "--remove_non_overlaps"), action="store_true", default=FALSE,
              help="remove QTL entries that show no overlap with DARs or CREs [default: %default]")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
qtl_in_loc <- NULL
qtl_out_loc <- NULL
variant_column <- NULL
feature_column <- NULL
chromosome_column <- NULL
position_column <- NULL
add_chrom <- NULL
gene_chunk_size <- NULL
remove_non_overlaps <- NULL
variant_window_size <- NULL

# load debug settings if set to debug mode
if (debug) {
  # qtl_in_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/monocyte/qtl_results_all_qval_allchroms_fdr005_significant.txt.gz'
  # qtl_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/monocyte/qtl_results_all_qval_allchroms_fdr005_significant_credar.txt.gz'
  qtl_in_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/eqtlgen_replications/cisonly_vs_transacting/qtl_tables/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping.txt'
  qtl_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/eqtlgen_replications/cisonly_vs_transacting/qtl_tables/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping_credar.txt.gz'
  # variant_column <- 'snp_id'
  variant_column <- 'variant'
  # feature_column <- 'feature_id'
  feature_column <- 'gene_name'
  # chromosome_column <- 'snp_chromosome'
  chromosome_column <- 'chromosome'
  # position_column <- 'snp_position'
  position_column <- 'bp'
  add_chrom <- 'chr'
  gene_chunk_size <- 100
  remove_non_overlaps <- T
  # variant_window_size <- 0
  variant_window_size <- 1000
  # let user know we are in debug mode
  warning('running in debug mode! parameters supplied will have no effect!')
} else {
  # there are some things we cannot allow
  if (is.null(opt[['in']])) {
    stop('input QTL file must be supplied')
  } else {
    qtl_in_loc <- opt[['in']]
  }
  if (is.null(opt[['out']])) {
    stop('output QTL file must be supplied')
  } else {
    qtl_out_loc <- opt[['out']]
  }
  # stop if we are overwriting our source file
  if (qtl_in_loc == qtl_out_loc) {
    stop('input and output are the same, do not overwrite your source file!')
  }
  # the others we can just fetch
  variant_column <- opt[['variant_column']]
  feature_column <- opt[['feature_column']]
  chromosome_column <- opt[['chromosome_column']]
  position_column <- opt[['position_column']]
  add_chrom <- opt[['add_chrom']]
  gene_chunk_size <- opt[['gene_chunk_size']]
  remove_non_overlaps <- opt[['remove_non_overlaps']]
  variant_window_size <- opt[['variant_window_size']]
}

# determine the left and right flank for the window around the variant
left_window <- floor(variant_window_size / 2)
right_window <- ceiling(variant_window_size / 2)

# read the qtl data
qtl_data <- fread(qtl_in_loc, header = T, sep = '\t')
# add prepend if supplied
if (!is.null(add_chrom)) {
  qtl_data[[chromosome_column]] <- paste0(add_chrom, qtl_data[[chromosome_column]])
}
# read the dars
message(paste('reading DAR data at', dars_output_loc))
dars <- fread(dars_output_loc, header = T, sep = '\t')
# read the scenic data
message(paste('reading scenic data at', scenic_output_loc))
scenic <- fread(scenic_output_loc, header = T, sep = '\t')

# we'll store openness data per cell type
openness_table_per_celltype <- list()
# read each of the openness tables
for (cell_type in openness_cell_types) {
  # we'll paste the path together
  cell_type_openness_loc <- paste0(openness_output_loc, '/', openness_prepend, cell_type, openness_append)
  # let the user know we are reading this data
  message(paste('reading openness file at', cell_type_openness_loc))
  # read the file
  openness_table_per_celltype[[cell_type]] <- fread(cell_type_openness_loc, header = T, sep = '\t')
}

# get all unique genes
qtl_genes <- unique(qtl_data[[feature_column]])
# save the result of each chunk
chunk_results <- list()
# go through the chunks
chunk_start <- 1
# message process
message(paste('processing', as.character(length(qtl_genes)), 'features in chunks of', as.character(gene_chunk_size)))
# keep checking each chunk
while(chunk_start < length(qtl_genes)) {
  # the end of the chunk
  chunk_end <- chunk_start + gene_chunk_size - 1
  # unless we don't have a full chunk left
  if (chunk_end > length(qtl_genes)) {
    chunk_end <- length(qtl_genes)
  }
  # message
  message(paste('processing chunk', chunk_start, 'to', chunk_end))
  # grab the genes of this chunk
  chunk_genes <- qtl_genes[chunk_start : chunk_end]
  # subset the data to those genes
  qtl_data_chunk <- qtl_data[!is.na(qtl_data[[feature_column]]) & qtl_data[[feature_column]] %in% chunk_genes, ]
  
  # overlap based on cell type openness
  for (cell_type in names(openness_table_per_celltype)) {
    qtl_data_chunk <- qtl_merge_with_openness(
      qtl_data_chunk, 
      openness_table_per_celltype[[cell_type]], 
      overlapping_region_column = paste0('openness_region_', cell_type), 
      overlapping_openness_column = paste0('openness_openness_', cell_type), 
      variant_window_left = left_window, 
      variant_window_right = right_window, 
      variant_column_qtls = variant_column, 
      chromosome_column_qtls = chromosome_column, 
      position_column_qtls = position_column
    )
  }
  
  # overlap based on DAR
  qtl_data_chunk <- qtl_merge_with_dars(
    qtl_data_chunk, dars, 
    variant_window_left = left_window, 
    variant_window_right = right_window, 
    variant_column_qtls = variant_column, 
    chromosome_column_qtls = chromosome_column, 
    position_column_qtls = position_column
  )
  
  # overlap based on scenic
  qtl_data_chunk <- qtl_merge_with_scenic(
    qtl_data_chunk, 
    scenic, 
    variant_window_left = left_window, 
    variant_window_right = right_window, variant_column_qtls = variant_column, 
    chromosome_column_qtls = chromosome_column, 
    position_column_qtls = position_column
  )
  
  # if requested, remove the entries that do not show any overlap (to conserve memory and disk)
  if (remove_non_overlaps) {
    # let the user know we are removing empty entries
    message('removing QTL entries without overlaps with DARs or CREs')
    # grab the cell type region column names
    cell_type_openness_region_columns <- colnames(qtl_data_chunk)[grep('openness_region', colnames(qtl_data_chunk))]
    # get the number of non-NA regions for the openness matched region names
    cell_type_openness_region_nonnacounts <- apply(qtl_data_chunk[, ..cell_type_openness_region_columns], 1, function(x){return(sum(!is.na(x)))})
    # remove entries that have no CRE, no DAR, and no openness matched regions
    qtl_data_chunk <- qtl_data_chunk[
      !is.na(qtl_data_chunk[['overlapping_cre']]) | !is.na(qtl_data_chunk[['overlapping_dar']]) | cell_type_openness_region_nonnacounts > 0, 
    ]
  }
  
  # put chunk in the list
  chunk_results[[paste(as.character(chunk_start), as.character(chunk_end), sep = '-')]] <- qtl_data_chunk
  
  # update chunk
  chunk_start <- chunk_start + gene_chunk_size
}

# message again
message(paste('merging', as.character(length(names(chunk_results))), 'chunks'))
# merge the chunks
qtl_data <- do.call('rbind', chunk_results)

# make the output location filehandle
output_loc_fh <- qtl_out_loc
# gz filehandle, if the output location ends with .gz
if (grepl('.gz$', output_loc_fh)) {
  # gzip if ends with .gz
  output_loc_fh <- gzfile(output_loc_fh)
}
# final step
message(paste('writing output', qtl_out_loc))
# write the resulting file
write.table(qtl_data, output_loc_fh, row.names = F, col.names = T, sep = '\t', quote = F)
# make a checksum as well
mdfiver::create_md5_for_file(qtl_out_loc)

# and let them know we are done
message('finished')
