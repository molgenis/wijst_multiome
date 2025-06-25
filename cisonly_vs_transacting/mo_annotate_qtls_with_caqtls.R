#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_annotate_qtls_with_caqtls.R
# Function: annotate QTL outputs with caQTLs from multiome
# Example: Rscript ./mo_annotate_qtls_with_caqtls.R \
# --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/eqtlgen_replications/cisonly_vs_transacting/qtl_tables/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping.txt.gz \
# --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/eqtlgen_replications/cisonly_vs_transacting/qtl_tables/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping_caqtls.txt.gz \
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




merge_qtl_outputs_horizontally <- function(qtl_output_folder, cell_types=NULL, filename='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz', feature_column='feature_id', add_cutoff=T, significance_column_nominal='p_value', significance_column_emperical='feature_q_value', significance_cutoff_emperical=0.05, columns_to_match=c('feature_id', 'snp_id', 'snp_chromosome', 'snp_position'), columns_to_keep=c('assessed_allele', 'p_value', 'beta', 'beta_se', 'empirical_feature_p_value', 'feature_q_value')) {
  # list the folders in the directory
  cell_type_folders <- list.dirs(qtl_output_folder, recursive = F, full.names = F)
  # filter on specific cell types if requested
  if (!is.null(cell_types)) {
    cell_type_folders <- intersect(cell_type_folders, cell_types)
  }
  # initialize the variable to store everything
  horizontal_qtl_table <- NULL
  # check each cell type
  for (cell_type in cell_type_folders) {
    # read the file
    qtl_file_loc <- paste(qtl_output_folder, cell_type, filename, sep = '/')
    qtl_file <- fread(qtl_file_loc, header = T, sep = '\t')
    
    # filter on minimal significance
    if (add_cutoff) {
      # sort by nominal significance
      qtl_file_top <- qtl_file[order(qtl_file[[significance_column_nominal]]), ]
      # select the top effects
      qtl_file_top <- qtl_file_top[!duplicated(qtl_file_top[[feature_column]]), ]
      # filter on significance
      qtl_file_top <- qtl_file_top[qtl_file_top[[significance_column_emperical]] < significance_cutoff_emperical, ]
      # and get the smallest significant effect
      max_nominal_p <- max(qtl_file_top[[significance_column_nominal]])
      # filter on this
      qtl_file <- qtl_file[qtl_file[[significance_column_nominal]] <= max_nominal_p, ]
    }
    # rename the fixed columns
    qtl_file_matching <- qtl_file[, c(..columns_to_match)]
    qtl_file_keep <- qtl_file[, c(..columns_to_keep)]
    colnames(qtl_file_keep) <- paste(colnames(qtl_file_keep), cell_type, sep = '_')
    # and merge back together
    qtl_file <- cbind(qtl_file_matching, qtl_file_keep)
    # check if we have a table already
    if (is.null(horizontal_qtl_table)) {
      horizontal_qtl_table <- qtl_file
    }
    # otherwise we need to merge
    else {
      horizontal_qtl_table <- merge(horizontal_qtl_table, qtl_file, all = T, by = columns_to_match)
    }
  }
  return(horizontal_qtl_table)
}


merge_qtls_to_caqtls <- function(qtl_input, caqtl_input, variant_column_qtls='snp_id', chromosome_column_qtls='snp_chromosome', position_column_qtls='snp_position', caqtl_column_variant='snp_id', caqtl_column_chromosome='snp_chromosome', caqtl_column_position='snp_position', variant_window_left=0, variant_window_right=0) {
  # get the chromosomes in the qtl data
  qtl_chroms <- unique(qtl_input[[chromosome_column_qtls]])
  # and in the dars
  caqtl_chroms <- unique(caqtl_input[[caqtl_column_chromosome]])
  # only do the ones present in both
  chroms_both <- intersect(qtl_chroms, caqtl_chroms)
  # we'll save the results in a list
  overlaps_per_chrom <- list()
  # and check each chromosome
  for (chrom in chroms_both) {
    # subset to this chrom
    caqtl_input_regions_chromosome <- caqtl_input[!is.na(caqtl_input[[caqtl_column_chromosome]]) & caqtl_input[[caqtl_column_chromosome]] == chrom, ]
    qtl_regions_chromosome <- qtl_input[!is.na(qtl_input[[chromosome_column_qtls]]) & qtl_input[[chromosome_column_qtls]] == chrom, ]
    # then subset to unique variants for the QTLs
    qtl_regions_chromosome_variants <- unique(qtl_regions_chromosome[, c(..variant_column_qtls, ..chromosome_column_qtls, ..position_column_qtls)])
    # turn into iranges objects
    caqtl_input_chromosome_iranges <- IRanges(start = caqtl_input_regions_chromosome[[caqtl_column_position]], end = caqtl_input_regions_chromosome[[caqtl_column_position]])
    # with the window supplied (making sure that we don't get a negative number)
    qtl_region_chromosome_variants_windowed_starts <- qtl_regions_chromosome_variants[[position_column_qtls]] - variant_window_left
    qtl_region_chromosome_variants_windowed_starts[qtl_region_chromosome_variants_windowed_starts < 1] <- 1
    qtl_chromosome_iranges <- IRanges(start = qtl_region_chromosome_variants_windowed_starts, end = qtl_regions_chromosome_variants[[position_column_qtls]] + variant_window_right)
    # find overlaps
    feature_chromosome_overlaps <- findOverlaps(caqtl_input_chromosome_iranges, qtl_chromosome_iranges)
    # extract overlapping ranges
    overlapping_ranges <- pintersect(qtl_chromosome_iranges[subjectHits(feature_chromosome_overlaps)], caqtl_input_chromosome_iranges[queryHits(feature_chromosome_overlaps)])
    # create a  table for the overlaps
    overlaps_table <- data.table(
      'variant_id' = qtl_regions_chromosome_variants[[variant_column_qtls]][subjectHits(feature_chromosome_overlaps)],
      'variant_id_caqtl' = caqtl_input_regions_chromosome[[caqtl_column_variant]][queryHits(feature_chromosome_overlaps)]
    )
    # add these positions
    qtl_regions_chromosome <- merge(qtl_regions_chromosome, overlaps_table, by.x = variant_column_qtls, by.y = 'variant_id', all.x = T, allow.cartesian=TRUE)
    # rename the columns in the qtl output
    colnames(caqtl_input_regions_chromosome) <- paste('caqtl', colnames(caqtl_input_regions_chromosome), sep = '_')
    # now also add the qtl data
    qtl_regions_chromosome <- merge(qtl_regions_chromosome, caqtl_input_regions_chromosome, by.x = 'variant_id_caqtl', by.y = paste('caqtl', caqtl_column_variant, sep = '_'), all.x = T, allow.cartesian = T)
    # add the variant ID from the caQTLs as the last column
    reordered_column_names <- colnames(qtl_regions_chromosome)[c(2:ncol(qtl_regions_chromosome), 1)]
    qtl_regions_chromosome <- qtl_regions_chromosome[, ..reordered_column_names]
    
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
debug <- T


####################
# Main code        #
####################

# location of the caQTL output
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'
# and the openness cell types
caqtl_cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

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
  qtl_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/eqtlgen_replications/cisonly_vs_transacting/qtl_tables/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping_caqtls.txt.gz'
  # variant_column <- 'snp_id'
  variant_column <- 'variant'
  # feature_column <- 'feature_id'
  feature_column <- 'gene_name'
  # chromosome_column <- 'snp_chromosome'
  chromosome_column <- 'chromosome'
  # position_column <- 'snp_position'
  position_column <- 'bp'
  add_chrom <- ''
  gene_chunk_size <- 100
  remove_non_overlaps <- T
  variant_window_size <- 0
  # variant_window_size <- 1000
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

# get the caqtl data
caqtl_data <- merge_qtl_outputs_horizontally(caqtl_output_loc)

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
  # add qtl data
  qtl_data_chunk <- merge_qtls_to_caqtls(qtl_data_chunk, 
                                   caqtl_input = caqtl_data, 
                                   variant_column_qtls = variant_column, 
                                   chromosome_column_qtls = chromosome_column, 
                                   position_column_qtls = position_column, 
                                   caqtl_column_variant='snp_id', 
                                   caqtl_column_chromosome='snp_chromosome', 
                                   caqtl_column_position='snp_position', 
                                   variant_window_left=0, 
                                   variant_window_right=0)
  # put chunk in the list
  chunk_results[[paste(as.character(chunk_start), as.character(chunk_end), sep = '-')]] <- qtl_data_chunk
  
  # update chunk
  chunk_start <- chunk_start + gene_chunk_size
}
# message again
message(paste('merging', as.character(length(names(chunk_results))), 'chunks'))
# merge the chunks
qtl_data <- do.call('rbind', chunk_results)
qtl_data <- unique(qtl_data)

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
