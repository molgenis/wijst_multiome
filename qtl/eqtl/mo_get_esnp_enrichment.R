#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_esnp_enrichment.R
# Function: check if eSNPs are enriched to be in open chromatin regions
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(IRanges)


####################
# Functions        #
####################


get_qtls_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column='pval_nominal_threshold_global', nominal_significance_column='p_value', verbose=T) {
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
    if (!is.null(significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by significance column', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[significance_column]]) &
          cell_type_output[[significance_column]] < significance_cutoff, 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by significance column', nrow(cell_type_output)))
      }
    }
    if (!is.null(nominal_cutoff_column) & !is.null(nominal_significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by nominal cutoff', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
          !is.na(cell_type_output[[nominal_cutoff_column]]) & 
          !is.na(cell_type_output[[nominal_significance_column]]) &
          cell_type_output[[nominal_significance_column]] <= cell_type_output[[nominal_cutoff_column]], 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by nominal cutoff', nrow(cell_type_output)))
      }
    }
    # add the celltype as a column
    cell_type_output[['cell_type']] <- cell_type
    # add to the list
    qtls_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(qtls_per_celltype)
}


add_top_effect_annotation <- function(qtls_per_celltype, feature_column='feature_id', variant_column='snp_id', significance_column='p_value', decreasing=F) {
  # go through each of the cell types
  for (cell_type in names(qtls_per_celltype)) {
    # extract that table
    qtls_celltype <- qtls_per_celltype[[cell_type]]
    # order by the significancce
    qtls_celltype_ordered <- qtls_celltype[order(qtls_celltype[[significance_column]], decreasing = decreasing), ]
    # keep only the top effect
    qtls_celltype_top <- qtls_celltype_ordered[!duplicated(qtls_celltype_ordered[[feature_column]]), ]
    # and annotate in the original table if the variant-feature combination was the top one
    qtls_celltype[['is_top_variant']] <- paste(qtls_celltype[[variant_column]], qtls_celltype[[feature_column]]) %in% paste(qtls_celltype_top[[variant_column]], qtls_celltype_top[[feature_column]])
    # put that back in the list
    qtls_per_celltype[[cell_type]] <- qtls_celltype
  }
  return(qtls_per_celltype)
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
    qtl_regions_chromosome_variants <- unique(qtl_regions_chromosome[, c(variant_column_qtls, chromosome_column_qtls, position_column_qtls)])
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
# Main code        #
####################

# location of the eQTL output
qtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# get all the QTL output
qtl_output <- get_qtls_per_celltype_limix(qtl_output_loc)
# add the top effect information
qtl_output <- add_top_effect_annotation(qtl_output)
# merge all the results
qtl_output_all <- do.call('rbind', qtl_output)


# location of openness files
openness_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/'
# prepend and append
openness_prepend <- 'mo_peaks_lane1to80_'
openness_append <- '.bed'
# and the openness cell types
openness_cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

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

# add 'chr' to the chromosome
qtl_output_all[['snp_chromosome']] <- paste0('chr', qtl_output_all[['snp_chromosome']])

# we'll do this in chunks, because the data gets quite big
gene_chunk_size <- 250
# get all unique genes
qtl_genes <- unique(qtl_output_all[['feature_id']])
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
  qtl_data_chunk <- qtl_output_all[!is.na(qtl_output_all[['feature_id']]) & qtl_output_all[['feature_id']] %in% chunk_genes, ]
  
  # overlap based on cell type openness
  for (cell_type in names(openness_table_per_celltype)) {
    qtl_data_chunk <- qtl_merge_with_openness(
      qtl_data_chunk, 
      openness_table_per_celltype[[cell_type]], 
      overlapping_region_column = paste0('openness_region_', cell_type), 
      overlapping_openness_column = paste0('openness_openness_', cell_type), 
      variant_window_left = 0, 
      variant_window_right = 0, 
      variant_column_qtls = 'snp_id', 
      chromosome_column_qtls = 'snp_chromosome', 
      position_column_qtls = 'snp_position'
    )
  }
  
  # put chunk in the list
  chunk_results[[paste(as.character(chunk_start), as.character(chunk_end), sep = '-')]] <- qtl_data_chunk
  
  # update chunk
  chunk_start <- chunk_start + gene_chunk_size
}

# merge all the chunks back
qtl_output_all <- do.call('rbind', chunk_results)

# add the openness of any region if present
qtl_output_all[['openness_region_any']] <- apply(qtl_output_all, 1, function(x) {
  # get the openness columns
  any_openness <- x[grepl('openness_region_', names(x))]
  # and then anything that is not NA
  any_openness <- any_openness[!is.na(any_openness)]
  # if any values are not NA, we'll return that
  if (length(any_openness) > 0) {
    return(any_openness[1])
  }
  else {
    return(NA)
  }
})


# TESTING STUFF
# get CD4T top effects
cd4t_top <- qtl_output_all[
    qtl_output_all[['cell_type']] == 'CD4T' & 
    qtl_output_all[['is_top_variant']] == T, ]
# get the monocyte effects
monocyte_top <- qtl_output_all[
  qtl_output_all[['cell_type']] == 'monocyte' & 
    qtl_output_all[['is_top_variant']] == T, ]
# check if it was also the top in mono
cd4t_top[['top_mono']] <- paste(cd4t_top[['snp_id']], cd4t_top[['feature_id']]) %in% paste(monocyte_top[['snp_id']], monocyte_top[['feature_id']])
cd4t_top[['any_mono']] <- paste(cd4t_top[['snp_id']], cd4t_top[['feature_id']]) %in% paste(qtl_output_all[qtl_output_all[['cell_type']] == 'monocyte', ][['snp_id']], qtl_output_all[qtl_output_all[['cell_type']] == 'monocyte', ][['feature_id']])

# check if in open chromatin for CD4T
cd4t_top[['open_cd4t']] <- !is.na(cd4t_top[['openness_openness_CD4T']]) & cd4t_top[['openness_openness_CD4T']] >= 0.001
# check if in open chromatin for monocyte
cd4t_top[['open_monocyte']] <- !is.na(cd4t_top[['openness_openness_monocyte']]) & cd4t_top[['openness_openness_monocyte']] >= 0.001
