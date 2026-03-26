#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_overlapping_qtl_effects_per_cs.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)



####################
# Functions        #
####################

get_coloc_output <- function(coloc_output_loc, append='_eqtl_caqtl_coloc.tsv.gz') {
  # list all of the files in the directory
  coloc_files <- list.files(coloc_output_loc, recursive = F, full.names = F)
  # filter to files that fit the pattern
  coloc_files <- coloc_files[grepl(paste0('.*', append, '$'), coloc_files)]
  # now store each resulting file into a list
  coloc_per_ct <- list()
  # check each file
  for (coloc_file in coloc_files) {
    # read the file
    coloc_file_loc <- paste(coloc_output_loc, coloc_file, sep = '/')
    coloc_file_output <- fread(coloc_file_loc, header = T, sep = '\t')
    # grab the cell type from the filename
    ct <- gsub(append, '', coloc_file)
    # put that into the table
    coloc_file_output <- cbind(data.table('cell_type' = rep(ct, times = nrow(coloc_file_output))), coloc_file_output)
    # put in the list
    coloc_per_ct[[ct]] <- coloc_file_output
  }
  # merge all
  coloc_all <- do.call('rbind', coloc_per_ct)
  return(coloc_all)
}


get_coloc_variants <- function(finemapped_colocs, cell_type_column='ct', feature_column='Gene', variant_column='hit1', PP_H4_abf_column='PP.H4.abf', PP_H4_abf_column_cutoff=0.8, link_split_char='_') {
  # subset to significant
  finemapped_colocs_colocing <- finemapped_colocs[finemapped_colocs[[PP_H4_abf_column]] >= PP_H4_abf_column_cutoff, ]
  # extract the variant
  finemapped_colocs_colocing_short <- data.frame('variant' = finemapped_colocs_colocing[[variant_column]])
  # extract the feature link
  finemapped_colocs_colocing_short[, c('gene', 'region')] <- str_split_fixed(finemapped_colocs_colocing[[feature_column]], link_split_char, 2)
  # and the cell type
  finemapped_colocs_colocing_short[['cell_type']] <- finemapped_colocs_colocing[[cell_type_column]]
  return(finemapped_colocs_colocing_short)
}


#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @param gene_column which column to use as the gene identifier
#' @returns a list with the output tables per cell type
#' 
get_output_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all_qval_allchroms_fdr005_significant_cs.tsv.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, verbose=T, add_nominal_cutoff=T, nominal_p_value_column='p_value') {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  egenes_per_celltype <- list()
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
        print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[significance_column]]) &
          cell_type_output[[significance_column]] < significance_cutoff, 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
      }
      # add nominal cutoff by first getting the top effects
      cell_type_output_top <- cell_type_output[order(cell_type_output[[nominal_p_value_column]]), ]
      cell_type_output_top <- cell_type_output_top[!duplicated(cell_type_output_top[[gene_column]]), ]
      # get the highest still significant p value
      max_sig_p <- max(cell_type_output_top[[nominal_p_value_column]])
      # and add that to the output
      cell_type_output[['nominal_p_value_cutoff']] <- max_sig_p
    }
    # add to the list
    egenes_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(egenes_per_celltype)
}


add_credible_sets_output <- function(output_table_per_celltype, feature_column='feature_id', cs_column='CS', item_column='cell_type') {
  # filter on pval_nominal_threshold_global
  eqtl_outputs_filtered <- list()
  for (ct in names(output_table_per_celltype)) {
    # extract this eQTL output
    ct_eqtl_output <- output_table_per_celltype[[ct]]
    # # filter on the pval threshold
    # ct_eqtl_output <- ct_eqtl_output[ct_eqtl_output[['p_value']] < ct_eqtl_output[['pval_nominal_threshold_global']], ]
    # get the unique combinations of features and credible sets
    ct_eqtl_feature_cs <- unique(ct_eqtl_output[, c(feature_column, cs_column)])
    # get features with NA values
    ct_eqtl_feature_cs_na <- ct_eqtl_feature_cs[is.na(ct_eqtl_feature_cs[[cs_column]]), feature_column]
    # and the ones without
    ct_eqtl_feature_cs_nona <- ct_eqtl_feature_cs[!is.na(ct_eqtl_feature_cs[[cs_column]]), feature_column]
    # get the ones which are only with na
    ct_eqtl_feature_cs_onlyna <- setdiff(ct_eqtl_feature_cs_na, ct_eqtl_feature_cs_nona)
    if (length(ct_eqtl_feature_cs_onlyna)) {
      # for those, set the CS to simply be L1
      ct_eqtl_output[ct_eqtl_output[[feature_column]] %in% ct_eqtl_feature_cs_onlyna, ][[cs_column]] <- 'L1'
    }
    # get the ones which are partly na
    ct_eqtl_feature_cs_somena <- intersect(ct_eqtl_feature_cs_nona, ct_eqtl_feature_cs_na)
    if (length(ct_eqtl_feature_cs_somena) > 0) {
      ct_eqtl_output[(ct_eqtl_output[[feature_column]] %in% ct_eqtl_feature_cs_somena) & is.na(ct_eqtl_output[[cs_column]]), ][[cs_column]] <- 'L0'
    }
    # add ct
    ct_eqtl_output[[item_column]] <- ct
    # put back in list
    eqtl_outputs_filtered[[ct]] <- ct_eqtl_output
  }
  return(eqtl_outputs_filtered)
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

# location of the finemapped results
finemapped_colocs_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/colocalization/eqtl_caqtl/ut_and_24hca_significant/'
# read the tabel
finemapped_colocs <- get_coloc_output(finemapped_colocs_loc)
# add the way that we merged these
finemapped_colocs[['method']] <- 'coloc'
# filter the colocs
finemapped_colocs <- finemapped_colocs[!is.na(finemapped_colocs[['PP.H4.abf']]) & finemapped_colocs[['PP.H4.abf']] >= 0.75, ]

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
eqtl_outputs <- get_output_per_celltype_limix(eqtl_output_loc)
# location of the QTL outputs
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/finemapping/caqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
caqtl_outputs <- get_output_per_celltype_limix(caqtl_output_loc)

# add cs info
eqtl_outputs_cs <- add_credible_sets_output(eqtl_outputs)
caqtl_outputs_cs <- add_credible_sets_output(caqtl_outputs)

# merge all of them
eqtl_outputs_cs_all <- rbindlist(eqtl_outputs_cs, fill = T)
caqtl_outputs_cs_all <- rbindlist(caqtl_outputs_cs, fill = T)
# remove variants not at threshold
eqtl_outputs_cs_all <- eqtl_outputs_cs_all[eqtl_outputs_cs_all[['p_value']] <= eqtl_outputs_cs_all[['pval_nominal_threshold_global']], ]
caqtl_outputs_cs_all <- caqtl_outputs_cs_all[caqtl_outputs_cs_all[['p_value']] <= caqtl_outputs_cs_all[['pval_nominal_threshold_global']], ]
# get variants for each
eqtl_outputs_cs_all_variants <- eqtl_outputs_cs_all[['snp_id']]
caqtl_outputs_cs_all_variants <- caqtl_outputs_cs_all[['snp_id']]
# and cell types
eqtl_outputs_cs_all_cts <- eqtl_outputs_cs_all[['cell_type']]
caqtl_outputs_cs_all_cts <- caqtl_outputs_cs_all[['cell_type']]
# remove the columns
eqtl_outputs_cs_all[['snp_id']] <- NULL
caqtl_outputs_cs_all[['snp_id']] <- NULL
eqtl_outputs_cs_all[['cell_type']] <- NULL
caqtl_outputs_cs_all[['cell_type']] <- NULL
# rename so that the e and ca are added to the names
colnames(eqtl_outputs_cs_all) <- paste('e', colnames(eqtl_outputs_cs_all), sep = '_')
colnames(caqtl_outputs_cs_all) <- paste('ca', colnames(caqtl_outputs_cs_all), sep = '_')
# add variants back
eqtl_outputs_cs_all <- cbind(data.table('cell_type' = eqtl_outputs_cs_all_cts, 'variant_id' = eqtl_outputs_cs_all_variants), eqtl_outputs_cs_all)
caqtl_outputs_cs_all <- cbind(data.table('cell_type' = caqtl_outputs_cs_all_cts, 'variant_id' = caqtl_outputs_cs_all_variants), caqtl_outputs_cs_all)
# merge the outputs based on overlapping variant (for that cell type)
qtl_outputs_cs_overlapping <- merge(x = caqtl_outputs_cs_all, y = eqtl_outputs_cs_all, by = c('cell_type', 'variant_id'))

# add the CS info
finemapped_colocs[['e_CS']] <- eqtl_outputs_cs_all[
  match(
    paste(finemapped_colocs[['hit2']], finemapped_colocs[['trait2']]), 
    paste(eqtl_outputs_cs_all[['variant_id']], eqtl_outputs_cs_all[['e_feature_id']])
  ), 
][['e_CS']]

# get the unique CS-region-gene
finemapped_colocs_unique_overlap <- unique(finemapped_colocs[, c('trait1', 'trait2', 'e_CS')])
qtl_outputs_cs_unique_overlapping <- unique(qtl_outputs_cs_overlapping[, c('ca_feature_id', 'e_feature_id', 'e_CS')])
# set the same column names
colnames(finemapped_colocs_unique_overlap) <- c('region', 'gene', 'cs')
colnames(qtl_outputs_cs_unique_overlapping) <- c('region', 'gene', 'cs')
# rbind
colocs_unique_overlap_both <- unique(rbind(
  finemapped_colocs_unique_overlap, 
  qtl_outputs_cs_unique_overlapping
))
# for now, just keep the L0 as L1
colocs_unique_overlap_both[['cs']] <- gsub('L0', 'L1', colocs_unique_overlap_both[['cs']])
# do unique again
colocs_unique_overlap_both <- unique(colocs_unique_overlap_both)

# get unique overlaps
n_overlapping <- nrow(colocs_unique_overlap_both)

# get unique gene and CS combinations that are in the overlap
length(unique(paste(colocs_unique_overlap_both$gene, colocs_unique_overlap_both$cs)))
# 1753
# get the number of unique CS-gene pairs
eqtl_outputs_cs_all_nol0 <- eqtl_outputs_cs_all
eqtl_outputs_cs_all_nol0[['e_CS']] <- gsub('L0', 'L1', eqtl_outputs_cs_all_nol0[['e_CS']])
length(unique(paste(eqtl_outputs_cs_all_nol0$e_feature_id, eqtl_outputs_cs_all_nol0$e_CS)))
# 560
