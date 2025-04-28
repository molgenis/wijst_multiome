#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: funmeta_lookup_funmeta_variants.R
# Function: look up the variants from the fungal meta-analysis in the in-house datasets
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(qvalue)
library(GenomicRanges)
library(stringr)


####################
# Functions        #
####################

#' Add Q-Values to QTL Table
#'
#' This function adds q-values to a QTL table based on the specified p-value column. It ensures that only the first entry for each feature is considered.
#'
#' @param qtl_table Data table. The QTL table containing p-values and feature IDs.
#' @param p_value_column Character. The name of the column containing p-values. Default is `'empirical_feature_p_value'`.
#' @param feature_id_column Character. The name of the column containing feature IDs. Default is `'feature_id'`.
#'
#' @return Data table. A new table with feature IDs and their corresponding q-values.
#'
#' @examples
#' qtl_table <- data.table(feature_id = c('gene1', 'gene2', 'gene1', 'gene3'),
#'                         empirical_feature_p_value = c(0.01, 0.02, 0.03, 0.04))
#' add_qvalue(qtl_table)
add_qvalue <- function(qtl_table, p_value_column='empirical_feature_p_value', feature_id_column='feature_id') {
  # sort by significance
  qtl_table <- qtl_table[order(qtl_table[[p_value_column]]), ]
  # now only get the first entry
  qtl_table_unique <- qtl_table[!duplicated(qtl_table[[feature_id_column]]), ]
  # extract the features
  features <- qtl_table_unique[[feature_id_column]]
  # and p values
  p_values <- qtl_table_unique[[p_value_column]]
  # get q values for those
  q_values <- qvalue(p_values)$qvalues
  # make into new table
  q_value_table <- data.table('features' = features, 'qval' = q_values)
  colnames(q_value_table) <- c(feature_id_column, 'qval')
  return(q_value_table)
}


#' Get QTLs from Files
#'
#' This function reads QTL files from a specified location, filters them based on given variants, and optionally adds q-values. It returns a combined data table of QTLs for all cell types.
#'
#' @param qtls_loc Character. The location of the QTL files.
#' @param variants Character vector. The variants to filter the QTLs by.
#' @param qtl_prepend Character. The prefix for QTL file names. Default is `''`.
#' @param qtl_append Character. The suffix for QTL file names. Default is `'.tsv.gz'`.
#' @param variant_column Character. The name of the column containing variant IDs. Default is `'snp_id'`.
#' @param add_qvalue Logical. Whether to add q-values to the QTLs. Default is `TRUE`.
#' @param p_value_column Character. The name of the column containing p-values. Default is `'empirical_feature_p_value'`.
#' @param feature_id_column Character. The name of the column containing feature IDs. Default is `'feature_id'`.
#' @param alpha_param_clean Logical. Whether to clean entries based on alpha parameters. Default is `TRUE`.
#' @param alpha_param_column Character. The name of the column containing alpha parameters. Default is `'alpha_param'`.
#' @param use_snp_pos_as_name Logical. Whether to use SNP positions as variant names. Default is `FALSE`.
#' @param snp_chrom_column Character. The name of the column containing SNP chromosome information. Default is `'snp_chromosome'`.
#' @param snp_pos_column Character. The name of the column containing SNP position information. Default is `'snp_position'`.
#'
#' @return Data table. A combined table of QTLs for all cell types, with optional q-values added.
#'
#' @examples
#' qtls_loc <- "path/to/qtls"
#' variants <- c("rs123", "rs456")
#' get_qtls_from_files(qtls_loc, variants)
get_qtls_from_files <- function(qtls_loc, variants, qtl_prepend='', qtl_append='.tsv.gz', variant_column='snp_id', add_qvalue=T, p_value_column='empirical_feature_p_value', feature_id_column='feature_id', alpha_param_clean=T, alpha_param_column='alpha_param', use_snp_pos_as_name=F, snp_chrom_column='snp_chromosome', snp_pos_column='snp_position') {
  # we'll store a results per cell type
  qtls_celltype <- list()
  # list all the files
  qtl_files <- list.files(qtls_loc, full.names = F, recursive = F, include.dirs = F)
  # paste expression together
  file_pattern <- paste0(qtl_prepend, '.*', qtl_append)
  # and filter those files
  qtl_files <- qtl_files[grep(pattern = file_pattern, x = qtl_files)]
  # read each file
  for (qtl_file in qtl_files) {
    # get the full path
    qtl_file_full <- paste0(qtls_loc, '/', qtl_file)
    # read the file
    qtls <- fread(qtl_file_full, header = T, sep = '\t')
    # filter to only include the variants we care about
    if (add_qvalue) {
      # remove entries with weird alpha params
      if(alpha_param_clean) {
        qtls <- qtls[!(qtls[[alpha_param_column]] < .2 | qtls[[alpha_param_column]] > 5), ]
      }
      # get qvalues
      qvalue_table <- add_qvalue(qtl_table = qtls, p_value_column = p_value_column, feature_id_column = feature_id_column)
      # add the tables together
      qtls <- merge(qtls, qvalue_table, by = feature_id_column, all = T)
    }
    # rename if requested
    if (use_snp_pos_as_name) {
      qtls[[variant_column]] <- paste(qtls[[snp_chrom_column]], qtls[[snp_pos_column]], sep = ':')
    }
    # filter on only the variants we care about
    qtls <- qtls[qtls[[variant_column]] %in% variants, ]
    # get the cell type from the name
    cell_type <- gsub(qtl_prepend, '', qtl_file)
    cell_type <- gsub(qtl_append, '', cell_type)
    # add cell type as column
    qtls <- cbind(data.table('cell_type' = rep(cell_type, times = nrow(qtls))), qtls)
    # put in the list
    qtls_celltype[[cell_type]] <- qtls
  }
  # merge all together
  qtls_all <- do.call('rbind', qtls_celltype)
  return(qtls_all)
}


#' Get QTLs from Foldered Files
#'
#' This function reads QTL files from specified folders, filters them based on given variants and significance filters, and returns a combined data table of QTLs for all cell types.
#'
#' @param qtls_loc Character. The location of the QTL folders.
#' @param variants Character vector. The variants to filter the QTLs by.
#' @param qtl_file_prepend Character. The prefix for QTL file names. Default is `'qtl_results_all_qval_'`.
#' @param qtl_file_append Character. The suffix for QTL file names. Default is `'.txt.gz'`.
#' @param variant_column Character. The name of the column containing variant IDs. Default is `'snp_id'`.
#' @param significance_filters List. A list of significance filters to apply. Default is `list('p_value' = 0.05)`.
#' @param use_snp_pos_as_name Logical. Whether to use SNP positions as variant names. Default is `FALSE`.
#' @param snp_chrom_column Character. The name of the column containing SNP chromosome information. Default is `'snp_chromosome'`.
#' @param snp_pos_column Character. The name of the column containing SNP position information. Default is `'snp_position'`.
#' @param folder_include Character vector. A vector of folder names to include. Default is `NULL`.
#'
#' @return Data table. A combined table of QTLs for all cell types, with optional significance filtering applied.
#'
#' @examples
#' qtls_loc <- "path/to/qtls"
#' variants <- c("rs123", "rs456")
#' get_qtls_foldered(qtls_loc, variants)
get_qtls_foldered <- function(qtls_loc, variants, qtl_file_prepend='qtl_results_all_qval_', qtl_file_append='.txt.gz', variant_column='snp_id', significance_filters=list('p_value' = 0.05), use_snp_pos_as_name=F, snp_chrom_column='snp_chromosome', snp_pos_column='snp_position', folder_include=NULL) {
  # list the folders in this directory
  qtl_folders <- list.dirs(qtls_loc, recursive = F, full.names = F)
  # if chosen, filter them
  if (!is.null(folder_include)) {
    qtl_folders <- intersect(qtl_folders, folder_include)
  }
  # create a list to save all the results
  results_per_celltype <- list()
  # check each folder
  for (qtl_folder in qtl_folders) {
    # paste together the folder we are looking at
    qtl_folder_celltype <- paste(qtls_loc, qtl_folder, sep = '/')
    # list all the files
    qtl_celltype_files <- list.files(qtl_folder_celltype, full.names = F)
    # paste together a pattern
    qtl_celltype_pattern <- paste0('^', qtl_file_prepend, '\\d+', qtl_file_append, '$')
    # filter by that cell type pattern
    qtl_celltype_files <- qtl_celltype_files[grep(qtl_celltype_pattern, qtl_celltype_files)]
    # check each file
    for (qtl_celltype_file in qtl_celltype_files) {
      # paste together the path to the file
      qtl_celltype_file <- paste(qtl_folder_celltype, qtl_celltype_file, sep = '/')
      # read the file
      qtl_celltype_file_contents <- fread(qtl_celltype_file, header = T, sep = '\t')
      # rename if requested
      if (use_snp_pos_as_name) {
        qtl_celltype_file_contents[[variant_column]] <- paste(qtl_celltype_file_contents[[snp_chrom_column]], qtl_celltype_file_contents[[snp_pos_column]], sep = ':')
      }
      # filter if requested
      if (!is.null(significance_filters)) {
        # check each of the significance filters
        for (significance_column in names(significance_filters)) {
          # now filter using that column (list key) and value (list value)
          qtl_celltype_file_contents <- qtl_celltype_file_contents[qtl_celltype_file_contents[[significance_column]] < significance_filters[[significance_column]], ]
        }
      }
      # filter on only the variants we care about
      qtl_celltype_file_contents <- qtl_celltype_file_contents[qtl_celltype_file_contents[[variant_column]] %in% variants, ]
      # add the celltype as a column
      qtl_celltype_file_contents <- cbind(data.table('cell_type' = rep(qtl_folder, times = nrow(qtl_celltype_file_contents))), qtl_celltype_file_contents)
      # add to the list
      results_per_celltype[[paste(qtl_folder, qtl_celltype_file, sep = '/')]] <- qtl_celltype_file_contents
    }
  }
  # merge all together
  results_all <- do.call('rbind', results_per_celltype)
  return(results_all)
}


#' Get QTLs from Foldered Files (Full)
#'
#' This function reads QTL files from specified folders, filters them based on given variants and significance filters, optionally adds q-values, and returns a combined data table of QTLs for all cell types.
#'
#' @param qtls_loc Character. The location of the QTL folders.
#' @param variants Character vector. The variants to filter the QTLs by.
#' @param qtl_file Character. The name of the QTL file. Default is `'qtl_results_all.txt.gz'`.
#' @param variant_column Character. The name of the column containing variant IDs. Default is `'snp_id'`.
#' @param add_qvalue Logical. Whether to add q-values to the QTLs. Default is `TRUE`.
#' @param p_value_column Character. The name of the column containing p-values. Default is `'empirical_feature_p_value'`.
#' @param feature_id_column Character. The name of the column containing feature IDs. Default is `'feature_id'`.
#' @param alpha_param_clean Logical. Whether to clean entries based on alpha parameters. Default is `TRUE`.
#' @param alpha_param_column Character. The name of the column containing alpha parameters. Default is `'alpha_param'`.
#' @param significance_filters List. A list of significance filters to apply. Default is `list('p_value' = 0.05)`.
#' @param use_snp_pos_as_name Logical. Whether to use SNP positions as variant names. Default is `FALSE`.
#' @param snp_chrom_column Character. The name of the column containing SNP chromosome information. Default is `'snp_chromosome'`.
#' @param snp_pos_column Character. The name of the column containing SNP position information. Default is `'snp_position'`.
#' @param folder_include Character vector. A vector of folder names to include. Default is `NULL`.
#' @param sep Character. The field separator character. Default is `'\t'`.
#' @param verbose Logical. Whether to print messages during processing. Default is `FALSE`.
#'
#' @return Data table. A combined table of QTLs for all cell types, with optional q-values and significance filtering applied.
#'
#' @examples
#' qtls_loc <- "path/to/qtls"
#' variants <- c("rs123", "rs456")
#' get_qtls_foldered_full(qtls_loc, variants)
get_qtls_foldered_full <- function(qtls_loc, variants, qtl_file='qtl_results_all.txt.gz', variant_column='snp_id', add_qvalue=T, p_value_column='empirical_feature_p_value', feature_id_column='feature_id', alpha_param_clean=T, alpha_param_column='alpha_param', significance_filters=list('p_value' = 0.05), use_snp_pos_as_name=F, snp_chrom_column='snp_chromosome', snp_pos_column='snp_position', folder_include=NULL, sep = '\t', verbose=F) {
  # list the folders in this directory
  qtl_folders <- list.dirs(qtls_loc, recursive = F, full.names = F)
  # if chosen, filter them
  if (!is.null(folder_include)) {
    qtl_folders <- intersect(qtl_folders, folder_include)
  }
  # create a list to save all the results
  results_per_celltype <- list()
  # check each folder
  for (qtl_folder in qtl_folders) {
    # paste together the path to the file
    qtl_celltype_file <- paste(qtls_loc, qtl_folder, qtl_file, sep = '/')
    if (verbose) {
      message(paste('reading', qtl_celltype_file))
    }
    # read the file
    qtl_celltype_file_contents <- fread(qtl_celltype_file, header = T, sep = sep)
    # filter to only include the variants we care about
    if (add_qvalue) {
      # remove entries with weird alpha params
      if(alpha_param_clean) {
        if (verbose) {
          message(paste('cleaning alpha param for', qtl_celltype_file))
        }
        qtl_celltype_file_contents <- qtl_celltype_file_contents[!(qtl_celltype_file_contents[[alpha_param_column]] < .2 | qtl_celltype_file_contents[[alpha_param_column]] > 5), ]
      }
      if (verbose) {
        message(paste('adding qvalue for', qtl_celltype_file))
      }
      # get qvalues
      qvalue_table <- add_qvalue(qtl_table = qtl_celltype_file_contents, p_value_column = p_value_column, feature_id_column = feature_id_column)
      # add the tables together
      qtl_celltype_file_contents <- merge(qtl_celltype_file_contents, qvalue_table, by = feature_id_column, all = T)
    }
    # rename if requested
    if (use_snp_pos_as_name) {
      qtl_celltype_file_contents[[variant_column]] <- paste(qtl_celltype_file_contents[[snp_chrom_column]], qtl_celltype_file_contents[[snp_pos_column]], sep = ':')
    }
    # filter if requested
    if (!is.null(significance_filters)) {
      # check each of the significance filters
      for (significance_column in names(significance_filters)) {
        # now filter using that column (list key) and value (list value)
        qtl_celltype_file_contents <- qtl_celltype_file_contents[qtl_celltype_file_contents[[significance_column]] < significance_filters[[significance_column]], ]
      }
    }
    # filter on only the variants we care about
    qtl_celltype_file_contents <- qtl_celltype_file_contents[qtl_celltype_file_contents[[variant_column]] %in% variants, ]
    # add the celltype as a column
    qtl_celltype_file_contents <- cbind(data.table('cell_type' = rep(qtl_folder, times = nrow(qtl_celltype_file_contents))), qtl_celltype_file_contents)
    # add to the list
    results_per_celltype[[paste(qtl_folder, qtl_celltype_file, sep = '/')]] <- qtl_celltype_file_contents
  }
  # merge all together
  results_all <- do.call('rbind', results_per_celltype)
  return(results_all)
}


check_overlaps_beds <- function(beds_loc, variants, bed_prepend='mo_peaks_lane1to80_wrna_', bed_append='.bed', conditions=c('UT', '24hCA'), variant_chromosome_field='chromosome', variant_start_field='start', variant_stop_field='stop') {
  # store the results in a list
  overlaps_per_bed <- list()
  # convert the variants to genomic ranges
  variants_ranges <- GenomicRanges::makeGRangesFromDataFrame(variants_all, 
                                          ignore.strand = T, 
                                          seqnames.field = c(variant_chromosome_field), 
                                          start.field = variant_start_field, 
                                          end.field = variant_stop_field, 
                                          seqinfo = NULL)
  # check each condition
  for (condition in conditions) {
    # paste together the path
    ranges_condition_pattern <- paste('^', bed_prepend, condition, '_.*', bed_append, sep = '')
    # get each file
    files_condition <- list.files(beds_loc, pattern = ranges_condition_pattern, full.names = F, recursive = F, include.dirs = F)
    # now check each file
    for (cell_type_file in files_condition) {
      # read the file
      cell_type_contents <- read.table(paste(beds_loc, cell_type_file, sep = '/'), header = T, sep = '\t', comment.char = '', check.names = F)
      # turn it into a genomic ranges object
      cell_type_ranges <- GenomicRanges::makeGRangesFromDataFrame(df = cell_type_contents, ignore.strand = T, seqnames.field = c('#chrom'), start.field = 'start', end.field = 'end', seqinfo = NULL)
      # get the overlap between the two
      cell_type_overlaps <- GenomicRanges::findOverlaps(variants_ranges, cell_type_ranges, ignore.strand = T)
      cell_type_overlaps <- data.frame(cell_type_overlaps)
      # extract the cell type
      cell_type <- gsub(paste0('^', bed_prepend, condition, '_'), '', cell_type_file)
      cell_type <- gsub(paste0(bed_append, '$'), '', cell_type)
      # combine the results
      all_overlap <- do.call('cbind', list(
        variants[cell_type_overlaps[['queryHits']], ],
        data.table('cell_type' = rep(cell_type, times = nrow(cell_type_overlaps)), 'condition' = rep(condition, times = nrow(cell_type_overlaps))), 
        cell_type_contents[cell_type_overlaps[['subjectHits']], ]
      ))
      # put in the list
      overlaps_per_bed[[paste(condition, cell_type, sep = '_')]] <- all_overlap
    }
  }
  # combine everything
  all_combined <- do.call('rbind', overlaps_per_bed)
  return(all_combined)
}


####################
# Main Code        #
####################

# we will look up the variants, and their proxy variants
variant1_ld_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/proxy_variants/proxy_rs3807184.txt'
variant2_ld_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/proxy_variants/proxy_rs413524.txt'
variant3_ld_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/proxy_variants/proxy_rs113472423.txt'

# sc-eQTLgen output loc
sceqtlgen_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/sceqtlgen_eqtls.tsv.gz'
# multiomics eQTL output loc
multiome_eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/multiome_eqtls.tsv.gz'
# multiomics caQTL output loc
multiome_caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/multiome_caqtls.tsv.gz'
# overlapping regions loc
multiome_peak_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/multiome_peaks.tsv.gz'

# here are the sc-eQTLgen results
sceqtlgen_overlap_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta/Out_202406/'

# here are the multiome eQTLs
multiome_eqtl_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/UT/'
multiome_eqtl_24hca_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/24hCA/'

# and here are the multiome caQTLs
multiome_caqtl_ut_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/combined_output_50kb/L1/UT/'
multiome_caqtl_24hca_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/combined_output_50kb/L1/24hCA/'

# here are our bed files
peak_beds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/signac_peaks/output/'


# with a cutoff of 0.8
ld_r2_cutoff <- 0.8

# get the tables
variant1_ld <- fread(variant1_ld_loc, header = T, sep = '\t')
variant2_ld <- fread(variant2_ld_loc, header = T, sep = '\t')
variant3_ld <- fread(variant3_ld_loc, header = T, sep = '\t')

# let's filter first
variant1_ld <- variant1_ld[variant1_ld[['R2']] >= ld_r2_cutoff, ]
variant2_ld <- variant2_ld[variant2_ld[['R2']] >= ld_r2_cutoff, ]
variant3_ld <- variant3_ld[variant3_ld[['R2']] >= ld_r2_cutoff, ]

# add the original variant
variant1_ld <- cbind(data.table('original_snp' = rep('rs3807184', times = nrow(variant1_ld))), variant1_ld)
variant2_ld <- cbind(data.table('original_snp' = rep('rs413524', times = nrow(variant2_ld))), variant2_ld)
variant3_ld <- cbind(data.table('original_snp' = rep('rs113472423', times = nrow(variant3_ld))), variant3_ld)

# we'll keep a set of columns
ld_columns_to_keep <- c('original_snp', 'RS_Number', 'Coord', 'Alleles', 'Distance', 'R2', 'Correlated_Alleles')
# and subset to those
variant1_ld <- variant1_ld[, ..ld_columns_to_keep]
variant2_ld <- variant2_ld[, ..ld_columns_to_keep]
variant3_ld <- variant3_ld[, ..ld_columns_to_keep]
# and we'll give them better names
ld_columns_to_rename <- c('original_snp', 'proxy_snp', 'coordinate', 'alleles', 'distance', 'R2', 'cor_alleles')
# which we'll use
colnames(variant1_ld) <- ld_columns_to_rename
colnames(variant2_ld) <- ld_columns_to_rename
colnames(variant3_ld) <- ld_columns_to_rename

# remove the 'chr' from the coordinages
variant1_ld[['coordinate']] <- gsub('^chr', '', variant1_ld[['coordinate']])
variant2_ld[['coordinate']] <- gsub('^chr', '', variant2_ld[['coordinate']])
variant3_ld[['coordinate']] <- gsub('^chr', '', variant3_ld[['coordinate']])


# get qtls for all variants
sceqtlgen_eqtls <- get_qtls_from_files(sceqtlgen_loc, unique(c(variant1_ld[['coordinate']], variant2_ld[['coordinate']], variant3_ld[['coordinate']])), use_snp_pos_as_name = T)
# there is no need to keep all columns, so let's reduce a bit
qtl_columns_to_deep <- c('cell_type', 'feature_id', 'snp_id', 'p_value', 'beta', 'beta_se', 'assessed_allele', 'maf', 'qval')
sceqtlgen_eqtls <- sceqtlgen_eqtls[, ..qtl_columns_to_deep]
# merge all the variants together
variants_all <- do.call(rbind, list(variant1_ld, variant2_ld, variant3_ld))
# and merge that onto the sc-eqtlgen data
variants_sceqtlgen <- merge(x = variants_all, y = sceqtlgen_eqtls, by.x = 'coordinate', by.y = 'snp_id', all = T)
# order on coordinate
variants_sceqtlgen <- variants_sceqtlgen[order(variants_sceqtlgen[['original_snp']], variants_sceqtlgen[['coordinate']]), ]
# and keep only nominally significant
variants_sceqtlgen <- variants_sceqtlgen[!is.na(variants_sceqtlgen[['p_value']]) & variants_sceqtlgen[['p_value']] < 0.05, ]
# write the result
write.table(variants_sceqtlgen, gzfile(sceqtlgen_output_loc), row.names = F, col.names = T, sep = '\t')
# create the checksum
mdfiver::create_md5_for_file(sceqtlgen_output_loc)

# get the eQTLs in multiome
multiome_eqtls_ut <- get_qtls_foldered(multiome_eqtl_ut_loc, unique(c(variant1_ld[['coordinate']], variant2_ld[['coordinate']], variant3_ld[['coordinate']])), use_snp_pos_as_name = T)
# filter columns
qtl_columns_to_keep_mo <- c('cell_type', 'feature_id', 'snp_id', 'p_value', 'beta', 'beta_se', 'assessed_allele', 'maf', 'feature_q_value')
multiome_eqtls_ut <- multiome_eqtls_ut[, ..qtl_columns_to_keep_mo]
# add the condition specifically as a column
multiome_eqtls_ut <- cbind(data.table('condition' = rep('UT', times = nrow(multiome_eqtls_ut))), multiome_eqtls_ut)
# now do 24hCA as well
multiome_eqtls_24hca <- get_qtls_foldered_full(multiome_eqtl_24hca_loc, unique(c(variant1_ld[['coordinate']], variant2_ld[['coordinate']], variant3_ld[['coordinate']])), use_snp_pos_as_name = T)
# harmonize name of column names
multiome_eqtls_24hca <- cbind(multiome_eqtls_24hca, data.table('feature_q_value' = multiome_eqtls_24hca[['qval']]))
multiome_eqtls_24hca <- multiome_eqtls_24hca[, ..qtl_columns_to_keep_mo]
multiome_eqtls_24hca <- cbind(data.table('condition' = rep('24hCA', times = nrow(multiome_eqtls_24hca))), multiome_eqtls_24hca)
# combine the conditions
multiome_eqtls <- rbind(multiome_eqtls_ut, multiome_eqtls_24hca)
# merge onto our variants again
variants_multiome_eqtls <- merge(x = variants_all, y = multiome_eqtls, by.x = 'coordinate', by.y = 'snp_id', all = T)
# keep only significant ones, in this case removing LD variants with no match in the eQTL data
variants_multiome_eqtls <- variants_multiome_eqtls[!is.na(variants_multiome_eqtls[['p_value']]) & variants_multiome_eqtls[['p_value']] < 0.05, ]
# and write the result
write.table(variants_multiome_eqtls, gzfile(multiome_eqtl_output_loc), row.names = F, col.names = T, sep = '\t')
# create the checksum
mdfiver::create_md5_for_file(multiome_eqtl_output_loc)

# get the caQTLs in multiome
multiome_caqtls_24hca <- get_qtls_foldered_full(multiome_caqtl_24hca_loc, unique(c(variant1_ld[['coordinate']], variant2_ld[['coordinate']], variant3_ld[['coordinate']])), use_snp_pos_as_name = T, sep = ',', verbose = T)
# keep a subset of columns
multiome_caqtls_24hca <- multiome_caqtls_24hca[, ..qtl_columns_to_deep]
# add stimulation condition
multiome_caqtls_24hca <- cbind(data.table('condition' = rep('24hCA', times = nrow(multiome_caqtls_24hca))), multiome_caqtls_24hca)
# repeat for UT
multiome_caqtls_ut <- get_qtls_foldered_full(multiome_caqtl_ut_loc, unique(c(variant1_ld[['coordinate']], variant2_ld[['coordinate']], variant3_ld[['coordinate']])), use_snp_pos_as_name = T, sep = ',', verbose = T)
multiome_caqtls_ut <- multiome_caqtls_ut[, ..qtl_columns_to_deep]
multiome_caqtls_ut <- cbind(data.table('condition' = rep('UT', times = nrow(multiome_caqtls_ut))), multiome_caqtls_ut)
# combine the conditions
multiome_caqtls <- rbind(multiome_caqtls_ut, multiome_caqtls_24hca)
# merge onto our variants again
variants_multiome_caqtls <- merge(x = variants_all, y = multiome_caqtls, by.x = 'coordinate', by.y = 'snp_id', all = T)
# keep only significant ones, in this case removing LD variants with no match in the eQTL data
variants_multiome_caqtls <- variants_multiome_caqtls[!is.na(variants_multiome_caqtls[['p_value']]) & variants_multiome_caqtls[['p_value']] < 0.05, ]
# and write the result
write.table(variants_multiome_eqtls, gzfile(multiome_caqtl_output_loc), row.names = F, col.names = T, sep = '\t')
# create the checksum
mdfiver::create_md5_for_file(multiome_caqtl_output_loc)

# add the position columns
variants_positions <- data.frame('chromosome' = rep(NA, times = nrow(variants_all)), 'start' = rep(NA, times = nrow(variants_all)), 'stop' = rep(NA, times = nrow(variants_all)))
variants_positions[c('chromosome', 'start')] <- str_split_fixed(variants_all$coordinate, ':', 2)
variants_positions[['stop']] <- variants_positions[['start']]
variants_positions[['chrom_string']] <- paste0('chr', variants_positions[['chromosome']])
variants_all <- cbind(variants_all, data.table(variants_positions))
# get the overlaps
variant_overlaps <- check_overlaps_beds(peak_beds_loc, variants_all, variant_chromosome_field = 'chrom_string')
# rename the columns somewhat
colnames(variant_overlaps) <- c('original_snp', 'proxy_snp', 'coordinate', 'alleles', 'distance', 'R2', 'cor_alleles', 'chromosome', 'start', 'stop', 'chrom_string', 'cell_type', 'condition', 'peak_chrom', 'peak_start', 'peak_end', 'peak_name', 'log10reads', 'avg_open', 'ncell', 'frac_wread')
# and write the result
write.table(variant_overlaps, gzfile(multiome_peak_output_loc), row.names = F, col.names = T, sep = '\t')
# create the checksum
mdfiver::create_md5_for_file(multiome_peak_output_loc)
