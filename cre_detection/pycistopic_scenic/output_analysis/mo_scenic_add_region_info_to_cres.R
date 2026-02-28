#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_scenic_add_region_info_to_cres.R
# Function: add information about the regions to the cre output of scenic
############################################################################################################################

####################
# libraries        #
####################

library(r2r)
library(data.table)
library(mdfiver)
library(stringr)


####################
# Functions        #
####################

#' Get Regions to Topics Mapping
#'
#' This function creates a hashmap that maps genomic regions to topics based on BED files located in a specified directory.
#'
#' @param topic_dar_bed_loc Character. The directory location containing the topic BED files.
#' @param region_column Character. The name of the column containing region information. Default is 'Region'.
#' @param chrom_sep Character. The separator between chromosome and start position. Default is ':'.
#' @param region_sep Character. The separator between start and end positions. Default is '-'.
#' @param interested_regions Character vector. A vector of regions to filter the results. Default is NULL.
#'
#' @return A hashmap where keys are regions and values are topics.
#'
#' @examples
#' \dontrun{
#' topic_dar_bed_loc <- "path/to/bed/files"
#' region_to_topics <- get_regions_to_topics(topic_dar_bed_loc)
#' }
get_regions_to_topics <- function(topic_dar_bed_loc, region_column='Region', chrom_sep=':', region_sep='-', interested_regions=NULL) {
  # create hashmap of region to topics
  region_to_topics <- r2r::hashmap()
  # list the files in the topic directory
  topic_files <- list.files(topic_dar_bed_loc, recursive = F, full.names = F)
  # we will only keep the bed files
  topic_file_regex <- '.*\\.bed'
  # filter on those
  topic_files <- topic_files[grep(topic_file_regex, topic_files)]
  # check each file
  for (topic_file in topic_files) {
    # read the topic files
    topic_table <- fread(paste0(topic_dar_bed_loc, '/',  topic_file), header = F, sep = '\t')
    # create a new region name vector
    topic_table_regions <- paste0(topic_table[['V1']], chrom_sep, topic_table[['V2']], region_sep, topic_table[['V3']])
    # extract topic from file name
    topic_name <- gsub('\\.bed', '', topic_file)
    # subset regions we are interested in
    if (!is.null(interested_regions)) {
      topic_table_regions <- topic_table_regions[topic_table_regions %in% interested_regions]
    }
    # check each region
    for (region in topic_table_regions) {
      # set if this does not exist
      if (!r2r::has_key(region_to_topics, region)) {
        region_to_topics[[region]] <- topic_name
      }
      else {
        region_to_topics[[region]] <- paste(region_to_topics[[region]], topic_name, sep = ',')
      }
    }
  }
  return(region_to_topics)
}

#' Convert r2r Hashmap to Data Table
#'
#' This function converts an r2r hashmap to a data.table object.
#'
#' @param r2r_hasmap An r2r hashmap object.
#'
#' @return A data.table with two columns: 'r2r_key' and 'r2r_values'.
#'
#' @examples
#' \dontrun{
#' r2r_hasmap <- r2r::hashmap()
#' r2r_hasmap[['key1']] <- 'value1'
#' r2r_hasmap[['key2']] <- 'value2'
#' keyvalue_dt <- r2r_to_datatable(r2r_hasmap)
#' }
r2r_to_datatable <- function(r2r_hasmap) {
  # get the keys
  r2r_keys <- unlist(keys(r2r_hasmap))
  # we'll get the values as well
  r2r_values <- rep(NA, times = length(r2r_keys))
  # check each key
  for (i in 1 : length(r2r_keys)) {
    # grab the i
    r2r_key <- r2r_keys[i]
    # and set the corresponding value
    r2r_value <- r2r_hasmap[[r2r_key]]
    r2r_values[i] <- r2r_value
  }
  # create a table
  keyvalue_dt <- data.table(data.frame('r2r_key' = r2r_keys, 'r2r_values' = r2r_values))
  return(keyvalue_dt)
}


#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @param gene_column which column to use as the gene identifier
#' @returns a list with the egenes per cell type
#' 
get_egenes_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all.txt.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, verbose=T) {
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
    }
    # get the unique genes in this file
    unique_genes <- unique(cell_type_output[[gene_column]])
    # add to the list
    egenes_per_celltype[[cell_type]] <- unique_genes
  }
  # turn into a dataframe
  return(egenes_per_celltype)
}


add_region_qtl_info <- function(table_to_annotate, group_to_regions_list, to_annotate_region_column='Region', column_to_add='caqtl_celltypes') {
  # extrac the unique regions
  regions <- table_to_annotate[[to_annotate_region_column]]
  # create table to put results
  region_to_groups <- data.frame('region' = regions)
  # add the groups to the table
  region_to_groups[[column_to_add]] <- apply(region_to_groups, 1, function(x) {
    # grab the region
    region <- gsub(':', '-', x[['region']])
    # initialize the results
    groups_string <- NA
    # check each of the vectors in the list
    for (group in names(group_to_regions_list)) {
      # check if the region is in there
      if (region %in% group_to_regions_list[[group]]) {
        # if the value was NA, it will be this one
        if (is.na(groups_string)) {
          groups_string <- group
        }
        # otherwise add it
        else {
          groups_string <- paste(groups_string, group, sep = ',')
        }
      }
    }
    return(groups_string)
  })
  # add this info to the table
  table_to_annotate[[column_to_add]] <- region_to_groups[match(table_to_annotate[[to_annotate_region_column]], region_to_groups[['region']]), ][[column_to_add]]
  return(table_to_annotate)
}


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


####################
# Main Code        #
####################

# location of the topic files for membership after binarization
binarized_topic_beds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/region_sets/topic_memberships/'
# location of the topic dars
dar_topic_beds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/region_sets/topic_dars/'

# location of the direct file scenic output
scenic_output_direct_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_direct.tsv'
# read the scenic output
scenic_output_direct <- read.table(scenic_output_direct_loc, header = T, sep = '\t')
# location of the extended file scenic output
scenic_output_extended_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulons_extended.tsv'
# read the scenic output
scenic_output_extended <- read.table(scenic_output_extended_loc, header = T, sep = '\t')

# add where the info came from
scenic_output_extended[['source']] <- 'extended'
scenic_output_direct[['source']] <- 'direct'
# merge them
scenic_output <- rbind(scenic_output_direct, scenic_output_extended)

# add the signatures that come from the direction of the region to tf and gene
scenic_output[['Gene_signature_direction']] <- str_extract(scenic_output[['Gene_signature_name']], '\\+\\/\\+|\\-\\/\\-|\\+\\/\\-|\\-\\/\\+')
scenic_output[['Region_signature_direction']] <- str_extract(scenic_output[['Region_signature_name']], '\\+\\/\\+|\\-\\/\\-|\\+\\/\\-|\\-\\/\\+')

# get the regions to the topics of the binarization
regions_to_topic_memberships <- get_regions_to_topics(topic_dar_bed_loc = binarized_topic_beds_loc, interested_regions = scenic_output[['Region']])
# make into data.table
regions_to_topic_memberships_dt <- r2r_to_datatable(regions_to_topic_memberships)
# join onto the table
scenic_output[['topics_membership']] <- regions_to_topic_memberships_dt[match(scenic_output[['Region']], regions_to_topic_memberships_dt[['r2r_key']]), 'r2r_values'][[1]]

# get the regions to the topics of the binarization
regions_to_topics_dars <- get_regions_to_topics(topic_dar_bed_loc = dar_topic_beds_loc, interested_regions = scenic_output[['Region']])
# make into data.table
regions_to_topics_dars_dt <- r2r_to_datatable(regions_to_topics_dars)
# join onto the table
scenic_output[['topics_dar']] <- regions_to_topics_dars_dt[match(scenic_output[['Region']], regions_to_topics_dars_dt[['r2r_key']]), 'r2r_values'][[1]]

# location of eQTL
eqtl_output_combined_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# location of caQTL output
caqtl_output_combined_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'
# get regions and genes that have eGenes or ePeaks
eqtl_egenes_combined <- get_egenes_per_celltype_limix(eqtl_output_combined_loc, output_file = 'qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')
caqtl_epeaks_combined <- get_egenes_per_celltype_limix(caqtl_output_combined_loc, output_file = 'qtl_results_all_qval_allchroms_fdr005_significant.txt.gz')

# add the QTL info
scenic_output <- add_region_qtl_info(table_to_annotate = scenic_output, group_to_regions_list = caqtl_epeaks_combined, to_annotate_region_column = 'Region', column_to_add = 'caqtl_celltype')
scenic_output <- add_region_qtl_info(table_to_annotate = scenic_output, group_to_regions_list = eqtl_egenes_combined, to_annotate_region_column = 'Gene', column_to_add = 'eqtl_celltype')

# get all the QTL output
eqtl_output <- get_qtls_per_celltype_limix(eqtl_output_combined_loc)
# merge all the results
eqtl_output_all <- do.call('rbind', eqtl_output)
# filter on significance
eqtl_output_all_sig <- eqtl_output_all[eqtl_output_all[['feature_q_value']] < 0.05 &
                                        eqtl_output_all[['p_value']] < eqtl_output_all[['pval_nominal_threshold_global']], ]
# get the locations of the genetic variants
qtl_variants_all_cpeaks_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
# and read that
qtl_variants_all_cpeaks <- fread(qtl_variants_all_cpeaks_loc, header = T, sep = '\t')
# now add the region info for each variant
eqtl_output_all_sig[['variant_region']] <- qtl_variants_all_cpeaks[match(eqtl_output_all_sig[['snp_id']], qtl_variants_all_cpeaks[['snp_id']]), ][['overlapping_feature']]
# make a list to store each region-gene pair per cell type
eqtl_r2g_per_celltype <- list()
# check each cell type
for (cell_type in unique(eqtl_output_all_sig[['cell_type']])) {
  # extract unique entries for cell type
  eqtl_output_celltype <- eqtl_output_all_sig[eqtl_output_all_sig[['cell_type']] == cell_type &
                                                !is.na(eqtl_output_all_sig[['variant_region']]), ]
  # get the unique region-gene pairs
  eqtl_r2g_per_celltype[[cell_type]] <- unique(paste(eqtl_output_celltype[['variant_region']], eqtl_output_celltype[['feature_id']], sep = '_'))
}
# add r2g
scenic_output[['r2g']] <- paste(scenic_output[['Region']], scenic_output[['Gene']], sep = '_')
# and add the variant-region info
scenic_output <- add_region_qtl_info(table_to_annotate = scenic_output, group_to_regions_list = eqtl_r2g_per_celltype, to_annotate_region_column = 'r2g', column_to_add = 'eqtl_varregion_celltype')

# save this file somewhere
scenic_output_region_info <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
write.table(scenic_output, gzfile(scenic_output_region_info), row.names = F, col.names = T, sep = '\t')
# and make a checksum
mdfiver::create_md5_for_file(scenic_output_region_info)

# check now many effects are under genetic regulation
n_cre_genreg <- nrow(scenic_output[!is.na(scenic_output[['caqtl_celltype']]) | !is.na(scenic_output[['eqtl_varregion_celltype']]), ])
# and the proportion
prop_cre_genreg <- n_cre_genreg / nrow(scenic_output)
# [1] 0.1375937

# get the total number of eGenes
n_egenes <- length(unique(do.call('c', eqtl_egenes_combined)))
n_ePeaks <- length(unique(do.call('c', caqtl_epeaks_combined)))

# make a table
scenic_qtl_overlap <- data.frame(
  'set' = c('cre_total', 'cre_with_caqtl', 'cre_with_eqtl', 'cre_with_any_qtl', 'cre_with_both_qtl', 'egene_total', 'epeak_total'), 
  'number' = c(
    nrow(scenic_output), 
    nrow(scenic_output[!is.na(scenic_output[['caqtl_celltype']]), ]), 
    nrow(scenic_output[!is.na(scenic_output[['eqtl_celltype']]), ]), 
    nrow(scenic_output[!is.na(scenic_output[['eqtl_celltype']]) | !is.na(scenic_output[['caqtl_celltype']]), ]), 
    nrow(scenic_output[!is.na(scenic_output[['eqtl_celltype']]) & !is.na(scenic_output[['caqtl_celltype']]), ]), 
    n_egenes, 
    n_ePeaks
  )
)
# now if we do this without considering the TFs
scenic_output_onlyregions <- unique(scenic_output[, c('Region', 'Gene', 'caqtl_celltype', 'eqtl_celltype')])
scenic_qtl_overlap_notf <- data.frame(
  'set' = c('cre_total', 'cre_with_caqtl', 'cre_with_eqtl', 'cre_with_any_qtl', 'cre_with_both_qtl', 'egene_total', 'epeak_total'), 
  'number' = c(
    nrow(scenic_output_onlyregions), 
    nrow(scenic_output_onlyregions[!is.na(scenic_output_onlyregions[['caqtl_celltype']]), ]), 
    nrow(scenic_output_onlyregions[!is.na(scenic_output_onlyregions[['eqtl_celltype']]), ]), 
    nrow(scenic_output_onlyregions[!is.na(scenic_output_onlyregions[['eqtl_celltype']]) | !is.na(scenic_output_onlyregions[['caqtl_celltype']]), ]), 
    nrow(scenic_output_onlyregions[!is.na(scenic_output_onlyregions[['eqtl_celltype']]) & !is.na(scenic_output_onlyregions[['caqtl_celltype']]), ]), 
    n_egenes, 
    n_ePeaks
  )
)
