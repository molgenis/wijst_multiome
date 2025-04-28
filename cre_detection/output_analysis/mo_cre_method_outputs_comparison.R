#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_cre_method_outputs_comparison.R
# Function: compare the cre outputs of the different methods
############################################################################################################################

####################
# libraries        #
####################

library(r2r)
library(data.table)
library(mdfiver)

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


####################
# Main Code        #
####################

# location of the topic files
topic_beds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_limma_001/region_sets/model_20_topics/Topics_otsu/'
# location of the file scenic output
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_limma_001/output/eRegulon_direct.tsv'
# read the scenic output
scenic_outpu <- read.table(scenic_output_loc, header = T, sep = '\t')
# get the regions to the topics
regions_to_topics <- get_regions_to_topics(topic_dar_bed_loc = topic_beds_loc, interested_regions = scenic_outpu[['Region']])
# make into data.table
regions_to_topics_dt <- r2r_to_datatable(regions_to_topics)
# join onto the table
scenic_outpu[['topics']] <- regions_to_topics_dt[match(scenic_outpu[['Region']], regions_to_topics_dt[['r2r_key']]), 'r2r_values'][[1]]
# location of the coeqtl-esque method
coeqtl_cre_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/pseudobulk_replication/betas_ps/monocyte/meta_analysis_results.txt.gz'
coeqtl_cre_output <- read.table(coeqtl_cre_output_loc, header = T, sep = '\t')
# see what replicated
scenic_outpu[['coeqtl_rep']] <- paste(scenic_outpu[['Gene']], gsub(':', '-', scenic_outpu[['Region']]), sep = '-') %in% rownames(coeqtl_cre_output[coeqtl_cre_output[['meta_q_value']] < 0.05, ])
# location of the significant pseudobulk CREs
significant_pseudobulks_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/coeqtl/trial_run/cre_lists/monocyte_pseudobulk_cres.tsv.gz'
significant_pseudobulks <- read.table(significant_pseudobulks_loc, header = F, sep = '\t')
# see what replicated
scenic_outpu[['pseudo_rep']] <- paste(scenic_outpu[['Gene']], gsub(':', '-', scenic_outpu[['Region']]), sep = '-') %in% paste(significant_pseudobulks[['V1']], significant_pseudobulks[['V2']], sep = '-')
# store the result
scenic_annotated_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_limma_001/output/eRegulon_direct_overlaps.tsv.gz'
write.table(scenic_outpu, gzfile(scenic_annotated_loc), row.names = F, col.names = T, quote = F, sep = '\t')
mdfiver::create_md5_for_file(scenic_annotated_loc)
