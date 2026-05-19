#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_hybrid_cre_interaction_merge.R
# Function: merge perform interaction-eQTL analysis at single-cell level result files
# Example: 
# ~/start_Rscript.sh \
#   /groups/umcg-franke-scrna/tmp02/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_hybrid_cre_interaction_merge.R \
#   --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna/ \
#   --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna/merged/results_fdr.tsv.gz
# 
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)
# load command line parameters
library(optparse)


####################
# Functions        #
####################

write_empty_result <- function(output_loc) {
  # gz file ends with .gz
  if (grepl('.gz$', output_loc)) {
    # gzip if ends with .gz
    con <- gzfile(output_loc, 'w')
    close(con)
  }
  else {
    file.create(output_loc)
  }
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

# make command line options
option_list <- list(
  make_option(c("-i", "--in"), type="character", default=NULL, 
              help="input directory of chunks", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file", metavar="character"), 
  make_option(c("-p", "--pattern"), type="character", default='chr\\d+\\-\\d+\\-\\d+', 
              help="pattern for grabbing chunks", metavar="character")
)


# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# input directory
in_dir <- NULL
# location of the output
output_loc <- NULL
# get the chunk pattern
chunk_pattern <- NULL

if (debug) {
  # set all of the variables hardcoded for a testing debug run
  in_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor/'
  output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor/merged/results_fdr.tsv.gz'
  # get the chunk pattern
  chunk_pattern <- 'chr\\d+\\-\\d+\\-\\d+'
  
} else {
  # obligatory parameters without a default
  if (is.null(opt[['in']])) {
    error("i/--in is an obligatory parameter")
  } else {
    in_dir <- opt[['in']]
  }
  if (is.null(opt[['out']])) {
    error("o/--out is an obligatory parameter")
  } else {
    output_loc <- opt[['out']]
  }
  # get the chunk pattern
  chunk_pattern <- opt[['pattern']]
}

# get the chunk pattern
# chunk_pattern <- 'chr\\d+\\-\\d+\\-\\d+'
# and the name of each file
chunk_filename <-  'result.tsv.gz'

# list the directories present
dirs_present <- list.dirs(in_dir, recursive = F)
# get the base names
dirs_basenames <- basename(dirs_present)

# filter the dirs on chunks
chunk_fulldirs <- dirs_present[grep(chunk_pattern, dirs_basenames)]
# order for replication sake
chunk_fulldirs <- chunk_fulldirs[order(chunk_fulldirs)]
# add the full paths as well
chunk_fullpaths <- paste(chunk_fulldirs, chunk_filename, sep = '/')
# make a list to keep each of the tables
chunk_tables <- list()
# keep track of all the column
chunk_cols <- list()
# check each chunk
for (chunk_path in chunk_fullpaths) {
  # check if the file is empty
  if (length(count.fields(chunk_path)) > 1) {
    # read each chunk
    chunk_contents <- fread(chunk_path, header = T, sep = '\t')
    # check if we have rows
    if (nrow(chunk_contents) > 0) {
      # place in the list
      chunk_tables[[chunk_path]] <- chunk_contents
      # add to the list of columns
      chunk_cols[[chunk_path]] <- colnames(chunk_contents)
    }
  }
}
# get the outer join of those columns
all_columns <- unique(do.call('c', chunk_cols))
# make sure all tables have all columns, and in the same order
for (chunk_path in names(chunk_tables)) {
  # get the columns of this table
  chunk_cols_i <- chunk_cols[[chunk_path]]
  # get the missing columns
  missing_cols <- setdiff(all_columns, chunk_cols_i)
  # add those missing columns as NA
  for (missing_col in missing_cols) {
    chunk_tables[[chunk_path]][[missing_col]] <- NA
  }
  # order the columns in the same way for all tables
  chunk_tables[[chunk_path]] <- chunk_tables[[chunk_path]][, all_columns, with = F]
}
# merge all
chunks_merged <- rbindlist(chunk_tables, fill = T)

# write empty file if there are no results
if (is.null(chunks_merged)) {
  warning('merging chunks gave completely empty results')
  write_empty_result(output_loc)
} else {
  # extract each column that has a p value
  p_columns <- colnames(chunks_merged)[grep('_p$', colnames(chunks_merged))]
  # do MTC for each of these columns
  for (p_column in p_columns) {
    # check which rows have a value for this p
    p_valid_i <- !is.na(chunks_merged[[p_column]])
    # add MTC column
    bh_column <- gsub('_p$', '_bh', p_column)
    chunks_merged[[bh_column]] <- NA
    # extract valid p values, and calculate B&H
    bhs <- p.adjust(chunks_merged[[p_column]][p_valid_i], method = 'BH')
    # then place those BHs
    chunks_merged[[bh_column]][p_valid_i] <- bhs
    # add bf column as well
    bf_column <- gsub('_p$', '_bf', p_column)
    # adjust the valid p values
    bfs <- p.adjust(chunks_merged[[p_column]][p_valid_i], method = 'bonferroni')
    # add the bonferroni column
    chunks_merged[[bf_column]] <- NA
    # and place the bonferroni adjusted p values
    chunks_merged[[bf_column]][p_valid_i] <- bfs
  }
  # set output loc as the tsv
  output_loc_full <- output_loc
  # gz file ends with .gz
  if (grepl('.gz$', output_loc)) {
    # gzip if ends with .gz
    output_loc_full <- gzfile(output_loc)
  }
  # remove any columns that are fully NA
  chunks_merged <- chunks_merged[, which(unlist(lapply(chunks_merged, function(x) !all(is.na(x))))), with = F]
  # write table
  write.table(chunks_merged, output_loc_full, row.names = F, col.names = T, sep = '\t', quote = F)
  # make checksum
  mdfiver::create_sha256_for_file(output_loc)
}

