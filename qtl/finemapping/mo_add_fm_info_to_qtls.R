#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_add_fm_info_to_qtls.R
# Function: add finemapping information to the qtl output file
# Example: 
# Rscript ~/mo_add_fm_info_to_qtls.R \
#  --qtl_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/monocyte/qtl_results_all_qval_allchroms_fdr005_significant.txt.gz \
#  --finemapping_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_caqtl/sc-eqtlgen/output/combined_significant/L1/monocyte_finemapped.tsv.gz \
#  --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/caqtl/sc-eqtlgen/combined_with_qtl/combined/L1/qtl_results_all_qval_allchroms_fdr005_significant_cs.tsv.gz \
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(optparse)
library(mdfiver)


####################
# Functions        #
####################


check_inputs <- function(qtl_table, finemapping_table, qtl_variant_column='snp_id', finemapping_variant_column='variant', qtl_feature_column='feature', finemapping_feature_column='feature', columns_to_add=c('CS')) {
  # check if the columns are present in the qtl table
  qtl_colnames <- colnames(qtl_table)
  if (!(qtl_variant_column %in% qtl_colnames)) {
    stop(paste('qtl_variant_column', qtl_variant_column, 'not present in column names of QTL table. Columns present are', paste(qtl_colnames, collapse=',')))
  }
  if (!(qtl_feature_column %in% qtl_colnames)) {
    stop(paste('qtl_feature_column', qtl_feature_column, 'not present in column names of QTL table. Columns present are', paste(qtl_colnames, collapse=',')))
  }
  # and check if columns are present in the finemapping table
  finemapping_colnames <- colnames(finemapping_table)
  if (!(finemapping_variant_column %in% finemapping_colnames)) {
    stop(paste('finemapping_variant_column', finemapping_variant_column, 'not present in column names of finemapping table. Columns present are', paste(finemapping_colnames, collapse=',')))
  }
  if (!(finemapping_feature_column %in% finemapping_colnames)) {
    stop(paste('finemapping_feature_column', finemapping_feature_column, 'not present in column names of finemapping table. Columns present are', paste(finemapping_colnames, collapse=',')))
  }
  # and the ones to add as well
  for (column_to_add in columns_to_add) {
    if (!(column_to_add %in% finemapping_colnames)) {
      stop(paste('one of the columns to add', column_to_add, 'not present in column names of finemapping table. Columns present are', paste(finemapping_colnames, collapse=',')))
    }
  }
  return(0)
}


####################
# Settings        #
####################


####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-q", "--qtl_file"), type="character", default=NULL,
              help="QTL mapping output file", metavar="character"),
  make_option(c("-f", "--finemapping_file"), type="character", default=NULL,
              help="finemapping output file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-v", "--qtl_variant_column"), type="character", default='snp_id',
              help="column describing the variant in the QTL data [default= %default]", metavar="character"),
  make_option(c("-s", "--finemapping_variant_column"), type="character", default='variant',
              help="column describing the variant in the finemapping data [default= %default]", metavar="character"),
  make_option(c("-g", "--qtl_feature_column"), type="character", default='feature_id',
              help="column describint the feature in the QTL data [default= %default]", metavar="character"),
  make_option(c("-t", "--finemapping_feature_column"), type="character", default='feature',
              help="column describing the feature in the finemapping data [default= %default]", metavar="character"), 
  make_option(c("-c", "--credible_set_columns"), type="character", default='CS,pip',
              help="comma separated string with column names in finemapping to add to QTL output [default= %default]", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
qtl_file <- NULL
finemapping_file <- NULL
out <- NULL

# stop if these are missing
if (is.null(opt[['qtl_file']])) {
  stop(paste('-q/--qtl_file is an obligatory parameter\n'))
} else {
  qtl_file <- opt[['qtl_file']]
}
if (is.null(opt[['finemapping_file']])) {
  stop(paste('-f/--finemapping_file is an obligatory parameter\n'))
} else {
  finemapping_file <- opt[['finemapping_file']]
}
if (is.null(opt[['out']])) {
  stop(paste('-o/--out is an obligatory parameter\n'))
} else {
  out <- opt[['out']]
}

# load options that have defaults
qtl_variant_column <- opt[['qtl_variant_column']]
finemapping_variant_column <- opt[['finemapping_variant_column']]
qtl_feature_column <- opt[['qtl_feature_column']]
finemapping_feature_column <- opt[['finemapping_feature_column']]
credible_set_columns_string <- opt[['credible_set_columns']]

# check if the directory for the output exists, otherwise we would fail at the very last step
out_dir <- dirname(out)
if (!(dir.exists(out_dir))) {
  stop(paste('directory for output file', out, 'does not exist. This would mean the final step writing to file would fail, please create the directory first\n'))
}

# read the first ten lines of the input files
qtl_table <- fread(qtl_file, header = T, sep = '\t', nrows=10)
finemapping_table <- fread(finemapping_file, header = T, sep = '\t', nrows=10)

# extract the credible columns
credible_set_columns <- strsplit(credible_set_columns_string, split = ',')[[1]]

# check the inputs
check_inputs(qtl_table, finemapping_table, qtl_variant_column, finemapping_variant_column, qtl_feature_column, finemapping_feature_column, credible_set_columns)

# if the files are okay, read everything
qtl_table <- fread(qtl_file, header = T, sep = '\t')
finemapping_table <- fread(finemapping_file, header = T, sep = '\t')

# paste together the variant and feature in the qtls
qtl_variant_feature <- paste(qtl_table[[qtl_variant_column]], qtl_table[[qtl_feature_column]], sep = '_')
# and in the finemapping
finemapping_variant_feature <- paste(finemapping_table[[finemapping_variant_column]], finemapping_table[[finemapping_feature_column]], sep = '_')

# extract the columns we care about for the finemapping table
finemapping_relevant_columns <- finemapping_table[, ..credible_set_columns, drop = F]
# get the matching positions of the variant-feature links
finemapping_matching_positions <- match(qtl_variant_feature, finemapping_variant_feature)
# now order the finemapping relevant columns by those positions
finemapping_relevant_columns <- finemapping_relevant_columns[finemapping_matching_positions, , drop = F]
# now add that to the qtl data
qtl_table <- cbind(qtl_table, finemapping_relevant_columns)

# rename the output loc
output_loc <- out
# gz file handle if it ends with .gz
if (grepl('.gz$', output_loc)) {
  # gzip if ends with .gz
  output_loc <- gzfile(output_loc)
}
write.table(qtl_table, output_loc, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum (need non-gzipped file handle for that)
mdfiver::create_md5_for_file(out)
