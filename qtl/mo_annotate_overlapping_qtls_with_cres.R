#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_annotate_overlapping_qtls_with_cres.R
# Function: annotate QTL outputs with DAR or CRE outputs
# Example: Rscript mo_annotate_overlapping_qtls_with_cres.R \
# --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl_caqtl_overlap/combined/L1/all/eqtl_caqtl_overlapping_variants.tsv.gz \
# --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl_caqtl_overlap/combined/L1/all/eqtl_caqtl_overlapping_variants_credar.tsv.gz
#
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
# scenic columns to add
scenic_columns_to_keep <- c('eRegulon_name', 'Gene_signature_name', 'Region_signature_name', 'Gene_signature_direction', 'Region_signature_direction')

# make command line options
option_list <- list(
  make_option(c("-i", "--in"), type="character", default=NULL, 
              help="input QTL file to add information to", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output QTL file to save", metavar="character"), 
  make_option(c("-r", "--region_column"), type="character", default='snp_id', 
              help="column denoting the variant [default: %default]", metavar="character"),
  make_option(c("-f", "--feature_column"), type="character", default='feature_id', 
              help="column denoting the feature [default: %default]", metavar="character"), 
  make_option(c("-e", "--remove_non_overlaps"), action="store_true", default=FALSE,
              help="remove QTL entries that show no overlap with DARs or CREs [default: %default]")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
qtl_in_loc <- NULL
qtl_out_loc <- NULL
region_column <- NULL
feature_column <- NULL
remove_non_overlaps <- NULL

# load debug settings if set to debug mode
if (debug) {
  qtl_in_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl_caqtl_overlap/combined/L1/all/eqtl_caqtl_overlapping_variants.tsv.gz'
  qtl_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl_caqtl_overlap/combined/L1/all/eqtl_caqtl_overlapping_variants_credar.tsv.gz'
  feature_column <- 'feature_eqtl'
  region_column <- 'feature_caqtl'
  remove_non_overlaps <- T
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
  feature_column <- opt[['feature_column']]
  region_column <- opt[['region_column']]
  remove_non_overlaps <- opt[['remove_non_overlaps']]
}

# message
message('loading QTL data')
# load the QTL data
qtl_data <- fread(qtl_in_loc, header = T, sep = '\t')

# message
message('merging with DARs')
# load the DAR data
dars <- fread(dars_output_loc, header = T, sep = '\t')
# check if the region was a DAR
qtl_data[['is_dar']] <- qtl_data[[region_column]] %in% gsub(':', '-', dars[dars[['Adjusted_pval']] < 0.05, ][['region']])
# clear memory
rm(dars)

# message
message('merging with openness data')
# check each cell type for the openness
for (cell_type in openness_cell_types) {
  # we'll paste the path together
  cell_type_openness_loc <- paste0(openness_output_loc, '/', openness_prepend, cell_type, openness_append)
  # let the user know we are reading this data
  message(paste('reading openness file at', cell_type_openness_loc))
  # read the file
  openness_table <- fread(cell_type_openness_loc, header = T, sep = '\t')
  # add openness of cell type in table
  qtl_data[[paste('openness', cell_type, sep = '_')]] <- openness_table[match(qtl_data[[region_column]], openness_table[['name']]), ][['pct_exp']]
}

# load the SCENIC+ data
scenic <- fread(scenic_output_loc, header = T, sep = '\t')

# message
message('merging with SCENIC+ data')

# rename regions in the scenic output
scenic[['Region']] <- gsub(':', '-', scenic[['Region']])

# add region and gene combinations for SCENIC
scenic[['region_to_gene']] <- paste(scenic[['Region']], scenic[['Gene']], sep = '_')
# do the same for the qtl_data
qtl_data[['region_to_gene']] <- paste(qtl_data[[region_column]], qtl_data[[feature_column]], sep = '_')

# subset scenic to the columns I think are important
scenic <- scenic[, c('region_to_gene', ..scenic_columns_to_keep)]
# merge onto the QTL data
qtl_data <- merge(x = qtl_data, y = scenic, by = 'region_to_gene', all.x = T, all.y = F, allow.cartesian=TRUE)

# clear memory
rm(scenic)

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
