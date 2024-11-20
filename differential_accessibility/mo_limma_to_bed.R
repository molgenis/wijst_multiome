#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_limma_to_bed.R
# Function: convert limma output to bed files (to use in for example scenic plus)
############################################################################################################################

####################
# libraries        #
####################

library(data.table)


####################
# Functions        #
####################

limma_to_bed <- function(limma_output_loc, bed_output_loc, feature_column='feature', significance_column='p.bonferroni', significance_cutoff=0.05, region_to_location=NULL, region_to_location_chrom_column='chrom', region_to_location_start_column='start', region_to_location_end_column='end', region_to_location_feature_column='feature', feature_to_location_split='-') {
  # read the limma file
  limma_output <- NULL
  # if the feature column is NULL, there will be row names
  if (is.null(feature_column)) {
    limma_output <- fread(limma_output_loc, header = T, sep = '\t', row.names = 1)
  }
  else {
    limma_output <- fread(limma_output_loc, header = T, sep = '\t')
  }
  # filter by significance
  if (!is.null(significance_column) & !is.null(significance_cutoff)) {
    limma_output <- limma_output[limma_output[[significance_column]] < significance_cutoff, ]
  }
  # get the significance features
  significant_features <- NULL
  # they will be the rownames if no feature column was supplied
  if (is.null(feature_column)) {
    significant_features <- rownames(limma_output)
  }
  else {
    significant_features <- limma_output[[feature_column]]
  }
  # turn into a bed file
  region_info <- NULL
  # depending on whether a region-to-location file was supplied, we can extract the location from the name, or need to do a lookup
  if (!is.null(region_to_location)) {
    region_info <- region_to_location[match(significant_features, region_to_location[[region_to_location_feature_column]]), c(region_to_location_chrom_column, region_to_location_start_column, region_to_location_end_column)]
  }
  # otherwise do based on character split
  else {
    region_info <- data.frame(do.call(rbind, strsplit(significant_features, split = feature_to_location_split)))
  }
  # set the column names
  colnames(region_info) <- c('#chrom', 'chromStart', 'chromEnd')
  # add the features with their name
  region_info[['name']] <- significant_features
  # write the result
  output_loc_gzipped <- bed_output_loc
  # zipped if necessary
  if (grepl('.gz$', bed_output_loc)) {
    output_loc_gzipped <- gzfile(bed_output_loc)
  }
  write.table(region_info, output_loc_gzipped, row.names = F, col.names = T, sep = '\t', quote = F)
  return(0)
}


limmas_to_beds <- function(limma_dir, bed_dir, limma_prepend='', limma_append='.nominal_significant.tsv.gz', bed_prepend='', bed_append='.bed', feature_column='feature', significance_column='p.bonferroni', significance_cutoff=0.05, region_to_location=NULL, region_to_location_chrom_column='chrom', region_to_location_start_column='start', region_to_location_end_column='end', region_to_location_feature_column='feature', feature_to_location_split='-') {
  # list the files in the directory
  limma_files <- list.files(limma_dir, recursive = F)
  # now get a regex of what the files should look like
  limma_file_pattern <- paste(limma_prepend, '*', limma_append, '$', sep = '')
  # filter the files by that pattern
  limma_files <- limma_files[grepl(limma_file_pattern, limma_files)]
  # go through each of these files
  for (limma_file in limma_files) {
    # remove the prepend and append
    limma_file_base <- gsub(limma_prepend, '', limma_file)
    limma_file_base <- gsub(limma_append, '', limma_file_base)
    # add the bed path, prepend and append
    bed_file_full <- paste(bed_dir, '/', bed_prepend, limma_file_base, bed_append, sep = '')
    # get the full path to the limma file
    limma_file_full <- paste(limma_dir, '/', limma_file, sep = '')
    # now do the actual conversion
    limma_to_bed(limma_output_loc = limma_file_full, 
                 bed_output_loc = bed_file_full, 
                 feature_column = feature_column, 
                 significance_column = significance_column, 
                 significance_cutoff = significance_cutoff, 
                 region_to_location = region_to_location, 
                 region_to_location_chrom_column = region_to_location_chrom_column, 
                 region_to_location_start_column = region_to_location_start_column, 
                 region_to_location_end_column = region_to_location_end_column, 
                 region_to_location_feature_column = region_to_location_feature_column, 
                 feature_to_location_split = feature_to_location_split)
      
  }
  return(0)
}


####################
# Settings         #
####################


####################
# Main Code        #
####################

# convert the monocyte data
limmas_to_beds(
  limma_dir = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/topics20_otsu_imputed_confined_ncell5/monocyte/merged_regions/', 
  bed_dir = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/topics20_otsu_imputed_confined_ncell5/monocyte/merged_regions/'
)
