#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_anonymize_signac_objects.R
# Function: create anonimyzed version of signac objects
############################################################################################################################

####################
# libraries        #
####################

# read the object
library(Seurat)
library(Signac)

####################
# Functions        #
####################

create_anonymized_mapping <- function(metadata, source_column, target_column, target_prepend='') {
  # get the unique values
  source_values <- unique(as.character(metadata[[source_column]]))
  # turn NA into 'unknown'
  source_values[is.na(source_values)] <- 'unknown'
  # randomly sort your values, by sampling all values
  source_values <- source_values[sample(1:length(source_values), length(source_values))]
  # create a mapping
  target_values <- paste(target_prepend, 1:length(source_values), sep = '')
  # now put that into a dataframe
  mapping_table <- data.frame(x = source_values, y = target_values)
  # now make the column names as expected
  colnames(mapping_table) <- c(source_column, target_column)
  return(mapping_table)
}


####################
# Main Code        #
####################

# get the seurat objects
signac_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/'
signac_object_prepend <- 'mo_cpeaks_filtered_'
signac_object_append <- '_wstatus_1_80_20240709.rds'
# list them
signac_objects <- c(
  'b', 'cd4t', 'cd8t', 'dc', 'monocyte', 'nk'
)
# we'll save the result the in the same location
signac_out_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/'

# check each file
for (ct in signac_objects) {
  # get the full location
  signac_object_full_loc <- paste(signac_object_loc, signac_object_prepend, ct, signac_object_append, sep = '')
  # read the object
  signac_object <- readRDS(signac_object_full_loc)
  # create anonimized mapping
  p_mapping <- create_anonymized_mapping(signac_object@meta.data, 'sample_final', 'sample_number', 's')
  # add these new columns
  signac_object@meta.data[['sample_final']] <- p_mapping[match(signac_object@meta.data[['sample_final']], p_mapping[['sample_final']]), 'sample_number']
  signac_object@meta.data[['best_match_sample']] <- p_mapping[match(signac_object@meta.data[['best_match_sample']], p_mapping[['sample_final']]), 'sample_number']
  signac_object@meta.data[['second_match_sample']] <- p_mapping[match(signac_object@meta.data[['second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
  signac_object@meta.data[['confined_best_match_sample']] <- p_mapping[match(signac_object@meta.data[['confined_second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
  signac_object@meta.data[['confined_second_match_sample']] <- p_mapping[match(signac_object@meta.data[['confined_second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
  signac_object@meta.data[['unconfined_best_match_sample']] <- p_mapping[match(signac_object@meta.data[['unconfined_best_match_sample']], p_mapping[['sample_final']]), 'sample_number']
  signac_object@meta.data[['unconfined_second_match_sample']] <- p_mapping[match(signac_object@meta.data[['unconfined_second_match_sample']], p_mapping[['sample_final']]), 'sample_number']
  # location of the new object
  signac_object_anon_loc <- paste(signac_out_loc, signac_object_prepend, ct, '_wstatus_1_80_20240709_anonymized.rds', sep = '')
  # write the sample mapping
  write.table(p_mapping, gzfile(paste(signac_object_anon_loc, '.sample_mapping.tsv.gz', sep = '')), sep = '\t', row.names = F, col.names = T)
  # remove all other IDs
  signac_object@meta.data[, c('realid')] <- NULL
  # write object
  saveRDS(signac_object, signac_object_anon_loc)
}
