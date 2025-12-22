#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lc_export_rna_objects.R
# Function: create an RNA Seurat object to use in LONG-COVID data
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(mdfiver)


####################
# functions        #
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
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 120 * 2000 * 2024^2)

# set seed
set.seed(7777)


####################
# Main Code        #
####################

# location of the expression object
mo_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds'
# read the object
mo <- readRDS(mo_object_loc)

# subset to only the UT samples
mo <- mo[, !is.na(mo@meta.data[['condition_final']]) & mo@meta.data[['condition_final']] == 'UT']
# and keep only the assigned long covid samples
mo <- mo[, !is.na(mo@meta.data[['LONG_COVID_method']]) & mo@meta.data[['LONG_COVID_method']] == 'assigned']
# remove some columns
mo@meta.data[, c('condition_imputed', 'celltype_imputed', 'final_condition', 'condition_sheet', 'condition_prev', 'LONG_COVID_method', 'sample_final')] <- NULL

# create a sample mapping
s_mapping <- create_anonymized_mapping(mo@meta.data, 'realid', 'snumber', 's')
# add these new columns
mo@meta.data[['snumber']] <- s_mapping[match(mo@meta.data[['realid']], s_mapping[['realid']]), 'snumber']
# remove the original annotation
mo@meta.data[['realid']] <- NULL

# add some info
mo@misc[['processed_by']] <- 'r.oelen@umcg.nl'
mo@misc[['generated_by']] <- 'm.g.p.van.der.wijst@umcg.nl'
mo@misc[['exported_at']] <- '2025-12-17'

# where to place the mappings
mo_sample_mapping_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_covid_sample_mapping.tsv.gz'
# and the object
mo_object_covid_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_rna_lc_study_only_20251217.rds'
# export
write.table(s_mapping, gzfile(mo_sample_mapping_loc), row.names = F, col.names = T, sep = '\t', quote = F)
saveRDS(mo, mo_object_covid_loc)
# and make sum checksum
mdfiver::create_sha256_for_file(mo_sample_mapping_loc)
mdfiver::create_sha256_for_file(mo_object_covid_loc)
