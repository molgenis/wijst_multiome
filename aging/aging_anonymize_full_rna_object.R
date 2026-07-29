#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: aging_anonymize_full_rna_object.R
# Function: 
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(mdfiver)
library(R.utils)


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
mo <- readRDS(mo_object_loc)# remove some columns
mo@meta.data[, c('realid', 'final_condition', 'condition_imputed', 'celltype_imputed', 'condition_sheet', 'condition_prev', 'LONG_COVID_method', 'best_match_sample', 'second_match_sample', 'confined_best_match_sample', 'confined_second_match_sample', 'unconfined_best_match_sample', 'unconfined_second_match_sample')] <- NULL

# create a sample mapping
s_mapping <- create_anonymized_mapping(mo@meta.data, 'sample_final', 'snumber', 's')
# add these new columns
mo@meta.data[['snumber']] <- s_mapping[match(mo@meta.data[['sample_final']], s_mapping[['sample_final']]), 'snumber']
# remove the original annotation
mo@meta.data[['sample_final']] <- NULL

# add some info
mo@misc[['processed_by']] <- 'r.oelen@umcg.nl'
mo@misc[['generated_by']] <- 'm.g.p.van.der.wijst@umcg.nl'
mo@misc[['exported_at']] <- '2026-07-20'

# write the mapping table
mapping_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid_smapping.tsv.gz'
write.table(s_mapping, gzfile(mapping_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(mapping_loc)

# now write the object
mo_object_anon_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid_anon.rds'
saveRDS(mo, mo_object_anon_loc)
mdfiver::create_sha256_for_file(mo_object_anon_loc)
