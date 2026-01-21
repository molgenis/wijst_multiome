#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: aging_subset_rna_object.R
# Function: subset the RNA object of the multiomics study to the subset to test topic modelling for age/sex
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(data.table)
library(mdfiver)


####################
# Functions        #
####################


####################
# Settings         #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

####################
# Main Code        #
####################

# object location for Seurat
mo_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds'
# load the object
mo <- readRDS(mo_object_loc)

# the location of the subset file
aging_subset_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/aging_subset_samples.tsv.gz'
# read the subset table
aging_subset <- fread(aging_subset_loc, header = F, sep = '\t')
# set column names
colnames(aging_subset) <- c('lane', 'participant', 'timepoint')

# subset the mo object
mo <- mo[, 
         paste(mo@meta.data[['lane']], mo@meta.data[['sample_final']], mo@meta.data[['condition_final']]) 
         %in%
         paste(aging_subset[['lane']], aging_subset[['participant']], aging_subset[['timepoint']])]

# set location for storing the object
mo_aging_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/aging/objects/rna/aging_mo_subset_20251021.rds'
# save object
saveRDS(mo, mo_aging_loc)
# make checksum
mdfiver::create_sha256_for_file(mo_aging_loc)
