#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Drew Neavin
# Name: mo_correlate_genotypes.R
# Function: correlate the souporcell output to gsa genotypes
############################################################################################################################


####################
# libraries        #
####################

# required do correlations
library(tidyr)
library(tidyverse)
library(dplyr)
library(vcfR)
# library to read the object
library(Seurat)
# libraries to do plots
library(ggplot2)

####################
# Functions        #
####################


pearson_correlation <- function(df, ref_df, clust_df){
  for (col in colnames(df)){
    for (row in rownames(df)){
      df[row,col] <- cor(as.numeric(pull(ref_df, col)), as.numeric(pull(clust_df, row)), method = "pearson", use = "complete.obs")
    }
  }
  return(df)
}

# method taken from https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/includes/Snakefile_souporcell.smk
correlate_genotypes <- function(ref_geno_loc, cluster_geno_loc){
  ref_geno <- read.vcfR(ref_geno_loc)
  cluster_geno <- read.vcfR(cluster_geno_loc)
  
  ########## Convert to tidy data frame ##########
  ref_geno_tidy <- as_tibble(extract.gt(element = "DS",ref_geno, IDtoRowNames =F))
  ref_geno_tidy$ID <- paste0(ref_geno@fix[,'CHROM'],":", ref_geno@fix[,'POS'],"_", ref_geno@fix[,'REF'], "_",ref_geno@fix[,'ALT'])
  ref_geno_tidy <- ref_geno_tidy[!(ref_geno_tidy$ID %in% ref_geno_tidy$ID[duplicated(ref_geno_tidy$ID)]),]
  
  cluster_geno_tidy <- as_tibble(extract.gt(element = "GT",cluster_geno, IDtoRowNames =F))
  cluster_geno_tidy <- as_tibble(lapply(cluster_geno_tidy, function(x) {gsub("0/0",0, x)}) %>%
                                   lapply(., function(x) {gsub("0/1",1, x)}) %>%
                                   lapply(., function(x) {gsub("1/0",1, x)}) %>%
                                   lapply(., function(x) {gsub("1/1",2, x)}))
  cluster_geno_tidy$ID <- paste0(cluster_geno@fix[,'CHROM'],":", cluster_geno@fix[,'POS'],"_", cluster_geno@fix[,'REF'], "_",cluster_geno@fix[,'ALT'])
  cluster_geno_tidy <- cluster_geno_tidy[colSums(!is.na(cluster_geno_tidy)) > 0]
  # cluster_geno_tidy <- cluster_geno_tidy[complete.cases(cluster_geno_tidy),]
  cluster_geno_tidy <- cluster_geno_tidy[!(cluster_geno_tidy$ID %in% cluster_geno_tidy$ID[duplicated(cluster_geno_tidy$ID)]),]
  
  
  ########## Get a unique list of SNPs that is in both the reference and cluster genotypes ##########
  locations  <- inner_join(ref_geno_tidy[,"ID"],cluster_geno_tidy[,"ID"])
  locations <- locations[!(locations$ID %in% locations[duplicated(locations),"ID"]),]
  
  ########## Keep just the SNPs that overlap ##########
  ref_geno_tidy <- left_join(locations, ref_geno_tidy)
  cluster_geno_tidy <- left_join(locations, cluster_geno_tidy)
  
  ########## Correlate all the cluster genotypes with the individuals genotyped ##########
  ##### Make a dataframe that has the clusters as the row names and the individuals as the column names #####
  pearson_correlations <- as.data.frame(matrix(nrow = (ncol(cluster_geno_tidy) -1), ncol = (ncol(ref_geno_tidy) -1)))
  colnames(pearson_correlations) <- colnames(ref_geno_tidy)[2:(ncol(ref_geno_tidy))]
  rownames(pearson_correlations) <- colnames(cluster_geno_tidy)[2:(ncol(cluster_geno_tidy))]
  pearson_correlations <- pearson_correlation(pearson_correlations, ref_geno_tidy, cluster_geno_tidy)
  
  return(pearson_correlations)
}


correlate_cluster_genotypes <- function(cluster_geno_loc1, cluster_geno_loc2){
  # read lane1
  cluster_geno1 <- read.vcfR(cluster_geno_loc1)
  cluster_geno1_tidy <- as_tibble(extract.gt(element = "GT",cluster_geno1, IDtoRowNames =F))
  cluster_geno1_tidy <- as_tibble(lapply(cluster_geno1_tidy, function(x) {gsub("0/0",0, x)}) %>%
                                    lapply(., function(x) {gsub("0/1",1, x)}) %>%
                                    lapply(., function(x) {gsub("1/0",1, x)}) %>%
                                    lapply(., function(x) {gsub("1/1",2, x)}))
  cluster_geno1_tidy$ID <- paste0(cluster_geno1@fix[,'CHROM'],":", cluster_geno1@fix[,'POS'],"_", cluster_geno1@fix[,'REF'], "_",cluster_geno1@fix[,'ALT'])
  cluster_geno1_tidy <- cluster_geno1_tidy[colSums(!is.na(cluster_geno1_tidy)) > 0]
  cluster_geno1_tidy <- cluster_geno1_tidy[!(cluster_geno1_tidy$ID %in% cluster_geno1_tidy$ID[duplicated(cluster_geno1_tidy$ID)]),]
  # read lane2
  cluster_geno2 <- read.vcfR(cluster_geno_loc2)
  cluster_geno2_tidy <- as_tibble(extract.gt(element = "GT",cluster_geno2, IDtoRowNames =F))
  cluster_geno2_tidy <- as_tibble(lapply(cluster_geno2_tidy, function(x) {gsub("0/0",0, x)}) %>%
                                    lapply(., function(x) {gsub("0/1",1, x)}) %>%
                                    lapply(., function(x) {gsub("1/0",1, x)}) %>%
                                    lapply(., function(x) {gsub("1/1",2, x)}))
  cluster_geno2_tidy$ID <- paste0(cluster_geno2@fix[,'CHROM'],":", cluster_geno2@fix[,'POS'],"_", cluster_geno2@fix[,'REF'], "_",cluster_geno2@fix[,'ALT'])
  cluster_geno2_tidy <- cluster_geno2_tidy[colSums(!is.na(cluster_geno2_tidy)) > 0]
  cluster_geno2_tidy <- cluster_geno2_tidy[!(cluster_geno2_tidy$ID %in% cluster_geno2_tidy$ID[duplicated(cluster_geno2_tidy$ID)]),]
  
  ########## Keep just the SNPs that overlap ##########
  locations  <- inner_join(cluster_geno1_tidy[,"ID"],cluster_geno2_tidy[,"ID"])
  locations <- locations[!(locations$ID %in% locations[duplicated(locations),"ID"]),]
  cluster_geno1_tidy <- left_join(locations, cluster_geno1_tidy)
  cluster_geno2_tidy <- left_join(locations, cluster_geno2_tidy)
  ##### Make a dataframe that has the clusters as the row names and the individuals as the column names #####
  pearson_correlations <- as.data.frame(matrix(ncol = (ncol(cluster_geno1_tidy) -1), nrow = (ncol(cluster_geno2_tidy) -1)))
  colnames(pearson_correlations) <- colnames(cluster_geno1_tidy)[2:(ncol(cluster_geno1_tidy))]
  rownames(pearson_correlations) <- colnames(cluster_geno2_tidy)[2:(ncol(cluster_geno2_tidy))]
  pearson_correlations <- pearson_correlation(pearson_correlations, cluster_geno1_tidy, cluster_geno2_tidy)
  return(pearson_correlations)
}

get_correlation_matrix_per_lane <- function(souporcell_output_loc, genotypes_loc, lanes, genotype_prepend='mo_', genotype_append='_maf005.vcf'){
  # we will save the correlations per lane
  correlations_per_lane <- list()
  # check each lane
  for(lane in lanes){
    print(paste('correlating lane', lane))
    # paste together the souporcell cluster loc
    cluster_loc <- paste(souporcell_output_loc, lane, '/cluster_genotypes.vcf', sep = '')
    # paste together the genotype loc
    geno_loc <- paste(genotypes_loc, genotype_prepend, lane, genotype_append, sep = '')
    # get the result
    correlations_lane <- correlate_genotypes(geno_loc, cluster_loc)
    # store the result in the list
    correlations_per_lane[[lane]] <- correlations_lane
  }
  return(correlations_per_lane)
}


get_best_correlations <- function(correlations_per_lane){
  # store the best correlation per lane
  best_correlation_per_lane <- NULL
  # check each lane
  for(lane in names(correlations_per_lane)){
    # get that specific data
    correlation_table <- correlations_per_lane[[lane]]
    # check each cluster
    for(cluster in rownames(correlation_table)){
      # get the values for that row
      correlations <- as.vector(unlist(correlation_table[cluster, ]))
      # order them
      order_correlations <- order(correlations, decreasing = T)
      # get the location of the best match
      best_loc <- order_correlations[1]
      # and the second best match
      second_loc <- order_correlations[2]
      # get the colnames, and thus matches for the best and second
      best_match <- colnames(correlation_table)[best_loc]
      second_match <- colnames(correlation_table)[second_loc]
      # and get the correlations we found for those
      best_correlation <- correlation_table[cluster, best_match]
      second_correlation <- correlation_table[cluster, second_match]
      # now get all the data together
      best_correlation_row <- data.frame(lane=c(lane), cluster=c(cluster), best_match_sample=c(best_match), second_match_sample=c(second_match), best_match_correlation=c(best_correlation), second_match_correlation=c(second_correlation))
      # add to the table
      if(is.null(best_correlation_per_lane)){
        best_correlation_per_lane <- best_correlation_row
      }
      else{
        best_correlation_per_lane <- rbind(best_correlation_per_lane, best_correlation_row)
      }
    }
  }
  return(best_correlation_per_lane)
}


add_donor_assignments <- function(seurat_object, best_correlations_table, cluster_column='soup_assignment', lane_column='lane'){
  # extract the metadata
  metadata <- seurat_object@meta.data
  # add an extra column for the combination of the lane and the cluster
  metadata[['lane_cluster']] <- paste(metadata[[lane_column]], metadata[[cluster_column]], sep = '_')
  best_correlations_table[['lane_cluster']] <- paste(best_correlations_table[['lane']], best_correlations_table[['cluster']], sep = '_')
  # add the assignments
  metadata[, c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')] <- best_correlations_table[match(metadata[['lane_cluster']], best_correlations_table[['lane_cluster']]), c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')]
  # now add it to the seurat object
  seurat_object <- AddMetaData(seurat_object, metadata[, c('best_match_sample', 'second_match_sample', 'best_match_correlation', 'second_match_correlation')])
  return(seurat_object)
}


check_nclusters_vs_ndonors <- function(metadata, lane_column='lane', cluster_column='soup_assignment', donor_column='best_match_sample'){
  # we will store the numbers per lane
  numbers_per_lane <- NULL
  # check each lane
  for(lane in unique(metadata[[lane_column]])){
    # subset to that lane
    metadata_lane <- metadata[metadata[[lane_column]] == lane, ]
    # check the unique number of clusters
    n_clusters <- length(unique(metadata_lane[[cluster_column]]))
    # chekc the number of unique donors
    n_donors <- length(unique(metadata_lane[[donor_column]]))
    # turn into a row
    numbers_row <- data.frame(lane = c(lane), n_clusters = c(n_clusters), n_donors = c(n_donors))
    # add to the rest of the table
    if(is.null(numbers_per_lane)){
      numbers_per_lane <- numbers_row
    }
    else{
      numbers_per_lane <- rbind(numbers_per_lane, numbers_row)
    }
  }
  # add a column denoting if the numbers match
  numbers_per_lane[['concordant']] <- numbers_per_lane[['n_clusters']] == numbers_per_lane[['n_donors']]
  return(numbers_per_lane)
}


add_inflammation_status <- function(seurat_object, sample_sheet, metadata_table, inflammation_column='inflammation_status', lane_column='lane', sample_column='donor_final'){
  # set the inflammation status empty
  seurat_object@meta.data[[inflammation_column]] <- NA
  # check each lane
  for(lane in lanes){
    # get the samples
    samples <- get_samples_in_lane(lane, sample_sheet, metadata_table)
    # get the inflammation by splitting on the comma, for that row
    inflammation <- strsplit(sample_sheet[sample_sheet[['lane']] == lane, 'inflammation_status'], ',')[[1]]
    # set a matching list
    inflammation_dict <- as.list(inflammation)
    names(inflammation_dict) <- samples
    # now get the donors that are in the in the Seurat metadata for that lane
    samples_in_metadata <- unique(seurat_object@meta.data[seurat_object@meta.data[[lane_column]] == lane, sample_column])
    # check each of the samples in the metadata
    for(sample in samples_in_metadata){
      # check if we have the status for this sample
      if(sample %in% names(inflammation_dict)){
        # get the inflammation status for that sample
        inflammation_status_sample <- inflammation_dict[[sample]]
        # add to the metadata
        seurat_object@meta.data[seurat_object@meta.data[[lane_column]] == lane &
                                  seurat_object@meta.data[[sample_column]] == sample, inflammation_column] <- inflammation_status_sample
      }
    }
  }
  return(seurat_object)
}


gsa_ids_to_sample_sheet <- function(sample_sheet, metadata_table, participants_sheet_column='samples', metadata_participant_column='Biopsy.storage.ID', metadata_gsa_column='GSA.sample.ID') {
  # create a new vector with th pasted GSA IDs
  gsa_rows <- rep(NA, times = nrow(sample_sheet))
  # now check each row
  for (i in 1:nrow(sample_sheet)) {
    # extract the biopsy IDs
    biopsy_ids_row_string <- sample_sheet[i, participants_sheet_column]
    # split them by the comma
    biopsy_ids_row <- strsplit(biopsy_ids_row_string, ',')[[1]]
    # now search all the GSA IDs for those
    gsa_ids_row <- metadata_table[match(biopsy_ids_row, metadata_table[[metadata_participant_column]]), metadata_gsa_column]
    # turn into a string
    gsa_ids_row_string <- paste(gsa_ids_row, collapse = ',')
    # put it into the vector
    gsa_rows[i] <- gsa_ids_row_string
  }
  # add the result to the table
  sample_sheet[['genoid']] <- gsa_rows
  return(sample_sheet)
}



####################
# Main Code        #
####################

# locations of the annotation files
souporcell_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/'
genotypes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/genotype_per_lane/'

lanes <- c('230105_lane1', '230105_lane2', '230105_lane3', '230105_lane4',
  '230105_lane5', '230105_lane6', '230105_lane7', '230105_lane8',
  '230112_lane1', '230112_lane2', '230112_lane3', '230112_lane4',
  '230112_lane5', '230112_lane6', '230112_lane7', '230112_lane8',
  '230120_lane1', '230120_lane2', '230120_lane3', '230120_lane4',
  '230120_lane5', '230120_lane6', '230120_lane7', '230120_lane8',
  '230127_lane1', '230127_lane2', '230127_lane3', '230127_lane4',
  '230127_lane5', '230127_lane6', '230127_lane7', '230127_lane8',
  '230202_lane1', '230202_lane2', '230202_lane3', '230202_lane4',
  '230202_lane5', '230202_lane6', '230202_lane7', '230202_lane8',
  '230209_lane1', '230209_lane2', '230209_lane3', '230209_lane4',
  '230209_lane5', '230209_lane6', '230209_lane7', '230209_lane8',
  '230216_lane1', '230216_lane2', '230216_lane3', '230216_lane4',
  '230216_lane5', '230216_lane6', '230216_lane7', '230216_lane8',
  '230223_lane1', '230223_lane2', '230223_lane3', '230223_lane4',
  '230223_lane5', '230223_lane6', '230223_lane7', '230223_lane8',
  '230302_lane1', '230302_lane2', '230302_lane3', '230302_lane4',
  '230302_lane5', '230302_lane6', '230302_lane7', '230302_lane8',
  '230316_lane1', '230316_lane2', '230316_lane3', '230316_lane4',
  '230316_lane5', '230316_lane6', '230316_lane7'
)

# get the correlations per lane
correlations_per_lane <- get_correlation_matrix_per_lane(souporcell_output_loc, genotypes_loc, lanes)
# get the best correlations
best_correlations <- get_best_correlations(correlations_per_lane)
