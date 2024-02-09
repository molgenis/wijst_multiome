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
library(cowplot)
# and get colors
library(RColorBrewer)


####################
# Functions        #
####################


#' calculate pearson correlations between genotypes in reference and cluster genotypes
#' 
#' @param df dataframe to add the correlations to
#' @param ref_df the dataframe that has the reference genotypes 
#' @param clust_df the dataframe that has the cluster genotypes
#' @returns the original dataframe, with the correlations added
#' 
pearson_correlation <- function(df, ref_df, clust_df){
  for (col in colnames(df)){
    for (row in rownames(df)){
      df[row,col] <- cor(as.numeric(pull(ref_df, col)), as.numeric(pull(clust_df, row)), method = "pearson", use = "complete.obs")
    }
  }
  return(df)
}

# method taken from https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/blob/master/Demultiplexing/includes/Snakefile_souporcell.smk
correlate_genotypes <- function(ref_geno_loc=NULL, cluster_geno_loc=NULL, ref_geno_vcfr=NULL, cluster_geno_vcfr=NULL){
  ref_geno <- NULL
  if (!is.null(ref_geno_vcfr)) {
    ref_geno <- ref_geno_vcfr
  }
  else if (!is.null(ref_geno_loc)) {
    ref_geno <- read.vcfR(ref_geno_loc)
  }
  else {
    stop('supply either reference VCF location, or vcfR object')
  }
  cluster_geno <- NULL
  if (!is.null(cluster_geno_vcfr)) {
    cluster_geno <- cluster_geno_vcfr
  }
  else if (!is.null(cluster_geno_loc)) {
    cluster_geno <- read.vcfR(cluster_geno_loc)
  }
  else {
    stop('supply either cluster VCF location, or vcfR object')
  }
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

#' get a correlate the genotypes of two souporcell outputs
#' 
#' @param cluster_geno_loc1 location of the first VCF
#' @param cluster_geno_loc2 location of the second VCF
#' @returns dataframe with the correlations between genotypes in the first and second genotype file
#' 
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

#' get a correlate the genotypes of souporcell and priori genotypes
#' 
#' @param souporcell_output_loc location of the souporcell output
#' @param genotypes_loc locations of the VCFs of genotypes per lane
#' @param lanes lanes to look at
#' @param genotype_prepend what the prior genotypes start with in the filename
#' @param genotype_append what the prior genotypes end with in the filename
#' @returns dataframe with the correlations between genotypes of the a priori genotypes and souporcell genotypes, in a list with the lanes as keys
#' 
get_correlation_matrix_per_lane <- function(souporcell_output_loc, genotypes_loc, lanes, genotype_prepend='mo_', genotype_append='_maf005.vcf.gz'){
  # we will save the correlations per lane
  correlations_per_lane <- list()
  # check each lane
  for(lane in lanes){
    print(paste('correlating lane', lane))
    # paste together the souporcell cluster loc
    cluster_loc <- paste(souporcell_output_loc, lane, '/cluster_genotypes.vcf', sep = '')
    # paste together the genotype loc
    geno_loc <- paste(genotypes_loc, genotype_prepend, lane, genotype_append, sep = '')
    if(file.exists(cluster_loc) & file.exists(geno_loc)) {
      # get the result
      correlations_lane <- correlate_genotypes(geno_loc, cluster_loc)
      # store the result in the list
      correlations_per_lane[[lane]] <- correlations_lane
    }
    else{
      warning(paste('skipped lane', lane, 'due to no input for the clusters or reference VCF'))
    }
  }
  return(correlations_per_lane)
}

#' get a correlate the genotypes of souporcell and priori genotypes
#' 
#' @param correlations_per_lane the correlations per lane, which is a list with the lanes as keys, and the correlation matrices per lane as values
#' @returns dataframe with the best correlating cluster and genotype and their correlation, as well as the second best matches
#' 
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

#' add donor assignments based on souporcell genotypes to the Seurat object
#' 
#' @param seurat_object the Seurat object to add the assignments to
#' @param best_correlations_table the table that has the best matching sample per souporcell clusters
#' @param cluster_column the column in the Seurat metadata that has the souporcell cluster
#' @param lane_column the column in the Seurat metadata that has the lane
#' @returns Seurat object with the assignment and correlation
#' 
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

#' check if the number of unique clusters matches the number of unique donors for each lane
#' 
#' @param metadata the Seurat object with the metadata
#' @param lane_column the column in the Seurat metadata that has the lane
#' @param cluster_column the column in the Seurat metadata that has the souporcell cluster
#' @param donor_column the column in the Seurat metadata that has the donor
#' @returns Seurat object with the assignment and correlation
#' 
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


#' get a vector of as distinct possible colours
#' 
#' @param number_of_colours how many colours to return
#' @param use_sampling whether or not to randomly extract the colours instead of grabbing the first n colours
#' @param color_indices (optional, not used by default) if specific colours are needed, supply the indices of the colours here. Use 'get_available_colours_grid' to get the colours and their indices
#' @returns a vector of colours
#' 
sample_many_colours <- function(number_of_colours, use_sampling=F, color_indices=NULL) {
  # get all colours from the 'quality' palettes
  quality_colour_palettes <- brewer.pal.info[brewer.pal.info[['category']] == 'qual', ]
  # save each palette
  colours_per_palette <- list()
  # apply over each palette
  for (i in 1:nrow(quality_colour_palettes)) {
    # get the name of the palette
    palette_name <- rownames(quality_colour_palettes)[i]
    # get the number of colours in the palette
    palette_max_colours <- quality_colour_palettes[i, 'maxcolors']
    # use brewer.pal to get all colours
    colours_palette <- brewer.pal(palette_max_colours, palette_name)
    # put result in the list
    colours_per_palette[[palette_name]] <- colours_palette
  }
  # merge all palettes
  all_colours <- do.call('c', colours_per_palette)
  # randomly get colours from the palette
  max_possible_colours <- length(all_colours)
  if (is.null(number_of_colours)) {
    message('no number of colors supplied, assuming color indices have been')
  }
  else if (number_of_colours > max_possible_colours) {
    message(paste('requesting more colours than is possible: ', as.character(number_of_colours), ' vs ', max_possible_colours, ', returning max possible', sep = ''))
    number_of_colours <- max_possible_colours
  }
  colours_to_return <- NULL
  # specific colours we like (the indices)
  if (!is.null(color_indices)) {
    colours_to_return <- all_colours[color_indices]
  }
  # or use sampling
  else if (use_sampling) {
    colours_to_return <- sample(all_colours, number_of_colours)
  }
  # or the first x colours
  else {
    colours_to_return <- all_colours[1 : number_of_colours]
  }
  return(colours_to_return)
}

#' get a vector of as distinct possible colours, but with more possibilities ()
#' 
#' @param number_of_colours how many colours to return
#' @param use_sampling whether or not to randomly extract the colours instead of grabbing the first n colours
#' @param color_indices (optional, not used by default) if specific colours are needed, supply the indices of the colours here. Use 'get_available_colours_grid' to get the colours and their indices
#' @returns a vector of colours
#' 
sample_tons_of_colors <- function(number_of_colours, use_sampling=F, color_indices=NULL) {
  # get colours available to device
  all_colours <- grDevices::colors()
  # remove gray
  all_colours <- all_colours[grep('gr(a|e)y', all_colours, invert = T)]
  # check how many are possible
  max_possible_colours <- length(all_colours)
  if (is.null(number_of_colours)) {
    message('no number of colors supplied, assuming color indices have been')
  }
  else if (number_of_colours > max_possible_colours) {
    message(paste('requesting more colours than is possible: ', as.character(number_of_colours), ' vs ', max_possible_colours, ', returning max possible', sep = ''))
    number_of_colours <- max_possible_colours
  }
  colours_to_return <- NULL
  # specific colours we like (the indices)
  if (!is.null(color_indices)) {
    colours_to_return <- all_colours[color_indices]
  }
  # or use sampling
  else if (use_sampling) {
    colours_to_return <- sample(all_colours, number_of_colours)
  }
  # or the first x colours
  else {
    colours_to_return <- all_colours[1 : number_of_colours]
  }
  return(colours_to_return)
}


#' get a grid showing the available colours and their indices
#' 
#' @param many use the 'many' method to get the colours
#' @param tons use the 'tons' method to get the colours
#' @returns a ggplot grid showing the available colours and their indices
#' 
get_available_colours_grid <- function(many=T, tons=F) {
  colours_possible <- NULL
  # get from the many method
  if (many) {
    # ask for unreasonable amount
    colours_possible <- sample_many_colours(1000)
  }
  else if(tons) {
    colours_possible <- sample_tons_of_colors(1000)
  }
  # get how many colours we actually have
  available_colours <- length(colours_possible)
  # we need to put that into a square grid, so we need to get the square root, to know how many rows and columns
  nrow_and_ncol <- sqrt(available_colours)
  # and we need to round that up of course
  nrow_and_ncol <- ceiling(nrow_and_ncol)
  # so we'll have a total number of blocks
  total_cells <- nrow_and_ncol * nrow_and_ncol
  # let's see how many colours we are off from that number of cells
  cells_no_colour_number <- total_cells - available_colours
  # we will just add white for those
  cells_no_colour <- rep('white', times = cells_no_colour_number)
  # add that to the colours we have
  colours_possible <- c(colours_possible, cells_no_colour)
  # create each combination of x and y
  indices_grid <- expand.grid(as.character(1 : nrow_and_ncol), as.character(1 : nrow_and_ncol))
  # add the index and colour name
  indices_grid[['index_colour']] <- paste(c(1:total_cells), colours_possible, sep = '\n')
  # make mapping of colours
  colours_to_use <- as.list(colours_possible)
  names(colours_to_use) <- indices_grid[['index_colour']]
  # now plot
  p <- ggplot(data = indices_grid, mapping = aes(x = Var1, y = Var2, fill = index_colour)) + 
    geom_tile() + 
    geom_text(aes(label=index_colour)) + 
    scale_fill_manual(values = colours_to_use) + 
    theme(legend.position = 'none')
  return(p)
}

#' add the rank of each participant of every lane, based on their correlations
#' 
#' @param correlation_table the table containing the correlations
#' @param correlation_column which column contains correlation number
#' @param lane_column which column contains the lane of the cell
#' @returns the input dataframe, with a new columm containing the rank of the participant in each lane
#' 
samples_to_rankings <- function(correlation_table, sample_column='best_match_sample', correlation_column='best_match_correlation', lane_column='lane') {
  # we'll first do this per lane
  ranked_correlations_per_lane <- list()
  # going through the lanes
  for (lane in unique(correlation_table[[lane_column]])) {
    # subset to that lane
    correlations_lane <- correlation_table[correlation_table[[lane_column]] == lane, ]
    # order this per correlation
    correlations_lane <- correlations_lane[order(correlations_lane[[correlation_column]]), ]
    # change to rank
    correlations_lane[[sample_column]] <- as.character(1:nrow(correlations_lane))
    # put in the list
    ranked_correlations_per_lane[[lane]] <- correlations_lane
  }
  # merge them all
  ranked_correlations <- do.call('rbind', ranked_correlations_per_lane)
  return(ranked_correlations)
}

#' plot the demuxlet assignments
#' 
#' @param correlations_table table with correlations and samples
#' @param sample_column the sample names
#' @param correlation_column column containing the correlation to plot
#' @param lane_column which column contains the lane of the cell (optional, 'lane' by default)
#' @param assignment_to_ranking use instead of the sample assignment, the rank of the sample (optional, T by default)
#' @param ylim supply the y limits as vector (optional)
#' @param pointless remove ticks on the bottom of plot (optional, default F)
#' @param legendless remove the legend (optional, default F)
#' @param paper_style add extra whitespace to plot (optional, default T)
#' @param angle_labels angle the x axis labels 90 degrees (optional, default T)
#' @param to_fractions plot fraction of cells belonging to sample in each lane instead of number of cells (optional, default F)
#' @param use_distinct_colours use the distinct colour palette instead of ggplot defaults (optional, default T)
#' @param use_sampling randomly sample colours from the distinct colour palette (optional, default F)
#' @param color_indices supply specific indices relating to colours to use from the distinct colour palette. Use 'get_available_colours_grid' to get the colours and their indices (optional, unused by default)
#' @returns plot with correlations in each lane
#' 
plot_correlations_per_lane <- function(correlations_table, sample_column='best_match_sample', correlation_column='best_match_correlation', lane_column='lane', assignment_to_ranking=T, ylim=NULL, pointless=F, legendless=F, paper_style=T, angle_labels=T, to_fractions=F, use_distinct_colours=T, use_sampling=F, color_indices=NULL) {
  # convert to ranks if requested
  if (assignment_to_ranking) {
    correlations_table <- samples_to_rankings(correlations_table, sample_column = sample_column, correlation_column = correlation_column, lane_column = lane_column)
  }
  # order the samples alphabetically
  correlations_table[[sample_column]] <- factor(correlations_table[[sample_column]], levels = sort(unique(correlations_table[[sample_column]])))
  # initialize the plot
  p <- ggplot(data = NULL, mapping = aes(x = correlations_table[[lane_column]], y = correlations_table[[correlation_column]], fill = correlations_table[[sample_column]])) + geom_bar(position='stack', stat='identity') + xlab(lane_column) + ylab('correlations')
  # use distinct colours if requested
  if (use_distinct_colours) {
    # get the unique possible assigments
    possible_assignments <- unique(correlations_table[[sample_column]])
    # get an equal amount of colours
    possible_colours <- NULL
    if(length(possible_assignments) > 74) {
      possible_colours <- sample_tons_of_colors(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    else {
      possible_colours <- sample_many_colours(length(possible_assignments), use_sampling = use_sampling, color_indices = color_indices)
    }
    # put into a list
    colour_mapping <- as.list(possible_colours)
    names(colour_mapping) <- possible_assignments
    # add to plot
    p <- p + scale_fill_manual(values = colour_mapping)
  }
  # add legend based on rank or assignment
  if (assignment_to_ranking) {
    p <- p + guides(fill=guide_legend(title='sample rank'))
  }
  else{
    p <- p + guides(fill=guide_legend(title='sample'))
  }
  # with some options
  if(!is.null(ylim)){
    p <- p + ylim(ylim)
  }
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p)
}

#' check which participants are missing in the best correlations, that should be present based on the sample sheets
#' 
#' @param best_correlations_table the table that has the best matching sample per souporcell clusters
#' @param participant_per_lane_loc location of the folder that has the lists of participants per lane
#' @param participant_per_lane_loc_prepend what the prior participants per lane start with in the filename
#' @param participant_per_lane_loc_apppend what the prior participants per lane start with in the filename
#' @returns Seurat object with the assignment and correlation
#' 
get_missing_participants_per_lane <- function(best_correlations, participant_per_lane_loc, participant_per_lane_loc_prepend='', participant_per_lane_loc_apppend='.txt', lane_column='lane', donor_column='best_match_sample', correlation_column='best_match_correlation') {
  # get the lanes from the correlation table
  lanes <- unique(best_correlations[[lane_column]])
  # and we need to know how many lanes we have
  nlanes <- length(lanes)
  # create a dataframe to store the participants that are
  participants_missing_per_lane = data.frame(lane = rep(NA, times = nlanes), missing = rep(NA, times = nlanes), present_instead = rep(NA, times = nlanes), present_correlations = rep(NA, times = nlanes))
  # now check each lane, by index
  for (i in 1 : nlanes) {
    # extract lane
    lane = lanes[i]
    # get the participants
    participants_matched <- unique(best_correlations[best_correlations[[lane_column]] == lane, donor_column])
    # get the location of the annotation that has the partipants per lane
    parts_per_lane_loc_lane <- paste(participant_per_lane_loc, participant_per_lane_loc_prepend, '/', lane, participant_per_lane_loc_apppend, sep = '')
    # get the participants
    parts_should_be_in_lane <- read.table(parts_per_lane_loc_lane, header = F)$V1
    # get what is missing
    parts_missing <- setdiff(parts_should_be_in_lane, participants_matched)
    # order to make more readable
    parts_missing <- parts_missing[order(parts_missing)]
    # get which were added instead
    present_instead <- setdiff(participants_matched, parts_should_be_in_lane)
    # order to make more readable
    present_instead <- present_instead[order(present_instead)]
    # and their correlations
    present_correlations <- best_correlations[best_correlations[[lane_column]] == lane, ][match(present_instead, best_correlations[best_correlations[[lane_column]] == lane, donor_column]), correlation_column]
    present_correlations <- round(present_correlations, digits = 2)
    # add to the dataframe
    participants_missing_per_lane[i, 'lane'] <- lane
    participants_missing_per_lane[i, 'missing'] <- paste(parts_missing, collapse = ',')
    participants_missing_per_lane[i, 'present_instead'] <- paste(present_instead, collapse = ',')
    participants_missing_per_lane[i, 'present_correlations'] <- paste(present_correlations, collapse = ',')
  }
  return(participants_missing_per_lane)
}

#' check which participants are missing in the best correlations, that should be present based on the sample sheets
#' 
#' @param best_correlations_table the table that has the best matching sample per souporcell clusters
#' @param participant_per_lane_loc location of the folder that has the lists of participants per lane
#' @param participant_per_lane_loc_prepend what the prior participants per lane start with in the filename
#' @param participant_per_lane_loc_apppend what the prior participants per lane ends with in the filename
#' @returns Seurat object with the assignment and correlation
#' 
create_assignment_per_barcode <- function(souporcell_output_loc, best_assignments, lanes) {
  # we will initially store per lane
  best_match_per_barcode_lane <- list()
  # we will check each lane
  for (lane in lanes) {
    # paste together the souporcell clusters file
    barcode_clusters_loc <- paste(souporcell_output_loc, '/', lane, '/clusters.tsv', sep = '')
    # check if the file exists
    if (file.exists(barcode_clusters_loc)) {
      # read the file
      barcode_clusters <- read.table(barcode_clusters_loc, header = T, sep = '\t')
      # add the lane as explicit column
      barcode_clusters_lane <- rep(lane, times = nrow(barcode_clusters))
      # add the bare barcode
      barcodes_short_cluster <- gsub('(-\\d+)', '', barcode_clusters[['barcode']])
      barcodes_lane_cluster <- paste(barcodes_short_cluster, rep(lane, times = length(barcodes_short_cluster)), sep = '_')
      # subset the best assignments to this lane
      barcodes_assignments <- best_assignments[best_assignments[['lane']] == lane, ]
      # check if we have those assignments
      if (nrow(barcodes_assignments) > 0) {
        # create the dataframe
        assignment_table_full <- data.frame(
          lane = barcode_clusters_lane,
          barcode_lane = barcodes_lane_cluster,
          barcode = barcodes_short_cluster,
          barcode_original = barcode_clusters[['barcode']],
          status = barcode_clusters[['status']],
          cluster = barcode_clusters[['assignment']]
        )
        # make the cluster a character vector, because in souporcell the cluster can be formatted as '0/7' for doublets
        barcodes_assignments[['cluster']] <- as.character(barcodes_assignments[['cluster']])
        # join this onto the assignments that we have
        assignment_table_full <- merge(assignment_table_full, barcodes_assignments, by = c('lane', 'cluster'), all = T)
        # now add back the info from the original souporcell output, that we didn't already have
        assignment_table_full <- cbind(assignment_table_full, # add to original table
                                       barcode_clusters[match(assignment_table_full[['barcode_original']], barcode_clusters[['barcode']]), # match to original by barcode
                                                        setdiff(colnames(barcode_clusters), c('barcode', 'status', 'assignment'))]) # get all columns, except the ones we specify
        # add to list
        best_match_per_barcode_lane[[lane]] <- assignment_table_full
      }
      else{
        warning(paste('skipping', lane, 'because correlations are missing'))
      }
    }
    else{
      warning(paste('skipping', lane, 'because file is missing at', barcode_clusters_loc))
    }
    
  }
  # merge everything
  full_assignments <- do.call('rbind', best_match_per_barcode_lane)
  # order by lane, then cluster
  full_assignments <- full_assignments[order(full_assignments[['lane']], full_assignments[['cluster']]), ]
  return(full_assignments)
}


wide_to_high_ggplot <- function(wide_table, variable_col_name='correlation', new_col_name='sample', new_row_name='cluster'){
  # init new table
  table_high <- NULL
  for(col in colnames(wide_table)){
    # get the variables in the column
    variables <- wide_table[[col]]
    # turn into dataframe
    table_high_rows <- data.frame(x=rep(col, times=length(variables)), y=rownames(wide_table), z=variables)
    # set column names
    colnames(table_high_rows) <- c(new_col_name, new_row_name, variable_col_name)
    # add to the rest of the table
    if(is.null(table_high)){
      table_high <- table_high_rows
    }
    else{
      table_high <- rbind(table_high, table_high_rows)
    }
  }
  return(table_high)
}



correlations_to_tile <- function(correlation_lane, plot_text=F, pointless=F, legendless=F, ylim=NULL, paper_style=T, angle_labels=T) {
  # first convert to a long table
  correlations_long <- wide_to_high_ggplot(correlation_lane)
  # round the value
  correlations_long[['value_rounded']] <- as.character(round(correlations_long[['correlation']], digits = 2))
  # create the ploc
  p <- ggplot(NULL, aes(x = correlations_long[['sample']], y = correlations_long[['cluster']], fill = correlations_long[['correlation']])) +
    geom_tile() +
    coord_fixed() + 
    scale_fill_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = .5, limits = c(0,1)) +
    labs(fill = 'correlation') +
    xlab('sample') +
    ylab('cluster')
  # plot the text if requested
  if (plot_text) {
    p <- p + geom_text(aes(label = correlations_long[['value_rounded']]), color = "black", size = text_size)
  }
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(),
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p)
}


plot_correlation_densities <- function(correlations, correlation_correlation_column='best_match_correlation', correlation_group_column='vs') {
  # get the fill categories
  fill_categories <- unique(correlations[[correlation_group_column]])
  # get the colors
  cols <- as.list(sample_many_colours(length(fill_categories)))
  # set the original categories as names
  names(cols) <- fill_categories
  # make the plot
  p <- ggplot(NULL, aes(correlations[[correlation_correlation_column]], fill = correlations[[correlation_group_column]])) + 
    geom_density(alpha = 0.2) + 
    scale_fill_manual(name = correlation_group_column, values = cols) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
    xlab(correlation_correlation_column) + 
    ylab('Density')
  return(p)
}


create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL, legendless=F){
  # init the table
  confusion_table <- NULL
  # check each truth
  for(truth in unique(assignment_table[[truth_column]])){
    # get these truths
    truth_rows <- assignment_table[assignment_table[[truth_column]] == truth, ]
    # check now many have this truth
    this_truth_number <- nrow(truth_rows)
    # check what was predicted for these truths
    #for(prediction in unique(truth_rows[[prediction_column]])){
    for(prediction in unique(assignment_table[[prediction_column]])){
      # check the number of this prediction
      this_prediction_number <- nrow(truth_rows[truth_rows[[prediction_column]] == prediction, ])
      # init variable
      fraction <- NULL
      # we can only calculate a fraction if the result is not zero
      if(this_prediction_number > 0){
        # calculate the fraction
        fraction <- this_prediction_number / this_truth_number
      }
      # otherwise we just set it to zero
      else{
        fraction <- 0
      }
      # turn into row
      this_row <- data.frame(truth=c(truth), prediction=c(prediction), freq=c(fraction), stringsAsFactors = F)
      # add this entry to the dataframe
      if(is.null(confusion_table)){
        confusion_table <- this_row
      }
      else{
        confusion_table <- rbind(confusion_table, this_row)
      }
    }
  }
  # round the frequency off to a sensible cutoff
  confusion_table$freq <- round(confusion_table$freq, digits=2)
  # turn into plot
  p <- ggplot(data=confusion_table, aes(x=truth, y=prediction, fill=freq)) + geom_tile() + scale_fill_gradient(low='red', high='blue') + geom_text(aes(label=freq))
  # some options
  if(!is.null(truth_column_label)){
    p <- p + xlab(truth_column_label)
  }
  if(!is.null(prediction_column_label)){
    p <- p + ylab(prediction_column_label)
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}


####################
# Main Code        #
####################

# locations of the annotation files
souporcell_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/gex/barcode_filtered/'
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
           '230316_lane5', '230316_lane6', '230316_lane7', '230316_lane8'
)

# get the correlations per lane
correlations_per_lane <- get_correlation_matrix_per_lane(souporcell_output_loc, genotypes_loc, lanes)
# Warning messages:
#   1: In get_correlation_matrix_per_lane(souporcell_output_loc, genotypes_loc,  :
#   skipped lane 230223_lane7 due to no input for the clusters or reference VCF
#   2: In get_correlation_matrix_per_lane(souporcell_output_loc, genotypes_loc,  :
#   skipped lane 230223_lane8 due to no input for the clusters or reference VCF
# get the best correlations
best_correlations <- get_best_correlations(correlations_per_lane)
# save the results
saveRDS(correlations_per_lane, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected.rds')
write.table(best_correlations, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_best_assignments.tsv', sep = '\t', row.names = F, col.names = T, quote = F)

# plot what the best correlations look like
plot_correlations_per_lane(best_correlations)

# add correlation data
correlation_mapping_per_barcode <- create_assignment_per_barcode(souporcell_output_loc, best_correlations, lanes = lanes)
write.table(correlation_mapping_per_barcode, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched.tsv', row.names = F,col.names = T, quote = F, sep = '\t')

# check how many samples per lane are assigned (should be 8 unique ones every time)
samples_per_lane <- unique(best_correlations[, c('lane', 'best_match_sample')])
nsample_per_lane <- data.frame(table(samples_per_lane[['lane']]))

# get which samples are missing
samples_missing_per_lane <- get_missing_participants_per_lane(best_correlations, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/participant_per_lane/')
write.table(samples_missing_per_lane[samples_missing_per_lane$missing != '', ], '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_missings.tsv', row.names = F,col.names = T, quote = F, sep = '\t')

# check some correlations
ref_geno_all_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/imputed_hg38_all_anc_mmaf005_chrprepend.vcf.gz'
ref_geno_all <- read.vcfR(ref_geno_all_loc)
# we'll save for each lane
vs_all_per_lane <- list()
# check lanes with missing participants
for (lane in lanes) {
  print(lane)
  # calculate correlations
  correlations_vs_all <- correlate_genotypes(
    ref_geno_vcfr = ref_geno_all, 
    cluster_geno_loc = paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/gex/barcode_filtered/', lane, '/cluster_genotypes.vcf', sep = ''),
    ref_geno_loc = NULL,
    cluster_geno_vcfr = NULL
  )
  # add to list
  vs_all_per_lane[[lane]] <- correlations_vs_all
}

# save the result
saveRDS(vs_all_per_lane, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_vsall.rds')
best_correlations_vs_all <- get_best_correlations(vs_all_per_lane)
write.table(best_correlations_vs_all, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_vsall_best_assignments.tsv', sep = '\t', row.names = F, col.names = T, quote = F)

# plot what the best correlations look like
plot_correlations_per_lane(best_correlations_vs_all)

# add correlation data
correlation_mapping_per_barcode_all <- create_assignment_per_barcode(souporcell_output_loc, best_correlations_vs_all, lanes = lanes)
write.table(correlation_mapping_per_barcode_all, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched_vs_all.tsv', row.names = F,col.names = T, quote = F, sep = '\t')

# get which samples are missing
samples_missing_per_lane_all <- get_missing_participants_per_lane(best_correlations_vs_all, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/participant_per_lane/')
write.table(samples_missing_per_lane_all, '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_corrected_sample_matched_vs_all_missings.tsv', sep = '\t', row.names = F, col.names = T)

# plot the per-lane and all-lane correlations per lane
plot_grid(
  plot_correlations_per_lane(best_correlations, legendless = T) + ggtitle('only sheet samples'),
  plot_correlations_per_lane(best_correlations_vs_all, legendless = T) + ggtitle('all samples'),
  get_legend(plot_correlations_per_lane(best_correlations, legendless = F)),
  rel_widths = c(3,3,1),
  nrow = 1
)

# combine the correlations
best_correlations[['vs']] <- 'lane samples'
best_correlations_vs_all[['vs']] <- 'all samples'
best_correlations_both <- rbind(best_correlations, best_correlations_vs_all)
# plot
plot_correlation_densities(best_correlations_both)

# plot the number of lanes with issues
nlanes_problem <- length(unique(best_correlations[best_correlations$best_match_correlation < 0.7, 'lane']))
nlanes_problem_vs_all <- length(unique(best_correlations_vs_all[best_correlations_vs_all$best_match_correlation < 0.7, 'lane']))
p_nlanes_problems <- ggplot(data = data.frame(lanes = c('total', 'lane_samples', 'all_samples'),  number = c(length(unique(best_correlations$lane)), nlanes_problem, nlanes_problem_vs_all)), mapping = aes(x = lanes, y = number, fill = lanes)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('darkred', 'darkblue', '#FFEA00')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# plot the number of samples with issues
ndonors_problem <- nrow(best_correlations[best_correlations$best_match_correlation < 0.7, ])
ndonors_problem_vs_all <- nrow(best_correlations_vs_all[best_correlations_vs_all$best_match_correlation < 0.7, ])
p_ndonors_problems <- ggplot(data = data.frame(samples = c('total', 'lane_samples', 'all_samples', 'all_samples_minus_missing'),  number = c(length(unique(best_correlations$best_match_sample)), ndonors_problem, ndonors_problem_vs_all, ndonors_problem_vs_all - 4)), mapping = aes(x = samples, y = number, fill = samples)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('darkred', 'purple', 'darkblue', '#FFEA00')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# plot the missing correlations
plot_correlation_densities(best_correlations_both[best_correlations_both[['best_match_correlation']] < 0.7, ])

# plot the number of lanes with issues
nlanes_problem_06 <- length(unique(best_correlations[best_correlations$best_match_correlation < 0.6, 'lane']))
nlanes_problem_vs_all_06 <- length(unique(best_correlations_vs_all[best_correlations_vs_all$best_match_correlation < 0.65, 'lane']))
p_nlanes_problems_06 <- ggplot(data = data.frame(lanes = c('total', 'lane_samples', 'all_samples'),  number = c(length(unique(best_correlations$lane)), nlanes_problem_06, nlanes_problem_vs_all_06)), mapping = aes(x = lanes, y = number, fill = lanes)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('darkred', 'darkblue', '#FFEA00')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  
# plot the number of samples with issues
ndonors_problem_06 <- nrow(best_correlations[best_correlations$best_match_correlation < 0.6, ])
ndonors_problem_vs_all_06 <- nrow(best_correlations_vs_all[best_correlations_vs_all$best_match_correlation < 0.6, ])
p_ndonors_problems_06 <- ggplot(data = data.frame(samples = c('total', 'lane_samples', 'all_samples', 'all_samples_minus_missing'),  number = c(length(unique(best_correlations$best_match_sample)), ndonors_problem_06, ndonors_problem_vs_all_06, ndonors_problem_vs_all_06 - 4)), mapping = aes(x = samples, y = number, fill = samples)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('darkred', 'purple', 'darkblue', '#FFEA00')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))

# plot 0.6 vs 0.7
plot_grid(
  p_ndonors_problems + ggtitle('0.7 cutoff') + theme(legend.position = 'none') + theme(axis.text.x=element_blank(), axis.ticks = element_blank()),
  p_ndonors_problems_06 + ggtitle('0.6 cutoff') + theme(legend.position = 'none') + theme(axis.text.x=element_blank(), axis.ticks = element_blank()),
  get_legend(p_ndonors_problems_06),
  nrow = 1
)
plot_grid(
  p_nlanes_problems + ggtitle('0.7 cutoff') + theme(legend.position = 'none') + theme(axis.text.x=element_blank(), axis.ticks = element_blank()),
  p_nlanes_problems_06 + ggtitle('0.6 cutoff') + theme(legend.position = 'none') + theme(axis.text.x=element_blank(), axis.ticks = element_blank()),
  get_legend(p_nlanes_problems),
  nrow = 1
)

# load the previous barcodes mapping
correlation_mapping_per_barcode_all_old <- read.table('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_sample_matched.tsv', sep = '\t', header = T)
# merge them based on what we care about
correlation_mappings_old_vs_new <- merge(
  x = correlation_mapping_per_barcode_all_old[, c('lane', 'cluster', 'barcode_lane')], 
  y = correlation_mapping_per_barcode_all[, c('lane', 'cluster', 'barcode_lane')], 
  by.x = 'barcode_lane', by.y = 'barcode_lane'
)
# replace all the / with 'multiplet'
correlation_mappings_old_vs_new[['cluster.x']] <- gsub('\\d+/\\d+', 'multiplet', correlation_mappings_old_vs_new[['cluster.x']])
correlation_mappings_old_vs_new[['cluster.y']] <- gsub('\\d+/\\d+', 'multiplet', correlation_mappings_old_vs_new[['cluster.y']])

plot_grid(
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[1], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[1]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[2], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[2]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[3], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[3]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[4], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[4]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[5], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[5]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[6], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[6]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[7], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[7]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[8], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[8]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[9], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[9]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[10], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[10]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[11], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[11]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[12], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[12]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[13], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[13]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[14], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[14]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[15], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[15]),
  create_confusion_matrix(correlation_mappings_old_vs_new[correlation_mappings_old_vs_new$lane.x == lanes[16], ], truth_column = 'cluster.x', prediction_column = 'cluster.y', truth_column_label = 'uncorrected', prediction_column_label = 'corrected', legendless = T) + ggtitle(lanes[16]),
  nrow = 4,
  ncol = 4
)
