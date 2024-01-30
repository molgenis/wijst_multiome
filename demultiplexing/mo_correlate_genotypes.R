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
# and get colors
library(RColorBrewer)


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

get_missing_participants_per_lane <- function(best_correlations, participant_per_lane_loc, participant_per_lane_loc_prepend='', participant_per_lane_loc_apppend='.txt', lane_column='lane', donor_column='best_match_sample') {
  # get the lanes from the correlation table
  lanes <- unique(best_correlations[[lane_column]])
  # and we need to know how many lanes we have
  nlanes <- length(lanes)
  # create a dataframe to store the participants that are
  participants_missing_per_lane = data.frame(lane = rep(NA, times = nlanes), missing = rep(NA, times = nlanes))
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
    # add to the dataframe
    participants_missing_per_lane[i, 'lane'] <- lane
    participants_missing_per_lane[i, 'missing'] <- paste(parts_missing, collapse = ',')
  }
  return(participants_missing_per_lane)
}


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

####################
# Main Code        #
####################

# locations of the annotation files
souporcell_output_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/'
genotypes_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/genotype/genotype_per_lane/'

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
# Warning messages:
#   1: In get_correlation_matrix_per_lane(souporcell_output_loc, genotypes_loc,  :
#   skipped lane 230223_lane7 due to no input for the clusters or reference VCF
#   2: In get_correlation_matrix_per_lane(souporcell_output_loc, genotypes_loc,  :
#   skipped lane 230223_lane8 due to no input for the clusters or reference VCF
# get the best correlations
best_correlations <- get_best_correlations(correlations_per_lane)
# save the results
saveRDS(correlations_per_lane, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected.rds')
write.table(best_correlations, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_best_assignments.tsv', sep = '\t', row.names = F, col.names = T, quote = F)

# plot what the best correlations look like
plot_correlations_per_lane(best_correlations)

# check how many samples per lane are assigned (should be 8 unique ones every time)
samples_per_lane <- unique(best_correlations[, c('lane', 'best_match_sample')])
nsample_per_lane <- data.frame(table(samples_per_lane[['lane']]))

# get which samples are missing
samples_missing_per_lane <- get_missing_participants_per_lane(best_correlations, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/participant_per_lane/')
write.table(samples_missing_per_lane[samples_missing_per_lane$missing != '', ], '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_missings.tsv', row.names = F,col.names = T, quote = F, sep = '\t')

# check some correlations
ref_geno_all_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/genotype/imputed_hg38_all_anc_mmaf005_chrprepend.vcf.gz'
ref_geno_all <- read.vcfR(ref_geno_all_loc)
# we'll save for each lane
vs_all_per_lane <- list()
# check lanes with missing participants
for (lane in lanes) {
  print(lane)
  # calculate correlations
  correlations_vs_all <- correlate_genotypes(
    ref_geno_vcfr = ref_geno_all, 
    cluster_geno_loc = paste('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/', lane, '/cluster_genotypes.vcf', sep = ''),
    ref_geno_loc = NULL,
    cluster_geno_vcfr = NULL
  )
  # add to list
  vs_all_per_lane[[lane]] <- correlations_vs_all
}

# save the result
saveRDS(vs_all_per_lane, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_vsall.rds')
best_correlations_vs_all <- get_best_correlations(vs_all_per_lane)
write.table(best_correlations_vs_all, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_vsall_best_assignments.tsv', sep = '\t', row.names = F, col.names = T, quote = F)


# add correlation data
correlation_mapping_per_barcode <- create_assignment_per_barcode(souporcell_output_loc, best_correlations, lanes = lanes)
write.table(correlation_mapping_per_barcode, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_souporcell_gex_uncorrected_sample_matched.tsv', row.names = F,col.names = T, quote = F, sep = '\t')

# try again with the gex data
correlations_per_lane <- get_correlation_matrix_per_lane('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/gex/cellbent/', genotypes_loc, lanes, genotype_prepend = '', genotype_append = '')

# try with filtered input
correlations_per_lane_barcodefilter <- get_correlation_matrix_per_lane('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/barcode_filtered/', genotypes_loc, lanes[1:8])
best_correlations_barcodefilter <- get_best_correlations(correlations_per_lane_barcodefilter)
# check comparison
best_correlations_compared <- merge(best_correlations, best_correlations_barcodefilter, by = c('lane', 'cluster'))

# check imputed vs unimputed
nonimputed_b38_loc <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/genotype/GSA2023_1044_025_V3/unimputed/GSA2022_1044_025_V3_rsid_noindels_nodots_onlychr_b38.vcf.gz'
# read the files
nonimputed_b38 <- read.vcfR(nonimputed_b38_loc)
# calculate correlations
correlations_imp_nonimp <- correlate_genotypes(
  ref_geno_vcfr = ref_geno_all, 
  cluster_geno_loc = NULL,
  ref_geno_loc = NULL,
  cluster_geno_vcfr = nonimputed_b38
)
saveRDS(correlations_imp_nonimp, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_imputed_vs_unimputed_GSA2022_1044_025_V3.rds')
best_correlations_imp_vs_nonimp <- get_best_correlations(list('GSA2022_1044_025_V3' = correlations_imp_nonimp))

# check batch1 as well
nonimputed_b38_batch1_loc <- '/groups/umcg-franke-scrna/tmp03/projects/venema-2022/processed/genotype/GSA2023_1009/unimputed/GSA2023_1009_rsid_noindels_nodots_onlychr_b38.vcf.gz'# calculate correlations
nonimputed_b38_batch1 <- read.vcfR(nonimputed_b38_batch1_loc)
correlations_imp_nonimp_batch1 <- correlate_genotypes(
  ref_geno_vcfr = ref_geno_all, 
  cluster_geno_loc = NULL,
  ref_geno_loc = NULL,
  cluster_geno_vcfr = nonimputed_b38_batch1
)
saveRDS(correlations_imp_nonimp_batch1, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/demultiplexing/souporcell/assignments/mo_imputed_vs_unimputed_GSA2021_1009.rds')
best_correlations_imp_vs_nonimp_batch1 <- get_best_correlations(list('GSA2021_1009' = correlations_imp_nonimp_batch1))
best_correlations_imp_vs_nonimp_batch1[grep('MO', best_correlations_imp_vs_nonimp_batch1$cluster), ]