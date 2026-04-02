#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_plot_pseudobulk_correlations.R
# Function: plot eQTL/caQTL correlating gene pairs
# Example: 
############################################################################################################################

####################
# libraries        #
####################

# for reading tables
library(data.table)
# for genotype data
library(snpStats)
# for plotting
library(ggplot2)
library(cowplot)


####################
# Functions        #
####################


#' Read Expression Files Per Cell Type
#'
#' This function reads expression files from a specified directory, filters them based on a given pattern, and returns a list of data frames, each corresponding to a different cell type.
#'
#' @param expression_dir A character string specifying the directory containing the expression files.
#' @param expression_file_prepend A character string to prepend to the file names (default is an empty string).
#' @param expression_file_append A character string to append to the file names (default is '.qtlInput.txt.gz').
#'
#' @return A list of data frames, where each data frame contains the expression data for a specific cell type.
#'
#' @examples
#' \dontrun{
#' expression_data <- read_expression_files_per_celltype("path/to/expression/files")
#' }
#'
read_expression_files_per_celltype <- function(expression_dir, expression_file_prepend='', expression_file_append='.qtlInput.txt.gz', row.names=1) {
  # list the files in the directory
  files_directory <- list.files(expression_dir, recursive = F, full.names = F, include.dirs = F)
  # make the pattern to match by
  files_pattern <- paste0('^', expression_file_prepend, '.*', expression_file_append, '$')
  # filter by that pattern
  files_directory <- files_directory[grep(files_pattern, files_directory)]
  # we will save per cell type
  expression_per_celltype <- list()
  # check each file
  for (qtl_file in files_directory) {
    # paste together the full path
    qtl_file_path_full <- paste(expression_dir, qtl_file, sep = '/')
    # read that file
    qtl_file_contents <- read.table(qtl_file_path_full, header = T, sep = '\t', row.names = row.names, check.names = F, comment.char = '')
    # extract the cell type from the file name
    cell_type <- gsub(paste0('^', expression_file_prepend), '', qtl_file)
    cell_type <- gsub(paste0(expression_file_append, '$'), '', cell_type)
    # put in the list under that celltype name
    expression_per_celltype[[cell_type]] <- qtl_file_contents
  }
  return(expression_per_celltype)
}


#' Create QTL Plot Table
#'
#' This function creates a table for QTL plotting by combining feature data, genotype data, and sample mapping information.
#'
#' @param feature A character string specifying the feature to extract from the feature table.
#' @param variant A character string specifying the variant to extract from the genotype data.
#' @param feature_table A data frame containing the feature data with phenotype IDs as column names.
#' @param genotypes A list containing genotype data and a map of alleles.
#' @param sample_mapping A data frame containing the mapping of phenotype IDs to genotype IDs.
#'
#' @return A data frame containing the combined feature data, genotype data, and allele information for plotting.
#'
#' @examples
#' \dontrun{
#' plot_table <- create_qtl_plot_table("feature1", "variant1", feature_table, genotypes, sample_mapping)
#' }
#'
create_qtl_plot_table <- function(feature, variant, feature_table, genotypes, sample_mapping) {
  # extract the feature
  feature_table <- data.frame(
    'phenotype_id' = colnames(feature_table), 
    'feature' = as.vector(unlist(feature_table[feature, ]))
  )
  # join onto the sample mapping table
  full_table <- merge(feature_table, sample_mapping, by = 'phenotype_id')
  # get the samples we have genotype data for
  samples_in_genotypes <- rownames(genotypes$genotypes)
  # subset the feature data to those samples
  full_table <- full_table[full_table[['genotype_id']] %in% samples_in_genotypes, ]
  # now extract the genotypes
  genotypes_samples <- data.frame(genotypes$genotypes[full_table[['genotype_id']], variant])
  # rename the column
  colnames(genotypes_samples) <- c('genotype')
  # convert to numeric
  genotypes_samples[['genotype']] <- as.numeric(genotypes_samples[['genotype']])
  # make a mapping of numeric genotypes to alleles
  alleles <- c(paste0(genotypes$map[variant, 'allele.1'], genotypes$map[variant, 'allele.1']), 
               paste0(genotypes$map[variant, 'allele.1'], genotypes$map[variant, 'allele.2']),
               paste0(genotypes$map[variant, 'allele.2'], genotypes$map[variant, 'allele.2'])
  )
  genotypes_samples[['alleles']] <- alleles[genotypes_samples[['genotype']]]
  # add this information to the feature table
  full_table <- cbind(full_table, genotypes_samples)
  return(full_table)
}


#' Get Color List
#'
#' This function generates a list of colors corresponding to a given vector of names. It ensures each unique name is assigned a unique color.
#'
#' @param vector_of_names A character vector containing the names for which colors are to be generated.
#' @param use_sampling A logical value indicating whether to use sampling when generating colors (default is FALSE).
#' @param color_indices An optional vector of color indices to use for sampling.
#'
#' @return A named list of colors, where each name in the input vector is assigned a unique color.
#'
#' @examples
#' \dontrun{
#' color_list <- get_color_list(c("apple", "banana", "cherry"))
#' }
#'
get_color_list <- function(vector_of_names, use_sampling=F, color_indices=NULL) {
  # get the unique entries
  vector_unique <- unique(vector_of_names)
  # remove any NA
  vector_unique <- vector_unique[!is.na(vector_unique)]
  # get some colors
  colors_to_use <- NULL
  if (length(vector_unique) > 74) {
    colors_to_use <- roycols::sample_tons_of_colors(length(vector_unique), use_sampling = use_sampling, color_indices = color_indices)
  }
  else{
    colors_to_use <- roycols::sample_many_colours(length(vector_unique), use_sampling = use_sampling, color_indices = color_indices)
  }
  # turn into a list
  colors_to_use_list <- as.list(colors_to_use)
  # and add the names
  names(colors_to_use_list) <- vector_unique
  return(colors_to_use_list)
}


#' Plot QTLs
#'
#' This function generates a QTL plot for a given feature and variant, using feature data, genotype data, and sample mapping information.
#'
#' @param feature A character string specifying the feature to plot.
#' @param variant A character string specifying the variant to plot.
#' @param feature_table A data frame containing the feature data with phenotype IDs as column names.
#' @param genotypes A list containing genotype data and a map of alleles.
#' @param sample_mapping A data frame containing the mapping of phenotype IDs to genotype IDs.
#' @param pointless A logical value indicating whether to remove x-axis text and ticks (default is FALSE).
#' @param legendless A logical value indicating whether to remove the legend (default is FALSE).
#' @param ylim A numeric vector specifying the y-axis limits (default is NULL).
#' @param paper_style A logical value indicating whether to apply a paper-style theme (default is TRUE).
#' @param angle_labels A logical value indicating whether to angle the x-axis labels (default is FALSE).
#' @param better_colors A logical value indicating whether to use a better color scheme (default is FALSE).
#'
#' @return A ggplot object representing the QTL plot.
#'
#' @examples
#' \dontrun{
#' qtl_plot <- plot_qtls("feature1", "variant1", feature_table, genotypes, sample_mapping)
#' }
#'
plot_qtls <- function(feature, variant, feature_table, genotypes, sample_mapping, pointless=F, legendless=F, ylim=NULL, paper_style=T, angle_labels=F, better_colors=F) {
  # get the plottable table
  qtl_table <- create_qtl_plot_table(feature, variant, feature_table, genotypes, sample_mapping)
  # turn into plot
  p <- ggplot(data = qtl_table, mapping = aes(x = alleles, y = feature, fill = alleles)) + 
    geom_boxplot(outlier.shape = NA) + 
    # and add jitter
    geom_jitter(size = 0.5, alpha = 0.5) +
    # and labels
    xlab(variant) + 
    ylab(feature) +
    ggtitle(paste(variant, 'affecting', feature))
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
  if (better_colors) {
    p <- p + scale_fill_manual(values = get_color_list(qtl_table[['alleles']]))
  }
  return(p)
}


#' Create Interaction Plot Table
#'
#' This function creates a table for interaction plotting by combining data from two features, genotype data, and sample mapping information.
#'
#' @param feature1 A character string specifying the first feature to extract from the first feature table.
#' @param feature2 A character string specifying the second feature to extract from the second feature table.
#' @param variant A character string specifying the variant to extract from the genotype data.
#' @param feature_table1 A data frame containing the first feature data with phenotype IDs as column names.
#' @param feature_table2 A data frame containing the second feature data with phenotype IDs as column names.
#' @param genotypes A list containing genotype data and a map of alleles.
#' @param sample_mapping1 A data frame containing the mapping of phenotype IDs to genotype IDs for the first feature.
#' @param sample_mapping2 A data frame containing the mapping of phenotype IDs to genotype IDs for the second feature.
#'
#' @return A data frame containing the combined feature data, genotype data, and allele information for interaction plotting.
#'
#' @examples
#' \dontrun{
#' interaction_plot_table <- create_interaction_plot_table("feature1", "feature2", "variant1", feature_table1, feature_table2, genotypes, sample_mapping1, sample_mapping2)
#' }
#'
create_interaction_plot_table <- function(feature1, feature2, variant, feature_table1, feature_table2, genotypes, sample_mapping1, sample_mapping2) {
  # extract the feature
  feature_table1 <- data.frame(
    'phenotype_id' = colnames(feature_table1), 
    'feature1' = as.vector(unlist(feature_table1[feature1, ]))
  )
  # the second one as well
  feature_table2 <- data.frame(
    'phenotype_id' = colnames(feature_table2), 
    'feature2' = as.vector(unlist(feature_table2[feature2, ]))
  )
  # combine them
  feature_table <- merge(feature_table1, feature_table2, by = 'phenotype_id')
  # combine the smf files
  sample_mapping <- unique(rbind(sample_mapping1, sample_mapping2))
  # join onto the sample mapping table
  full_table <- merge(feature_table, sample_mapping, by = 'phenotype_id')
  # get the samples we have genotype data for
  samples_in_genotypes <- rownames(genotypes$genotypes)
  # subset the feature data to those samples
  full_table <- full_table[full_table[['genotype_id']] %in% samples_in_genotypes, ]
  # now extract the genotypes
  genotypes_samples <- data.frame(genotypes$genotypes[full_table[['genotype_id']], variant])
  # rename the column
  colnames(genotypes_samples) <- c('genotype')
  # convert to numeric
  genotypes_samples[['genotype']] <- as.numeric(genotypes_samples[['genotype']])
  # make a mapping of numeric genotypes to alleles
  alleles <- c(paste0(genotypes$map[variant, 'allele.1'], genotypes$map[variant, 'allele.1']), 
               paste0(genotypes$map[variant, 'allele.1'], genotypes$map[variant, 'allele.2']),
               paste0(genotypes$map[variant, 'allele.2'], genotypes$map[variant, 'allele.2'])
  )
  genotypes_samples[['alleles']] <- alleles[genotypes_samples[['genotype']]]
  # add this information to the feature table
  full_table <- cbind(full_table, genotypes_samples)
  return(full_table)
}


#' Plot Interactions
#'
#' This function generates an interaction plot for two features and a variant, using feature data, genotype data, and sample mapping information.
#'
#' @param feature1 A character string specifying the first feature to plot.
#' @param feature2 A character string specifying the second feature to plot.
#' @param variant A character string specifying the variant to plot.
#' @param feature_table1 A data frame containing the first feature data with phenotype IDs as column names.
#' @param feature_table2 A data frame containing the second feature data with phenotype IDs as column names.
#' @param genotypes A list containing genotype data and a map of alleles.
#' @param sample_mapping1 A data frame containing the mapping of phenotype IDs to genotype IDs for the first feature.
#' @param sample_mapping2 A data frame containing the mapping of phenotype IDs to genotype IDs for the second feature.
#' @param pointless A logical value indicating whether to remove x-axis text and ticks (default is FALSE).
#' @param legendless A logical value indicating whether to remove the legend (default is FALSE).
#' @param ylim A numeric vector specifying the y-axis limits (default is NULL).
#' @param paper_style A logical value indicating whether to apply a paper-style theme (default is TRUE).
#' @param angle_labels A logical value indicating whether to angle the x-axis labels (default is FALSE).
#' @param better_colors A logical value indicating whether to use a better color scheme (default is FALSE).
#'
#' @return A ggplot object representing the interaction plot.
#'
#' @examples
#' \dontrun{
#' interaction_plot <- plot_interactions("feature1", "feature2", "variant1", feature_table1, feature_table2, genotypes, sample_mapping1, sample_mapping2)
#' }
#'
plot_interactions <- function(feature1, feature2, variant, feature_table1, feature_table2, genotypes, sample_mapping1, sample_mapping2, pointless=F, legendless=F, ylim=NULL, paper_style=T, angle_labels=F, better_colors=F) {
  # get the plottable table
  interaction_table <- create_interaction_plot_table(feature1, feature2, variant, feature_table1, feature_table2, genotypes, sample_mapping1, sample_mapping2)
  # turn into plot
  p <- ggplot(data = interaction_table, mapping = aes(x = feature1, y = feature2, colour = alleles)) + 
    geom_point() +
    # and labels
    xlab(feature1) + 
    ylab(feature2) +
    ggtitle(paste(variant, 'affecting\n', feature1, 'and\n', feature2))
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
  if (better_colors) {
    p <- p + scale_colour_manual(values = get_color_list(interaction_table[['alleles']]))
  }
  return(p)
}


get_color_coding_dict <- function() {
  # medhigh
  color_coding_dict <- list()
  color_coding_dict[["B"]] <- "#71BC4B"
  #color_coding_dict[['CD4_T_cells']] <- '#7FC97F'
  color_coding_dict[['CD4_T_cells']] <- '#153057'
  color_coding_dict[['CD4T']] <- '#153057'
  #color_coding_dict[['CD8_T_cells']] <- '#BEAED4'
  color_coding_dict[['CD8_T_cells']] <- '#009DDB'
  color_coding_dict[['CD8T']] <- '#009DDB'
  #color_coding_dict[['Dendritic_cells']] <- '#FDC086'
  color_coding_dict[['Dendritic_cells']] <- '#965EC8'
  color_coding_dict[['DC']] <- '#965EC8'
  color_coding_dict[['Endothelial_cells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia_cells']] <- '#F0027F'
  color_coding_dict[['Mast_cells']] <- '#BF5B17'
  color_coding_dict[['Mature_absorptive_enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature_secretory_enterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory_B']] <- '#D95F02'
  color_coding_dict[['Microfold_cell']] <- '#BEAED4'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Monocyte']] <- '#EDBA1B'
  color_coding_dict[['Naive_B_cells']] <- '#FDC086'
  color_coding_dict[['NK']] <- '#E64B50'
  #color_coding_dict[['Plasma_cells']] <- '#E7298A'
  color_coding_dict[['Plasma_cells']] <- '#DB8E00'
  color_coding_dict[['Stem_cells']] <- '#66A61E'
  color_coding_dict[['Stromal_cells']] <- '#8DD3C7'
  #color_coding_dict[['T_others']] <- '#A6761D'
  color_coding_dict[['T_others']] <- '#FF63B6'
  color_coding_dict[['Transit_amplifying_cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  #color_coding_dict[['CD4+ T cells']] <- '#7FC97F'
  color_coding_dict[['CD4+ T cells']] <- '#153057'
  color_coding_dict[['CD4+ T']] <- '#153057'
  #color_coding_dict[['CD8+ T cells']] <- '#BEAED4'
  color_coding_dict[['CD8+ T cells']] <- '#009DDB'
  color_coding_dict[['CD8+ T']] <- '#009DDB'
  #color_coding_dict[['Dendritic cells']] <- '#FDC086'
  color_coding_dict[['Dendritic cells']] <- '#965EC8'
  color_coding_dict[['Endothelial cells']] <- '#FFFFB3'
  color_coding_dict[['Endothelial\ncells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia cells']] <- '#F0027F'
  color_coding_dict[['MAST cells']] <- '#BF5B17'
  color_coding_dict[['Mature absorptive enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature\nabsorptive\nenterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature secretory enterocytes']] <- '#1B9E77'
  color_coding_dict[['Mature secretory\nenterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory B cells']] <- '#D95F02'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Microfold cells']] <- '#BEAED4'
  color_coding_dict[['Monocytes']] <- '#EDBA1B'
  color_coding_dict[['Naive B cells']] <- '#FDC086'
  #color_coding_dict[['Plasma cells']] <- '#E7298A'
  color_coding_dict[['Plasma cells']] <- '#DB8E00'
  color_coding_dict[['Stem cells']] <- '#66A61E'
  color_coding_dict[['Stromal cells']] <- '#8DD3C7'
  #color_coding_dict[['other T cells']] <- '#A6761D'
  color_coding_dict[['other T cells']] <- '#FF63B6'
  color_coding_dict[['Transit amplifying cells']] <- '#FF7F00'
  color_coding_dict[['Transit\namplifying cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  color_coding_dict[['UT']] <- 'gray'
  color_coding_dict[['24hCA']] <- 'darkgreen'
  # up and down regulation will be added to, we need a whitening percentage
  pct_whitening <- 40
  # then we will check each cell type
  for (cell_type in names(color_coding_dict)) {
    # the up color is the same as the regular one
    color_coding_dict[[paste(cell_type, 'up')]] <- color_coding_dict[[cell_type]]
    # but the down one will have a more faded colour
    color_coding_dict[[paste(cell_type, 'down')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    # we'll do something similiar when we have multiple conditions
    color_coding_dict[[paste(cell_type, 'combined')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'UT')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, '24hCA')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "black"))(100)[pct_whitening]
  }
  # general
  color_coding_dict[['AI']] <- 'darkblue'
  color_coding_dict[['NI']] <- 'darkred'
  color_coding_dict[['Actively Inflamed']] <- 'darkblue'
  color_coding_dict[['Non-Inflamed']] <- 'darkred'
  return(color_coding_dict)
}



####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- T


####################
# Main code        #
####################

# location of the expression data
exp_data_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/input/L1/combined/'
# and the sample mapping files
eqtl_smf_loc <- paste(exp_data_loc, '/smf.txt', sep = '/')
# read the sample mapping files
eqtl_smf <- read.table(eqtl_smf_loc, header = T, sep = '\t')
# read the expression data
eqtl_inputs <- read_expression_files_per_celltype(exp_data_loc)
# read the cell count data as well
eqtl_covariates <- read_expression_files_per_celltype(exp_data_loc, expression_file_append = '.covariates.txt.gz', row.names = NULL)
# read the proportion expressed as well
prop_expressed_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/frac_exp/'
eqtl_ct_props <- read_expression_files_per_celltype(prop_expressed_loc, expression_file_append = '_persample.tsv.gz', row.names = 1)
# read the metadata
mtdt_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
mtdt <- fread(mtdt_loc, header = T, sep = '\t')
# add to the smf the condition
eqtl_smf[['condition']] <- mtdt[match(eqtl_smf[['phenotype_id']], paste(mtdt[['sample_final']], mtdt[['lane']], sep = ';;'))][['condition_final']]
# which genotype to use
genotype_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/genotype/imputed_hg38_all_anc_qtl_tested_variants'
# read the genotypes, but only those in the file and in the confinement
genotypes <- read.plink(
  bed = paste(genotype_loc, '.bed', sep = ''),
  bim = paste(genotype_loc, '.bim', sep = ''),
  fam = paste(genotype_loc, '.fam', sep = '')
)
# CD4T UT
cor(
  x = as.vector(
    unlist(
      eqtl_inputs$CD4T['BLK', 
                       intersect(colnames(eqtl_inputs$CD4T), 
                                 eqtl_smf[
                                   !is.na(eqtl_smf[['condition']]) &
                                     eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
  y = as.vector(
    unlist(
      eqtl_inputs$CD4T['FAM167A', 
                       intersect(colnames(eqtl_inputs$CD4T), 
                                 eqtl_smf[
                                   !is.na(eqtl_smf[['condition']]) &
                                     eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
  method = 'spearman'
)
# -0.02010506
# CD4T 24hCA
cor(
  x = as.vector(
    unlist(
      eqtl_inputs$CD4T['BLK', 
                       intersect(colnames(eqtl_inputs$CD4T), 
                                 eqtl_smf[
                                   !is.na(eqtl_smf[['condition']]) &
                                     eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
  y = as.vector(
    unlist(
      eqtl_inputs$CD4T['FAM167A', 
                       intersect(colnames(eqtl_inputs$CD4T), 
                                 eqtl_smf[
                                   !is.na(eqtl_smf[['condition']]) &
                                     eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
  method = 'spearman'
)
# 0.0310517
# B UT
cor(
  x = as.vector(
    unlist(
      eqtl_inputs$B['BLK', 
                    intersect(colnames(eqtl_inputs$B), 
                              eqtl_smf[
                                !is.na(eqtl_smf[['condition']]) &
                                  eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
  y = as.vector(
    unlist(
      eqtl_inputs$B['FAM167A', 
                    intersect(colnames(eqtl_inputs$B), 
                              eqtl_smf[
                                !is.na(eqtl_smf[['condition']]) &
                                  eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
  method = 'spearman'
)
# -0.1068247
# B 24hCA
cor(
  x = as.vector(
    unlist(
      eqtl_inputs$B['BLK', 
                    intersect(colnames(eqtl_inputs$B), 
                              eqtl_smf[
                                !is.na(eqtl_smf[['condition']]) &
                                  eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
  y = as.vector(
    unlist(
      eqtl_inputs$B['FAM167A', 
                    intersect(colnames(eqtl_inputs$B), 
                              eqtl_smf[
                                !is.na(eqtl_smf[['condition']]) &
                                  eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
  method = 'spearman'
)
# -0.1382585

# create a frame to plot
plot_df <- rbind(
  data.frame(
    'BLK' = as.vector(
      unlist(
        eqtl_inputs$B['BLK', 
                      intersect(colnames(eqtl_inputs$B), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
    'FAM167A' = as.vector(
      unlist(
        eqtl_inputs$B['FAM167A', 
                      intersect(colnames(eqtl_inputs$B), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
    'phenotype_id' = intersect(colnames(eqtl_inputs$B), 
                               eqtl_smf[
                                 !is.na(eqtl_smf[['condition']]) &
                                   eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']]), 
    'condition' = '24hCA'
  ), 
  data.frame(
    'BLK' = as.vector(
      unlist(
        eqtl_inputs$B['BLK', 
                      intersect(colnames(eqtl_inputs$B), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
    'FAM167A' = as.vector(
      unlist(
        eqtl_inputs$B['FAM167A', 
                      intersect(colnames(eqtl_inputs$B), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
    'phenotype_id' = intersect(colnames(eqtl_inputs$B), 
                               eqtl_smf[
                                 !is.na(eqtl_smf[['condition']]) &
                                   eqtl_smf[['condition']] == 'UT', ][['phenotype_id']]), 
    'condition' = 'UT'
  )
)
# add genotype id
plot_df[['genotype_id']] <- eqtl_smf[match(plot_df[['phenotype_id']], eqtl_smf[['phenotype_id']]), ][['genotype_id']]
# add info for variant
genotype <- genotypes$genotypes[plot_df[['genotype_id']], '8:11491452:G:A']
# then to numeric
genotype_numeric <- as.vector(as(genotype, 'numeric'))
# add the genotype
plot_df[[gsub(':', '_', '8:11491452:G:A')]] <- genotype_numeric
# and as string
plot_df[[gsub(':', '_', '8:11491452:G:A_string')]] <- as.character(plot_df[[gsub(':', '_', '8:11491452:G:A')]])
# set the order
plot_df[['condition']] <- factor(plot_df[['condition']], c('UT', '24hCA'))
# add cell count as well
plot_df[['cell_count']] <- eqtl_covariates[['B']][match(plot_df[['phenotype_id']], eqtl_covariates[['B']][['Donor_Pool']]), ][['CellCount']]
# add the proportions
prot_df_props <- data.frame('phenotype_id' = colnames(eqtl_ct_props[['B']]), 'prop_FAM167A' = as.vector(unlist(eqtl_ct_props[['B']]['FAM167A', ])), 'prop_BLK' = as.vector(unlist(eqtl_ct_props[['B']]['BLK', ])))
plot_df <- merge(plot_df, prot_df_props, all.x = T, by = 'phenotype_id')

# and plot that
ggplot(data = plot_df, mapping = aes(x = BLK, y = FAM167A, colour = condition)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  # and colour of ncell
  scale_colour_manual(values = get_color_coding_dict()) + 
  xlab(paste('BLK', 'expression')) + 
  ylab(paste('FAM167A', 'expression')) + 
  labs(fill = 'condition', colour = 'condition') + 
  ggtitle(paste('BLK', 'FAM167A', 'co-expression')) +
  theme(legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
ggsave('~/multiome/plots/BLK_FAM167A_coexpression_B.pdf', width = 6, height = 6)
# and plot that
ggplot(data = plot_df, mapping = aes(x = BLK, y = FAM167A, colour = `8_11491452_G_A_string`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  # and colour of ncell
  scale_colour_manual(values = roycols::get_color_list(plot_df[['8_11491452_G_A_string']])) + 
  xlab(paste('BLK', 'expression')) + 
  ylab(paste('FAM167A', 'expression')) + 
  labs(fill = 'genotype', colour = 'genotype') + 
  ggtitle(paste('BLK', 'FAM167A', 'co-expression')) +
  theme(legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
plot_grid(
  ggplot(data = plot_df[plot_df[['condition']] == 'UT', ], mapping = aes(x = BLK, y = FAM167A, colour = `8_11491452_G_A_string`)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    # and colour of ncell
    scale_colour_manual(values = roycols::get_color_list(plot_df[['8_11491452_G_A_string']])) + 
    xlab(paste('BLK', 'expression')) + 
    ylab(paste('FAM167A', 'expression')) + 
    labs(fill = 'genotype', colour = 'genotype') + 
    ggtitle(paste('BLK', 'FAM167A', 'co-expression', 'UT')) +
    theme(legend.title = element_text(size=14), 
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(size=12)) + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")), 
  ggplot(data = plot_df[plot_df[['condition']] == '24hCA', ], mapping = aes(x = BLK, y = FAM167A, colour = `8_11491452_G_A_string`)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    # and colour of ncell
    scale_colour_manual(values = roycols::get_color_list(plot_df[['8_11491452_G_A_string']])) + 
    xlab(paste('BLK', 'expression')) + 
    ylab(paste('FAM167A', 'expression')) + 
    labs(fill = 'genotype', colour = 'genotype') + 
    ggtitle(paste('BLK', 'FAM167A', 'co-expression', '24hCA')) +
    theme(legend.title = element_text(size=14), 
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(size=12)) + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
)
ggplot(data = plot_df, mapping = aes(x = BLK, y = FAM167A, colour = `cell_count`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  # and colour of ncell
  scale_colour_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = 0) +
  xlab(paste('BLK', 'expression')) + 
  ylab(paste('FAM167A', 'expression')) + 
  labs(fill = 'genotype', colour = 'ncell') + 
  ggtitle(paste('BLK', 'FAM167A', 'co-expression')) +
  theme(legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
plot_grid(
  ggplot(data = plot_df, mapping = aes(x = BLK, y = FAM167A, colour = `prop_FAM167A`)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    # and colour of ncell
    scale_colour_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = 0) +
    xlab(paste('BLK', 'expression')) + 
    ylab(paste('FAM167A', 'expression')) + 
    labs(fill = 'genotype', colour = 'prop_FAM167A') + 
    ggtitle(paste('BLK', 'FAM167A', 'co-expression')) +
    theme(legend.title = element_text(size=14), 
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(size=12)) + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")), 
  ggplot(data = plot_df, mapping = aes(x = BLK, y = FAM167A, colour = `prop_BLK`)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    # and colour of ncell
    scale_colour_gradient2(low = 'darkblue', mid = 'white', high = 'darkred', midpoint = 0) +
    xlab(paste('BLK', 'expression')) + 
    ylab(paste('FAM167A', 'expression')) + 
    labs(fill = 'genotype', colour = 'prop_BLK') + 
    ggtitle(paste('BLK', 'FAM167A', 'co-expression')) +
    theme(legend.title = element_text(size=14), 
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(size=12)) + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
)

# create a frame to plot
plot_df_cd4t <- rbind(
  data.frame(
    'BLK' = as.vector(
      unlist(
        eqtl_inputs$CD4T['BLK', 
                      intersect(colnames(eqtl_inputs$CD4T), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
    'FAM167A' = as.vector(
      unlist(
        eqtl_inputs$CD4T['FAM167A', 
                      intersect(colnames(eqtl_inputs$CD4T), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']])])), 
    'phenotype_id' = intersect(colnames(eqtl_inputs$CD4T), 
                               eqtl_smf[
                                 !is.na(eqtl_smf[['condition']]) &
                                   eqtl_smf[['condition']] == '24hCA', ][['phenotype_id']]), 
    'condition' = '24hCA'
  ), 
  data.frame(
    'BLK' = as.vector(
      unlist(
        eqtl_inputs$CD4T['BLK', 
                      intersect(colnames(eqtl_inputs$CD4T), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
    'FAM167A' = as.vector(
      unlist(
        eqtl_inputs$CD4T['FAM167A', 
                      intersect(colnames(eqtl_inputs$CD4T), 
                                eqtl_smf[
                                  !is.na(eqtl_smf[['condition']]) &
                                    eqtl_smf[['condition']] == 'UT', ][['phenotype_id']])])), 
    'phenotype_id' = intersect(colnames(eqtl_inputs$CD4T), 
                               eqtl_smf[
                                 !is.na(eqtl_smf[['condition']]) &
                                   eqtl_smf[['condition']] == 'UT', ][['phenotype_id']]), 
    'condition' = 'UT'
  )
)
# add genotype id
plot_df_cd4t[['genotype_id']] <- eqtl_smf[match(plot_df_cd4t[['phenotype_id']], eqtl_smf[['phenotype_id']]), ][['genotype_id']]
# add info for variant
genotype <- genotypes$genotypes[plot_df_cd4t[['genotype_id']], '8:11491452:G:A']
# then to numeric
genotype_numeric <- as.vector(as(genotype, 'numeric'))
# add the genotype
plot_df_cd4t[[gsub(':', '_', '8:11491452:G:A')]] <- genotype_numeric
# and as string
plot_df_cd4t[[gsub(':', '_', '8:11491452:G:A_string')]] <- as.character(plot_df_cd4t[[gsub(':', '_', '8:11491452:G:A')]])
# set the order
plot_df_cd4t[['condition']] <- factor(plot_df_cd4t[['condition']], c('UT', '24hCA'))
# add cell count as well
plot_df_cd4t[['cell_count']] <- eqtl_covariates[['B']][match(plot_df_cd4t[['phenotype_id']], eqtl_covariates[['B']][['Donor_Pool']]), ][['CellCount']]
# add the proportions
prot_df_props <- data.frame('phenotype_id' = colnames(eqtl_ct_props[['B']]), 'prop_FAM167A' = as.vector(unlist(eqtl_ct_props[['B']]['FAM167A', ])), 'prop_BLK' = as.vector(unlist(eqtl_ct_props[['B']]['BLK', ])))
plot_df_cd4t <- merge(plot_df_cd4t, prot_df_props, all.x = T, by = 'phenotype_id')
# and plot that
ggplot(data = plot_df_cd4t, mapping = aes(x = BLK, y = FAM167A, colour = condition)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  # and colour of ncell
  scale_colour_manual(values = get_color_coding_dict()) + 
  xlab(paste('BLK', 'expression')) + 
  ylab(paste('FAM167A', 'expression')) + 
  labs(fill = 'condition', colour = 'condition') + 
  ggtitle(paste('BLK', 'FAM167A', 'co-expression')) +
  theme(legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
ggsave('~/multiome/plots/BLK_FAM167A_coexpression_CD4T.pdf', width = 6, height = 6)
