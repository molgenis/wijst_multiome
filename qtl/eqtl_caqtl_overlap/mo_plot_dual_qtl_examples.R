#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_plot_dual_qtl_examples.R
# Function: plot top eQTL/caQTL pairs nominated by Maryna
# Example: 
# Rscript 
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
read_expression_files_per_celltype <- function(expression_dir, expression_file_prepend='', expression_file_append='.qtlInput.txt.gz') {
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
    qtl_file_contents <- read.table(qtl_file_path_full, header = T, sep = '\t', row.names = 1, check.names = F, comment.char = '')
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

# the confinement
confinement_loc <- '~/tables/mo_dual_qtl_plot_tbl.tsv.gz'
# read the confinement file
confinement <- fread(confinement_loc, header = T, sep = '\t')

# the variant data
genotype_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/genotype/imputed_hg38_all_anc'
# get the variants from the confinement file
variants <- confinement[['variant']]
# read the bim
variants_in_gt <- fread(paste(genotype_loc, '.bim', sep = ''), header = F)[[2]]
# get overlapping variants
overlapping_variants <- intersect(variants, variants_in_gt)
# read the genotypes, but only those in the file and in the confinement
genotypes <- read.plink(
  bed = paste(genotype_loc, '.bed', sep = ''),
  bim = paste(genotype_loc, '.bim', sep = ''),
  fam = paste(genotype_loc, '.fam', sep = ''), 
  select.snps = overlapping_variants
)

# location of the expression data
exp_data_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/input/L1/combined/'
acc_data_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/'
# and the sample mapping files
eqtl_smf_loc <- paste(exp_data_loc, '/smf.txt', sep = '/')
caqtl_smf_loc <- paste(acc_data_loc, '/smf.txt', sep = '/')
# read the sample mapping files
eqtl_smf <- read.table(eqtl_smf_loc, header = T, sep = '\t')
caqtl_smf <- read.table(caqtl_smf_loc, header = T, sep = '\t')
# read the expression data
eqtl_inputs = read_expression_files_per_celltype(exp_data_loc)
caqtl_inputs = read_expression_files_per_celltype(acc_data_loc)

# we'll put the plots in a list
dual_qtl_plots <- list()
# check each of the celltype, variant, gene, region and celltype combination
for (i in 1 : nrow(confinement)) {
  # extract the variant, gene, region and cell type
  variant <- as.vector(unlist(confinement[i, 'variant'][[1]]))
  gene <- as.vector(unlist(confinement[i, 'gene'][[1]]))
  region <- as.vector(unlist(confinement[i, 'region'][[1]]))
  cell_type <- as.vector(unlist(confinement[i, 'cell_type'][[1]]))
  # make the gene plot
  p_gene <- plot_qtls(gene, 
                      variant, 
                      feature_table = eqtl_inputs[[cell_type]], 
                      genotypes = genotypes, 
                      sample_mapping = eqtl_smf, 
                      better_colors = T, legendless = T) + 
    ggtitle(paste(variant, 'affecting\n', gene, 'expression in', cell_type)) + xlab(paste(variant, 'genotype')) + ylab(paste(gene, 'expression'))
  # the accessibility plot
  p_region <- plot_qtls(region, 
                        variant, 
                        feature_table = caqtl_inputs[[cell_type]], 
                        genotypes = genotypes, 
                        sample_mapping = caqtl_smf, 
                        better_colors = T, legendless = T) + 
    ggtitle(paste(variant, 'affecting\n', region, 'accessibility in', cell_type)) + xlab(paste(variant, 'genotype')) + ylab(paste(region, 'accessibility'))
  # and both
  p_both <- plot_grid(
    p_region, p_gene,
    nrow = 1,
    ncol = 2
  )
  # put in the list
  dual_qtl_plots[[paste(variant, gene, region, cell_type, sep = '_')]] <- p_both
}

# let's try to save the plots
dual_plots_loc <- '~/plots/multiome/dual_qtls/'
# make that directory
dir.create(dual_plots_loc, recursive = T)
# now do each combination
for (comb in names(dual_qtl_plots)) {
  # extract the plot
  p_to_plot <- dual_qtl_plots[[comb]]
  # make a safer name
  p_name <- gsub('\\:', '.', comb)
  # save the plot
  ggsave(paste0(dual_plots_loc, '/', p_name, '.pdf'), plot = p_to_plot, width = 6, height = 6)
}
