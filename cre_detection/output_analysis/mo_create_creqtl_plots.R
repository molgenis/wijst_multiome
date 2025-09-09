#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_creqtl_plots.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

# for reading the data
library(data.table)
library(snpStats)
library(ggplot2)
library(roycols)
# transformation into gaussian normal distribution
library(bestNormalize)


####################
# Functions        #
####################

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
  # extract the feature location
  feature_locs <- feature_table[['feature']] == feature
  # extract the feature
  feature_table <- data.frame(
    'phenotype_id' = colnames(feature_table), 
    'feature' = as.numeric(as.vector(unlist(feature_table[feature_locs, ])))
  )
  # join onto the sample mapping table
  full_table <- merge(feature_table, sample_mapping, on = 'phenotype_id')
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
plot_qtls <- function(feature, variant, feature_table, genotypes, sample_mapping, pointless=F, legendless=F, ylim=NULL, paper_style=T, angle_labels=F, better_colors=F, plot_table=NULL) {
  # get the plottable table
  qtl_table <- NULL
  if (!is.null(plot_table)) {
    qtl_table <- plot_table
  }
  else {
    qtl_table <- create_qtl_plot_table(feature, variant, feature_table, genotypes, sample_mapping)
  }
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

gausnorm_independent_variable_matrix <- function(independent_variable_matrix) {
  # take the donor names, as they are the columns
  colnames_original <- colnames(independent_variable_matrix)
  # first make a list of all the transformed data
  transformed_list <- list()
  # and if they were numeric
  transformed_numeric <- list()
  # check each column
  for (column in colnames(independent_variable_matrix)) {
    # check if numeric
    if (is.numeric(independent_variable_matrix[[column]])) {
      # if so, do a gausnorm
      transformed_list[[column]] <- yeojohnson(independent_variable_matrix[[column]])$x.t
      transformed_numeric[[column]] <- T
    }
    else {
      # otherwise return as-is
      transformed_list[[column]] <- independent_variable_matrix[[column]]
      transformed_numeric[[column]] <- F
    }
  }
  # merge all data
  transformed_data <- data.frame(do.call('cbind', transformed_list))
  # make numeric what was numeric
  for (column in names(transformed_numeric)) {
    if (transformed_numeric[[column]]) {
      transformed_data[[column]] <- as.numeric(transformed_data[[column]])
    }
  }
  return(transformed_data)
}


load_sc_exp_acc <- function(chunk_loc, exp_file='expression.tsv.gz', acc_folder='accessibility.tsv.gz') {
  # read the expression
  exp_table <- fread(paste(chunk_loc, exp_file, sep = '/'), header = F, sep = '\t', skip = 1)
  # and the accessibility
  acc_table <- fread(paste(chunk_loc, acc_folder, sep = '/'), header = F, sep = '\t', skip = 1)
  # read the barcodes as well
  con_exp <- file(paste(chunk_loc, exp_file, sep = '/'),"r")
  exp_bcs <- readLines(con_exp,n=1)
  close(con_exp)
  con_acc <- file(paste(chunk_loc, acc_folder, sep = '/'),"r")
  acc_bcs <- readLines(con_acc,n=1)
  close(con_acc)
  # split into parts
  exp_bcs <- strsplit(exp_bcs, '\t')[[1]]
  acc_bcs <- strsplit(acc_bcs, '\t')[[1]]
  # set those column names
  colnames(exp_table) <- c('feature', exp_bcs)
  colnames(acc_table) <- c('feature', acc_bcs)
  # get the barcodes that we have in both of them
  bc_shared <- intersect(exp_bcs, acc_bcs)
  # store the features
  exp_features <- exp_table[['feature']]
  acc_features <- acc_table[['feature']]
  # subset the matrices to the shared barcodes
  exp_table <- exp_table[, ..bc_shared]
  acc_table <- acc_table[, ..bc_shared]
  # transpose the two matrices
  exp_table <- t(exp_table)
  acc_table <- t(acc_table)
  # make them into dataframes
  exp_table <- data.frame(exp_table)
  acc_table <- data.frame(acc_table)
  # set the features as column names
  colnames(exp_table) <- exp_features
  colnames(acc_table) <- acc_features
  # cbind them, with the barcodes as the first column
  full_table <- cbind(data.frame('bc' = bc_shared), exp_table, acc_table)
  return(full_table)
}


get_sc_acc_vs_exp <- function(feature_table, acc_column, exp_column, variant, genotypes, sample_mapping) {
  # join onto the sample mapping table
  full_table <- merge(data.table(feature_table), data.table(sample_mapping), by.x = 'bc', by.y = 'phenotype_id')
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
  full_table <- cbind(full_table, data.table(genotypes_samples))
  # get only the information we need
  exp_acc_table <- full_table[, c('bc', 'genotype_id', 'condition', 'alleles', ..acc_column, ..exp_column)]
  # rename columns
  colnames(exp_acc_table) <- c('bc', 'genotype_id', 'condition', 'alleles', 'accessibility', 'expression')
  return(exp_acc_table)
}

####################
# settings         #
####################

set.seed(7777)


####################
# debug code      #
####################


####################
# Main Code        #
####################

# location of the genotypes
genotype_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/genotype/imputed_hg38_all_anc'
# location of the beta file
beta_file_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/creqtl_replication/betas_ps/CD4T/betas.tsv.gz'
# location of the sample mapping
smf_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/creqtl_replication/matrices/CD4T/sample_mapping.tsv.gz'
# location of the ncells
ncell_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/creqtl_replication/matrices/CD4T/ncells.tsv.gz'

# read the genotype data
genotypes <- snpStats::read.plink(genotype_loc)

# read the number of cells matrix
ncells <- fread(ncell_loc, header = T, sep = '\t')

# read the beta file
betas <- fread(beta_file_loc, header = T, sep = '\t')
# add feature column
betas[['feature']] <- paste(betas[['region']], betas[['gene']], sep = '_')
# and remove original feature and gene
betas[['region']] <- NULL
betas[['gene']] <- NULL

# read sample mapping
smf <- fread(smf_loc, header = T, sep = '\t')
# rename to expected format
colnames(smf) <- c('phenotype_id', 'genotype_id', 'batch')

# get a plottable table
beta_plottable1 <- create_qtl_plot_table(feature = 'chr18-8415530-8417057_PTPRM', variant = '18:7780627:G:A', feature_table = betas, genotypes = genotypes, sample_mapping = smf)
beta_plottable2 <- create_qtl_plot_table(feature = 'chr2-6915372-6915767_NRIR', variant = '2:6849889:C:T', feature_table = betas, genotypes = genotypes, sample_mapping = smf)
beta_plottable3 <- create_qtl_plot_table(feature = 'chr10-62261744-62261937_RTKN2', variant = '10:62264836:AT:A', feature_table = betas, genotypes = genotypes, sample_mapping = smf)
beta_plottable4 <- create_qtl_plot_table(feature = 'chr11-376295-377113_PSMD13', variant = '11:214169:T:C', feature_table = betas, genotypes = genotypes, sample_mapping = smf)

# gaussnorm the data
beta_plottable1_gn <- gausnorm_independent_variable_matrix(beta_plottable1[!is.na(beta_plottable1[['feature']]), ])
beta_plottable2_gn <- gausnorm_independent_variable_matrix(beta_plottable2[!is.na(beta_plottable2[['feature']]), ])
beta_plottable3_gn <- gausnorm_independent_variable_matrix(beta_plottable3[!is.na(beta_plottable3[['feature']]), ])
beta_plottable4_gn <- gausnorm_independent_variable_matrix(beta_plottable4[!is.na(beta_plottable4[['feature']]), ])
# remove bigger than 3 Z
beta_plottable1_gn <- beta_plottable1_gn[abs(beta_plottable1_gn[['feature']]) < 3, ]
beta_plottable2_gn <- beta_plottable2_gn[abs(beta_plottable2_gn[['feature']]) < 3, ]
beta_plottable3_gn <- beta_plottable3_gn[abs(beta_plottable3_gn[['feature']]) < 3, ]
beta_plottable4_gn <- beta_plottable4_gn[abs(beta_plottable4_gn[['feature']]) < 3, ]
# and inverse rank
beta_plottable1_ir <- beta_plottable1
beta_plottable1_ir[['feature']] <- qnorm((rank(beta_plottable1_ir[['feature']],na.last="keep")-0.5)/sum(!is.na(beta_plottable1_ir[['feature']])))
beta_plottable2_ir <- beta_plottable2
beta_plottable2_ir[['feature']] <- qnorm((rank(beta_plottable2_ir[['feature']],na.last="keep")-0.5)/sum(!is.na(beta_plottable2_ir[['feature']])))
beta_plottable3_ir <- beta_plottable3
beta_plottable3_ir[['feature']] <- qnorm((rank(beta_plottable3_ir[['feature']],na.last="keep")-0.5)/sum(!is.na(beta_plottable3_ir[['feature']])))
beta_plottable4_ir <- beta_plottable4
beta_plottable4_ir[['feature']] <- qnorm((rank(beta_plottable4_ir[['feature']],na.last="keep")-0.5)/sum(!is.na(beta_plottable4_ir[['feature']])))


# plot this
plot_qtls(feature = 'chr18-8415530-8417057_PTPRM', variant = '18:7780627:G:A', feature_table = betas, genotypes = genotypes, sample_mapping = smf, better_colors = T, plot_table = beta_plottable1_gn)
plot_qtls(feature = 'chr2-6915372-6915767_NRIR', variant = '2:6849889:C:T', feature_table = betas, genotypes = genotypes, sample_mapping = smf, better_colors = T, plot_table = beta_plottable2_gn)
plot_qtls(feature = 'chr10-62261744-62261937_RTKN2', variant = '10:62264836:AT:A', feature_table = betas, genotypes = genotypes, sample_mapping = smf, better_colors = T, plot_table = beta_plottable3_gn)
plot_qtls(feature = 'chr11-376295-377113_PSMD13', variant = '11:214169:T:C', feature_table = betas, genotypes = genotypes, sample_mapping = smf, better_colors = T, plot_table = beta_plottable4_gn)

# filter to donors with many cells
columns_ncell300 <- c('feature', intersect(colnames(betas), ncells[ncells[['ncell']] >= 300, ][['sample']]))
beta_plottable_ncell300 <- betas[, ..columns_ncell300]
# and plot again
plot_qtls(feature = 'chr18-8415530-8417057_PTPRM', variant = '18:7780627:G:A', feature_table = beta_plottable_ncell300, genotypes = genotypes, sample_mapping = smf, better_colors = T)
plot_qtls(feature = 'chr2-6915372-6915767_NRIR', variant = '2:6849889:C:T', feature_table = beta_plottable_ncell300, genotypes = genotypes, sample_mapping = smf, better_colors = T)
plot_qtls(feature = 'chr10-62261744-62261937_RTKN2', variant = '10:62264836:AT:A', feature_table = beta_plottable_ncell300, genotypes = genotypes, sample_mapping = smf, better_colors = T)


# read the location of the expression metadata
sc_metadata_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
# read the metadata file
sc_metadata <- fread(sc_metadata_loc, header = T, sep = '\t')

# get the single-cell smf
smf_sc <- data.table('phenotype_id' = sc_metadata[['barcode_lane']], 'genotype_id' = sc_metadata[['sample_final']], 'condition' = sc_metadata[['final_condition']], 'lane' = sc_metadata[['lane']])

# read the information for a specific chunk
sc_exp_acc4 <- load_sc_exp_acc('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD4T/chr11-125926-590601/')
# add genotype info
sc_exp_acc4_plottable <- get_sc_acc_vs_exp(feature_table = sc_exp_acc4, acc_column = 'chr11-376295-377113', exp_column = 'PSMD13', variant = '11:214169:T:C', genotypes = genotypes, sample_mapping = smf_sc)
# remove any unknown conditions
sc_exp_acc4_plottable <- sc_exp_acc4_plottable[!is.na(sc_exp_acc4_plottable[['condition']])]
# and genotype
sc_exp_acc4_plottable[['condition_openness_genotype']] <- paste(sc_exp_acc4_plottable[['condition']], sc_exp_acc4_plottable[['accessibility']], sc_exp_acc4_plottable[['alleles']])
# show that
ggplot(data = sc_exp_acc4_plottable, mapping = aes(x = condition_openness_genotype, y = expression, fill = condition_openness_genotype)) + 
  geom_boxplot(outlier.shape = NA) + 
  # and add jitter
  geom_jitter(size = 0.5, alpha = 0.5) +
  # and labels
  xlab('11:214169:T:C effect in condition and accessibility') + 
  ylab('PSMD13 expression') + 
  # cuter plots
  theme(legend.position = 'none') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = get_color_list(sc_exp_acc4_plottable[['condition_openness_genotype']]))

