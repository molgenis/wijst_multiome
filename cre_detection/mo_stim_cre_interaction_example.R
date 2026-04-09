#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_stim_cre_interaction_example.R
# Function: get the interaction examples where stim status has an impact on CRE-gene relationship
############################################################################################################################


####################
# libraries        #
####################

library(ggplot2)
# for plink file loading
library(snpStats)
# for regression models
library(lme4)
# Seurat and Signac
library(Seurat)
library(Signac)
# data table
library(data.table)


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


normalize_mj <- function(seurat_object) {
  # get the count matrix where we have the correct cell type
  count_matrix <- GetAssayData(seurat_object, slot = "counts")
  # ignore genes that are never expressed
  count_matrix <-  count_matrix[which(rowSums(count_matrix) != 0), ]
  # create new object to store the counts in
  norm_count_matrix <- count_matrix
  # do mean sample-sum normalization
  sample_sum_info = colSums(norm_count_matrix)
  mean_sample_sum = mean(sample_sum_info)
  sample_scale = sample_sum_info / mean_sample_sum
  # divide each column by sample_scale
  norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
  if ('layers' %in% slotNames(seurat_object[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    seurat_object[['MJ']] <- CreateAssay5Object(data = norm_count_matrix)
    
  } else {
    print('using Seurat v3/4 style \'slot\'')
    seurat_object[['MJ']] <- CreateAssayObject(data = norm_count_matrix)
  }
  return(seurat_object)
}


binarize_matrix <- function(seurat_object, input_assay='peaks', output_assay='peaks_bin') {
  # add new assay
  seurat_object[[output_assay]] <- GetAssay(seurat_object, input_assay)
  # binarize
  BinarizeCounts(seurat_object, assay = output_assay)
  return(seurat_object)
}
  


####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- T


###################
# interested data #
###################

# get specific variant, region and gene
gene <- 'LY86'
region <- 'chr6-6587847-6589226'
variant <- '6:6586961:A:G'
cell_type <- 'B'

#########################
# Main code bulk method #
#########################

# location of the expression data
exp_data_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/input/L1/combined/'
acc_data_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/input/L1/'
# and the sample mapping files
eqtl_smf_loc <- paste(exp_data_loc, '/smf.txt', sep = '/')
caqtl_smf_loc <- paste(acc_data_loc, '/smf.txt', sep = '/')
# read the sample mapping files
eqtl_smf <- read.table(eqtl_smf_loc, header = T, sep = '\t')
caqtl_smf <- read.table(caqtl_smf_loc, header = T, sep = '\t')
# read the expression data
eqtl_inputs = read_expression_files_per_celltype(exp_data_loc)
caqtl_inputs = read_expression_files_per_celltype(acc_data_loc)

# get the location of the annotations for stimulation status
stim_status_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_sample_to_stim.tsv.gz'
# read the stimulation status
stim_status <- fread(stim_status_loc, header = T, sep = '\t')
# add the sample identier in the interaction format to the status table
stim_status[['sample_lane']] <- paste(stim_status[['donor']], stim_status[['lane']], sep = ';;')

# the variant data
genotype_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/imputed_hg38_all_anc_qtl_tested_variants'
# read the genotypes, but only those in the file and in the confinement
genotypes <- read.plink(
  bed = paste(genotype_loc, '.bed', sep = ''),
  bim = paste(genotype_loc, '.bim', sep = ''),
  fam = paste(genotype_loc, '.fam', sep = '')
)
# get eQTL
eqtl_data_gene <- create_qtl_plot_table(feature=gene, variant=variant, feature_table=eqtl_inputs[[cell_type]], genotypes=genotypes, eqtl_smf)
# get caQTL
caqtl_data_region <- create_qtl_plot_table(feature=region, variant=variant, feature_table=caqtl_inputs[[cell_type]], genotypes=genotypes, caqtl_smf)
# set names
colnames(eqtl_data_gene) <- c('phenotype_id', 'gene', 'genotype_id', 'genotype', 'alleles')
colnames(caqtl_data_region) <- c('phenotype_id', 'region', 'genotype_id', 'genotype', 'alleles')
# merge both
qtl_data_region_gene <- merge(eqtl_data_gene, caqtl_data_region[, c('phenotype_id', 'region')], by = 'phenotype_id')
# add stim status
qtl_data_region_gene[['stim_status']] <- stim_status[match(qtl_data_region_gene[['phenotype_id']], stim_status[['sample_lane']]), ][['condition']]
# create a formula
stim_cre_base_formula <- as.formula(
  'gene ~ region + stim_status + (1|genotype_id)'
)
# and for the interaction
stim_cre_interaction_formula <- as.formula(
  'gene ~ region + stim_status + (1|genotype_id) + region * stim_status'
)
# do the two models
stim_cre_base_model <- lmerTest::lmer(formula = stim_cre_base_formula, data = qtl_data_region_gene)
stim_cre_interaction_model <- lmerTest::lmer(formula = stim_cre_interaction_formula, data = qtl_data_region_gene)
# check if they are different
ftest_res <- anova(stim_cre_base_model, stim_cre_interaction_model, refit = FALSE, test = 'F')


################################
# Main code single-cell method #
################################

# location of the Seurat object
seurat_object_loc <- paste0('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_', tolower(cell_type), '_1_80_20240521.rds')
# read the object
seurat_object <- readRDS(seurat_object_loc)
# add MJ normalization
seurat_object <- normalize_mj(seurat_object)
# binarize atac
seurat_object <- binarize_matrix(seurat_object)
# make plot table
plot_table_sc <- data.frame(
  'gene' = as.vector(unlist(seurat_object@assays$MJ@data[gene, ])), 
  'region' = as.vector(unlist(seurat_object@assays$peaks_bin@counts[region, ])), 
  'phenotype_id' = rownames(seurat_object@meta.data)
)
# add the genotype id
plot_table_sc[['genotype_id']] <- seurat_object@meta.data[match(plot_table_sc[['phenotype_id']], seurat_object@meta.data[['barcode_lane']]), ][['sample_final']]
# add the lane
plot_table_sc[['lane']] <- seurat_object@meta.data[match(plot_table_sc[['phenotype_id']], seurat_object@meta.data[['barcode_lane']]), ][['lane']]
# add the stimulation status
plot_table_sc[['stim_status']] <- seurat_object@meta.data[match(plot_table_sc[['phenotype_id']], seurat_object@meta.data[['barcode_lane']]), ][['condition_final']]
# create a formula
stim_cre_base_formula_sc <- as.formula(
  'gene ~ region + stim_status + (1|genotype_id) + (1|lane)'
)
# and for the interaction
stim_cre_interaction_formula_sc <- as.formula(
  'gene ~ region + stim_status + (1|genotype_id) + (1|lane) + region * stim_status'
)
# do the two models
stim_cre_base_model_sc <- lmerTest::lmer(formula = stim_cre_base_formula_sc, data = plot_table_sc)
stim_cre_interaction_model_sc <- lmerTest::lmer(formula = stim_cre_interaction_formula_sc, data = plot_table_sc)
# check if they are different
ftest_res_sc <- anova(stim_cre_base_model_sc, stim_cre_interaction_model_sc, refit = FALSE, test = 'F')
