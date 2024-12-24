#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_perform_qtl_mediation_analysis.R
# Function: perform mediation analysis of eQTLs through caQTLs
############################################################################################################################

####################
# libraries        #
####################

# for tables, much faster than data.frame and matrix
library(data.table)
# for plink file loading
library(snpStats)
# for regression models
library(lme4)
# mediation dependencies for bootstrap methods
library(MEPS)
library(bda)
library(mediation)
# transformation into gaussian normal distribution
library(bestNormalize)
# for automatic md5 file creation
library(mdfiver)
# for multithreading
# library(parallel)
# library(doParallel)
# show progress
library(pbapply)
# parse command line arguments
library(optparse)

####################
# Functions        #
####################

#' Create a formula for mixed-effects models
#'
#' This function generates a formula for mixed-effects models based on the specified variable of interest, fixed effects, and random effects.
#'
#' @param var_of_interest A character string representing the dependent variable.
#' @param fixed_effects A character vector of fixed effect variables.
#' @param random_effects A character vector of random effect variables.
#' @return A formula object for use in mixed-effects models.
#' @examples
#' get_formula("y", c("x1", "x2"), c("group"))
get_formula <- function(var_of_interest, fixed_effects, random_effects) {
  # make the formula
  formula_string <- paste(var_of_interest, '~ 0 ', sep = ' ')
  # do the random effects
  if (!is.null(random_effects) & length(random_effects) > 0) {
    # style of R formulas
    random_effects_formula_style <- paste('(1|', random_effects, ')', sep = '')
    formula_string <- paste(formula_string, paste(random_effects_formula_style, collapse = '+'), sep = '+')
  }
  # do the fixed effects
  if (!is.null(fixed_effects) & length(fixed_effects) > 0) {
    formula_string <- paste(formula_string, paste(fixed_effects, collapse = '+'), sep = '+')
  }
  # turn into formula
  form <- as.formula(formula_string)
  return(form)
}


#' Create a mediation table for QTL analysis
#'
#' This function generates a mediation table for quantitative trait loci (QTL) analysis by combining accessibility, expression, and genotype data with metadata.
#'
#' @param accessibility_table A data.table containing accessibility data.
#' @param expression_table A data.table containing expression data.
#' @param metadata A data.table containing metadata with sample information.
#' @param genotypes A list containing genotype data.
#' @param variant A character string representing the variant of interest.
#' @param region A character string representing the region of interest.
#' @param gene A character string representing the gene of interest.
#' @return A data.frame containing the mediation table for QTL analysis.
#' @examples
#' create_mediation_table(accessibility_table, expression_table, metadata, genotypes, "rs12345", "chr1:1000-2000", "GENE1")
create_mediation_table <- function(accessibility_table, expression_table, metadata, genotypes, variant, region, gene) {
  # extract the expression
  expression <- expression_table[expression_table$feature == gene, .SD, .SDcols = metadata[['sample']]]
  # extract the accessibility
  accessibility <- accessibility_table[accessibility_table$feature == region, .SD, .SDcols = metadata[['sample']]]
  # extract genotypes
  genotype <- genotypes$genotypes[metadata[['donor']], variant]
  # then to numeric
  genotype_numeric <- as(genotype, 'numeric')
  
  # finally as vector
  genotype_vector <- as.vector(genotype_numeric[, 1])
  expression_vector <- as.vector(unlist(expression))
  accessibility_vector <- as.vector(unlist(accessibility))
  
  # make into qtl table
  qtl_table <- data.frame('genotype' = genotype_vector, expression = expression_vector, accessibility = accessibility_vector)
  # get the complete cases
  qtl_complete_cases <- complete.cases(qtl_table)
  # merge complete qtl and metadata table
  mediation_table <- cbind(qtl_table[qtl_complete_cases, ], metadata[qtl_complete_cases, ])
  
  return(mediation_table)
}


#' Perform mediation analysis
#'
#' This function conducts a mediation analysis using either linear models or mixed-effects models, depending on the presence of random effects in the formula.
#'
#' @param mediation_table A data.frame containing the mediation data.
#' @param formula_indirect1 A formula for the first part of the indirect effect model.
#' @param formula_mediation A formula for the direct and indirect effect model.
#' @param sims An integer specifying the number of simulations for the mediation analysis. Default is 1000.
#' @return An object of class \code{mediate} containing the results of the mediation analysis.
#' @examples
#' perform_mediation_analysis(mediation_table, formula_indirect1, formula_mediation, sims = 1000)
perform_mediation_analysis <- function(mediation_table, formula_indirect1, formula_mediation, sims=1000) {
  # set result
  mediation_results <- NULL
  # build random
  if (grepl('\\|', Reduce(paste, deparse(formula_indirect1)))) {
    ## First-part indirect effect model:
    indirect_model1 <- do.call(what = 'lmer', list(formula = formula_indirect1, data = mediation_table))
    ## Second-part direct + indirect effect model:
    indirect_model2 <- do.call(what = 'lmer', list(formula = formula_mediation, data = mediation_table))
    # warn about boostrap
    warning(paste('cannot currently boostrap with random effects'))
    ## Mediation analysis with 1000 simulations
    mediation_results <- mediate(indirect_model1, indirect_model2, treat = 'genotype', mediator = 'accessibility', boot = F, sims = sims)
  }
  # or only linear
  else {
    ## First-part indirect effect model:
    indirect_model1 <- do.call(what = 'glm', list(formula = formula_indirect1, data = mediation_table, family = gaussian(link = "identity")))
    ## Second-part direct + indirect effect model:
    indirect_model2 <- do.call(what = 'glm', list(formula = formula_mediation, data = mediation_table, family = gaussian(link = "identity")))
    ## Mediation analysis with 1000 simulations
    mediation_results <- mediate(indirect_model1, indirect_model2, treat = 'genotype', mediator = 'accessibility', boot = T, sims = sims)
  }
  return(mediation_results)
}


#' Mediate all effects for a set of variants, regions, and genes
#'
#' This function performs mediation analysis for a set of variants, regions, and genes by combining accessibility, expression, and genotype data with metadata.
#'
#' @param accessibility A data.table containing accessibility data.
#' @param expression A data.table containing expression data.
#' @param metadata A data.table containing metadata with sample information.
#' @param genotypes A list containing genotype data.
#' @param confinement A data.table or matrix where each row specifies a variant, region, and gene.
#' @param form_indirect1 A formula for the first part of the indirect effect model.
#' @param form_indirect2 A formula for the direct and indirect effect model.
#' @param cluster An optional parameter for parallel processing (default is NULL).
#' @return A list of mediation analysis results for each set of variant, region, and gene.
#' @examples
#' mediate_all_effects(accessibility, expression, metadata, genotypes, confinement, form_indirect1, form_indirect2)
mediate_all_effects <- function(accessibility, expression, metadata, genotypes, confinement, form_indirect1, form_indirect2, cluster=NULL) {
  # save result per set
  res_per_set <- pbapply::pbapply(confinement, 1, function(x) {
    # get from the row
    variant <- x[1][[1]]
    region <- x[2][[1]]
    gene <- x[3][[1]]
    # get the mediation table
    medation_table <- create_mediation_table(accessibility_table = accessibility, expression_table = expression, metadata = metadata, genotypes = genotypes, variant = variant, region = region, gene = gene)
    # do the mediation analysis
    med_result <- perform_mediation_analysis(mediation_table = medation_table, formula_indirect1 = form_indirect1, formula_mediation = form_indirect2)
    # add extra information
    med_result[['variant']] <- variant
    med_result[['region']] <- region
    med_result[['gene']] <- gene
    return(med_result)
  })
  return(res_per_set)
}


medatiate_mer_to_table <- function(mediation_mer_object) {
  # collect effects
  coef_indirect <- mediation_mer_object$d1
  coef_direct <- mediation_mer_object$z1
  total_effect <- mediation_mer_object$tau.coef
  prop_med <- mediation_mer_object$n1
  # and p values
  p_indirect <- mediation_mer_object$d1.p
  p_direct <- mediation_mer_object$z1.p
  p_total <- mediation_mer_object$tau.p
  p_mediated <- mediation_mer_object$n1.p
  
  # get indirect formulas
  formula_ind1 <- Reduce(paste, deparse(mediation_mer_object$model.m$formula))
  formula_ind2 <- Reduce(paste, deparse(mediation_mer_object$model.y$formula))
  # remove whitespace to make eventual tsv smaller
  formula_ind1 <- gsub(' +', '', formula_ind1)
  formula_ind2 <- gsub(' +', '', formula_ind2)
  
  # turn into a single-row table
  row_mediation <- data.frame(
    'variant' = c(mediation_mer_object[['variant']]), 
    'accessibility' = c(mediation_mer_object[['region']]), 
    'gene' = c(mediation_mer_object[['gene']]), 
    'mediator' = c(mediation_mer_object$mediator), 
    'p_mediated' = c(p_mediated), 
    'prop_med' = c(prop_med), 
    'p_total' = c(p_total), 
    'total_effect' = c(total_effect), 
    'p_direct' = c(p_direct), 
    'coef_direct' = c(coef_direct), 
    'p_indirect' = c(p_indirect), 
    'coef_indirect' = c(coef_indirect), 
    'form_indirect1' = c(formula_ind1), 
    'form_indirect2' = c(formula_ind2)
  )
  return(row_mediation)
}


#' Convert mediation analysis results to a table
#'
#' This function converts the results of a mediation analysis into a single-row data.frame for easy export and further analysis.
#'
#' @param mediation_mer_object An object containing the results of the mediation analysis.
#' @return A data.frame containing the mediation analysis results, including coefficients, p-values, and formulas.
#' @examples
#' medatiate_mer_to_table(mediation_mer_object)
medatiate_to_table <- function(mediation_object) {
  # collect effects
  coef_indirect <- mediation_object$d1
  coef_direct <- mediation_object$z1
  total_effect <- mediation_object$tau.coef
  prop_med <- mediation_object$n1
  # and p values
  p_indirect <- mediation_object$d1.p
  p_direct <- mediation_object$z1.p
  p_total <- mediation_object$tau.p
  p_mediated <- mediation_object$n1.p
  
  # get indirect formulas
  formula_ind1 <- Reduce(paste, deparse(summary(mediation_object$model.m)$call$formula))
  formula_ind2 <- Reduce(paste, deparse(summary(mediation_object$model.y)$call$formula))
  # remove whitespace to make eventual tsv smaller
  formula_ind1 <- gsub(' +', '', formula_ind1)
  formula_ind2 <- gsub(' +', '', formula_ind2)
  
  # turn into a single-row table
  row_mediation <- data.frame(
    'variant' = c(mediation_object[['variant']]), 
    'accessibility' = c(mediation_object[['region']]), 
    'gene' = c(mediation_object[['gene']]), 
    'mediator' = c(mediation_object$mediator), 
    'p_mediated' = c(p_mediated), 
    'prop_med' = c(prop_med), 
    'p_total' = c(p_total), 
    'total_effect' = c(total_effect), 
    'p_direct' = c(p_direct), 
    'coef_direct' = c(coef_direct), 
    'p_indirect' = c(p_indirect), 
    'coef_indirect' = c(coef_indirect), 
    'form_indirect1' = c(formula_ind1), 
    'form_indirect2' = c(formula_ind2)
  )
  return(row_mediation)
}


#' Convert a list of mediation results to a single table
#'
#' This function processes a list of mediation analysis results and converts them into a single data.frame for easy export and further analysis.
#'
#' @param mediation_results_list A list of mediation analysis results.
#' @return A data.frame containing the combined mediation analysis results.
#' @examples
#' mediation_to_tables(mediation_results_list)
mediation_to_tables <- function(mediation_results_list) {
  # create a new list to store the rows of the converted objects
  mediation_table_list <- list()
  # check each mediation result
  for (i in 1: length(mediation_results_list)) {
    # get the mediation effect
    med_effect <- mediation_results_list[[i]]

    # init variable
    row_mediation <- NULL
    # depending on the class, we use a different function
    if (attributes(med_effect)$class == 'mediate.mer') {
      row_mediation <- medatiate_to_table(med_effect)
    }
    else if(attributes(med_effect)$class == 'mediate') {
      row_mediation <- medatiate_mer_to_table(med_effect)
    }
    
    # put into list
    mediation_table_list[[i]] <- row_mediation
  }
  # make into one table
  mediation_table <- do.call('rbind', mediation_table_list)
  return(mediation_table)
}


#' Transform Independent Variable Matrix Using Yeo-Johnson Transformation
#'
#' This function applies the Yeo-Johnson transformation to the numeric columns of an independent variable matrix.
#' The feature ID column is preserved and reattached to the transformed data.
#'
#' @param independent_variable_matrix A data.table containing the independent variables. The columns represent different donors.
#' @param feature_id_column A character string specifying the column name that contains the feature IDs. Default is 'feature'.
#'
#' @return A data.table with the transformed numeric columns and the feature ID column reattached.
#'
#' @examples
#' \dontrun{
#' library(data.table)
#' library(car)
#' dt <- data.table(feature = c('A', 'B', 'C'), donor1 = c(1, 2, 3), donor2 = c(4, 5, 6))
#' transformed_dt <- gausnorm_independent_variable_matrix(dt, 'feature')
#' print(transformed_dt)
#' }
#'
gausnorm_independent_variable_matrix <- function(independent_variable_matrix, feature_id_column='feature') {
  # take the features
  features <- independent_variable_matrix[[feature_id_column]]
  # remove the feature ID
  independent_variable_matrix[[feature_id_column]] <- NULL
  # take the donor names, as they are the columns
  colnames_original <- colnames(independent_variable_matrix)
  # transpose the matrix, as we'll do this on a per-column basis
  transformed_data <- as.data.table(
    lapply(independent_variable_matrix, function(x) {
      if (is.numeric(x)) {
        yeojohnson(x)$x.t
      }
      else {
        x
      }
    })
  )
  # add back the donor names
  colnames(transformed_data) <- colnames_original
  # make the features as a data.table as well
  features_column <- data.table(x = features)
  # with the right column name
  colnames(features_column) <- feature_id_column
  # and merge the feature column back onto the data
  transformed_data <- cbind(features_column, transformed_data)
  return(transformed_data)
}


#' Filter input data for mediation analysis
#'
#' This function filters the input data to ensure that only samples with complete data across expression, accessibility, and genotype datasets are retained.
#'
#' @param inputs A list containing the input data: expression, accessibility, metadata, and genotypes.
#' @return A list with the filtered input data.
#' @examples
#' filtered_inputs <- filter_inputs(inputs)
filter_inputs <- function(inputs) {
  # load expression
  expression <- inputs[['expression']]
  # load accessibility
  accessibility <- inputs[['accessibility']]
  # load metadata
  metadata <- inputs[['metadata']]
  # load genotype data
  genotypes_data <- inputs[['genotypes']]
  # get the table of sample to genotype
  sample_to_donor <- metadata[, c('sample', 'donor')]
  # get whom we have the genotypes for
  samples_genotypes <- intersect(rownames(genotypes_data$genotypes), sample_to_donor[['donor']])
  # subset
  sample_to_donor <- sample_to_donor[sample_to_donor[['donor']] %in% samples_genotypes, ]
  # get whom we have expression data for
  samples_expression <- intersect(colnames(expression), sample_to_donor[['sample']])
  # subset
  sample_to_donor <- sample_to_donor[sample_to_donor[['sample']] %in% samples_expression, ]
  # get whom we have the atac data for
  samples_accessibility <- intersect(colnames(accessibility), sample_to_donor[['sample']])
  # subset
  sample_to_donor <- sample_to_donor[sample_to_donor[['sample']] %in% samples_accessibility, ]
  
  # subset genotypes
  genotypes_data$genotypes <- genotypes_data$genotypes[sample_to_donor[['donor']], ]
  genotypes_data$fam <-  genotypes_data$fam[rownames(genotypes_data$genotypes), ]
  
  # filter expression
  expression <- expression[, .SD, .SDcols=c('feature', sample_to_donor[['sample']])]
  # filter accessibility
  accessibility <- accessibility[, .SD, .SDcols=c('feature', sample_to_donor[['sample']])]
  # and finally metadata
  metadata <- metadata[metadata[['sample']] %in% sample_to_donor[['sample']], ]
  
  # do Yeo-Johnson transformation to get a gaussian distribution per feature, if requested
  if (inputs[['caqtl_gausnorm']]) {
    message('performing Yeo-Johnson transformation on accessibility data')
    accessibility <- gausnorm_independent_variable_matrix(accessibility)
  }
  if (inputs[['eqtl_gausnorm']]) {
    message('performing Yeo-Johnson transformation on expression data')
    expression <- gausnorm_independent_variable_matrix(expression)
  }

  # put back into the list
  inputs[['expression']] <- expression
  inputs[['accessibility']] <- accessibility
  inputs[['genotypes']] <- genotypes_data
  inputs[['metadata']] <- metadata
  
  return(inputs)
}


#' Load and preprocess input data for mediation analysis
#'
#' This function reads and preprocesses the input data files required for mediation analysis, including expression, accessibility, metadata, and genotype data.
#'
#' @param options A list of options specifying the file paths and other parameters.
#' @return A list containing the preprocessed input data.
#' @examples
#' options <- list(
#'   eqtl_file = "path/to/eqtl_file.txt",
#'   caqtl_file = "path/to/caqtl_file.txt",
#'   confinement_list = "path/to/confinement_list.txt",
#'   genotype_file = "path/to/genotype_file",
#'   metadata = "path/to/metadata.txt",
#'   fixed_effects = "effect1,effect2",
#'   random_effects = "effect3,effect4",
#'   out = "output_directory"
#' )
#' inputs <- load_inputs(options)
load_inputs <- function(options) {
  # read the eqtl file
  expression <- read.table(options[['eqtl_file']], header = T, sep = '\t', check.names = F, row.names = 1)
  # set the feature name
  expression <- cbind(data.frame('feature' = rownames(expression)), expression)
  # and turn into data.table
  expression <- data.table::data.table(expression)
  
  # read the caqtl file
  accessibility <- read.table(options[['caqtl_file']], header = T, sep = '\t', check.names = F, row.names = 1)
  # set the feature name
  accessibility <- cbind(data.frame('feature' = rownames(accessibility)), accessibility)
  # and turn into data.table
  accessibility <- data.table::data.table(accessibility)
  
  # read the confinement file
  confinement <- data.table::fread(options[['confinement_list']], header = F, sep = '\t', check.names = F)

  # subset the expression and accessibility
  accessibility <- accessibility[accessibility[['feature']] %in% confinement[[2]], ]
  expression <- expression[expression[['feature']] %in% confinement[[3]], ]
  
  # read the bim
  variants_in_gt <- fread(paste(options[['genotype_file']], '.bim', sep = ''), header = F)[[2]]
  
  # get overlapping variants
  overlapping_variants <- intersect(confinement[[1]], variants_in_gt)
  # read the genotypes, but only those in the file and in the confinement
  genotypes <- read.plink(
    bed = paste(options[['genotype_file']], '.bed', sep = ''),
    bim = paste(options[['genotype_file']], '.bim', sep = ''),
    fam = paste(options[['genotype_file']], '.fam', sep = ''), 
    select.snps = overlapping_variants
  )
  # filter the confinement on the variants we have in the genotype data as well
  confinement <- confinement[confinement[[1]] %in% overlapping_variants, ]
  
  # read the metadata
  metadata <- fread(options[['metadata']], header = T, sep  = '\t')
  
  # get the fixed and random effects
  fixed_effects <- c()
  if (!is.null(options[['fixed_effects']])) {
    fixed_effects <- strsplit(options[['fixed_effects']], split = ',')[[1]]
  }
  random_effects <- c()
  if (!is.null(options[['random_effects']])) {
    random_effects <- strsplit(options[['random_effects']], split = ',')[[1]]
  }
  # create the formula
  #form_indirect1 <- get_formula('accessibility', c(fixed_effects, 'genotype'), random_effects)
  form_indirect1 <- get_formula('accessibility', c('genotype'), random_effects)
  form_indirect2 <- get_formula('expression', c('accessibility', fixed_effects, 'genotype'), random_effects)
  # return a list with each of the inputs
  inputs <- list(
    'expression' = expression, 
    'accessibility' = accessibility, 
    'metadata' = metadata, 
    'confinement' = confinement, 
    'genotypes' = genotypes, 
    'form_indirect1' = form_indirect1, 
    'form_indirect2' = form_indirect2, 
    'out' = options[['out']], 
    'caqtl_gausnorm' = options[['caqtl_gausnorm']], 
    'eqtl_gausnorm' = options[['eqtl_gausnorm']]
  )
  return(inputs)
}


#' Run the full mediation analysis pipeline
#'
#' This function runs the entire mediation analysis pipeline, including loading inputs, filtering data, performing mediation analysis, summarizing results, and writing the output.
#'
#' @param options A list of options specifying the file paths and other parameters.
#' @param verbose A logical value indicating whether to print progress messages. Default is TRUE.
#' @return An integer value indicating the success of the function (0 for success).
#' @examples
#' options <- list(
#'   eqtl_file = "path/to/eqtl_file.txt",
#'   caqtl_file = "path/to/caqtl_file.txt",
#'   confinement_list = "path/to/confinement_list.txt",
#'   genotype_file = "path/to/genotype_file",
#'   metadata = "path/to/metadata.txt",
#'   fixed_effects = "effect1,effect2",
#'   random_effects = "effect3,effect4",
#'   out = "output_directory",
#'   threads = 4
#' )
#' run_full_analysis(options, verbose = TRUE)
run_full_analysis <- function(options, verbose=T) {
  # start multicore
  # cl <- makeCluster(options[['threads']])
  # registerDoParallel(cl)
  # get the inputs
  if (verbose) {
    message('loading inputs...')
  }
  inputs <- load_inputs(options)
  # filter the inputs
  if (verbose) {
    message('filtering inputs...')
  }
  inputs <- filter_inputs(inputs)
  # run each mediation
  if (verbose) {
    message('running mediation...')
  }
  mediation_results <- mediate_all_effects(
    accessibility = inputs$accessibility, 
    expression = inputs$expression, 
    metadata = inputs$metadata, 
    genotypes = inputs$genotypes, 
    confinement = inputs$confinement, 
    form_indirect1 = inputs$form_indirect1, 
    form_indirect2 = inputs$form_indirect2, 
    # cluster = cl)
    cluster = NULL)
  # summarize results
  if (verbose) {
    message('summarizing results...')
  }
  mediation_table <- mediation_to_tables(mediation_results)
  # write results
  if (verbose) {
    message('writing results...')
  }
  output_loc <- options[['out']]
  # gz file ends with .gz
  if (grepl('.gz$', output_loc)) {
    # gzip if ends with .gz
    output_loc <- gzfile(output_loc)
  }
  write.table(mediation_table, output_loc, sep = '\t', row.names = F, col.names = T, quote = F)
  # also make a checksum
  mdfiver::create_md5_for_file(options[['out']])
  # say where done
  if (verbose) {
    message('done')
  }
  return(0)
}


#' Run a debug version of the full mediation analysis pipeline
#'
#' This function sets up and runs a debug version of the full mediation analysis pipeline with predefined options.
#'
#' @return An integer value indicating the success of the function (0 for success).
#' @examples
#' do_debug()
do_debug <- function() {
  options_debug <- list()
  options_debug[['eqtl_file']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/L1/UT/monocyte.qtlInput.txt.gz'
  options_debug[['caqtl_file']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/monocyte.qtlInput.txt.gz'
  options_debug[['genotype_file']] <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr7'
  options_debug[['confinement_list']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/confinements/monocyte.confinement.tsv.gz'
  options_debug[['metadata']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/metadata/monocyte.metadata.tsv.gz'
  options_debug[['fixed_effects']] <- 'RNA_UT_PC1,RNA_UT_PC2,RNA_UT_PC3,RNA_UT_PC4,RNA_UT_PC5,RNA_UT_PC6,RNA_UT_PC7,RNA_UT_PC8,RNA_UT_PC9,RNA_UT_PC10'
  options_debug[['random_effects']] <- 'donor'
  #options_debug[['random_effects']] <- ''
  options_debug[['out']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/output/atac_to_expression/UT/monocyte.chr7.tsv.gz'
  options_debug[['threads']] <- 4
  options_debug[['caqtl_gausnorm']] <- F
  options_debug[['eqtl_gausnorm']] <- T
  run_full_analysis(options_debug)
}


####################
# Settings        #
####################

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)

####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-e", "--eqtl_file"), type="character", default=NULL,
              help="tsv file of the aggregated expression", metavar="character"),
  make_option(c("-c", "--caqtl_file"), type="character", default=NULL,
              help="tsv file of the aggregated accessibility [default= %default]", metavar="character"),
  make_option(c("-g", "--genotype_file"), type="character", default=NULL,
              help="genotype file in plink1 format, without extension [default= %default]", metavar="character"),
  make_option(c("-l", "--confinement_list"), type="character", default=NULL,
              help="tsv with in order variant-atacregion-gene to test for interactions [default= %default]", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="tsv file of the (aggregated) metadata [default= %default]", metavar="character"),
  make_option(c("-f", "--fixed_effects"), type="character", default=NULL,
              help="comman separated list of fixed effects to correct for [default= %default]", metavar="character"),
  make_option(c("-r", "--random_effects"), type="character", default=NULL,
              help="comman separated list of random effects to correct for [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output location of analysis [default= %default]", metavar="character"), 
  make_option(c("-a", "--caqtl_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on aggregated accessibility before modelling [default: %default]"), 
  make_option(c("-q", "--eqtl_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on aggregated expression before modelling [default: %default]")
  # make_option(c("-t", "--threads"), type="numeric", default=1,
  #             help="number of threads to use [default= %default]", metavar="numeric")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# there are some things we cannot allow
if (is.null(opt[['eqtl_file']])) {
  stop('eQTL expression matrix file must be supplied')
}
if (is.null(opt[['caqtl_file']])) {
  stop('caQTL expression matrix file must be supplied')
}
if (is.null(opt[['genotype_file']])) {
  stop('genotype file must be supplied')
}
if (is.null(opt[['confinement_list']])) {
  stop('confinement file must be supplied')
}
if (is.null(opt[['metadata']])) {
  stop('metadata file must be supplied')
}
if (is.null(opt[['out']])) {
  stop('output file must be supplied')
}

# do the pipeline
run_full_analysis(opt)
