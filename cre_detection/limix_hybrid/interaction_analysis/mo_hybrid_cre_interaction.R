#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_hybrid_cre_interaction.R
# Function: 
# Example: 
# Rscript 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
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




####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- F


####################
# Main code        #
####################

# make command line options
option_list <- list(
  make_option(c("-i", "--in"), type="character", default=NULL, 
              help="input directory of chunks", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output directory", metavar="character"), 
  make_option(c("-c", "--confinement"), type="character", default=NULL, 
              help="confinement file of tf-region-gene triplets to test", metavar="character"),
  make_option(c("-s", "--smf_loc"), type="character", default=NULL, 
              help="sample mapping file", metavar="character"), 
  make_option(c("-e", "--expression_file"), type="character", default='expression.tsv.gz', 
              help="expression filename for chunk", metavar="character"), 
  make_option(c("-a", "--accessibility_file"), type="character", default='accessibility.tsv.gz', 
              help="accessibility filename for chunk", metavar="character"), 
  ake_option(c("-v", "--covariates_file"), type="character", default=NULL, 
             help="accessibility filename for chunk", metavar="character"), 
  make_option(c("-f", "--fixed_effects"), type="character", default=NULL,
              help="comman separated list of fixed effects to correct for [default= %default]", metavar="character"),
  make_option(c("-r", "--random_effects"), type="character", default=NULL,
              help="comma separated list of random effects to correct for [default= %default]", metavar="character"),
  make_option(c("-g", "--accessibility_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on accessibility before modelling [default: %default]"), 
  make_option(c("-y", "--expression_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on expression before modelling [default: %default]"), 
  make_option(c("-t", "--interaction_terms"), type="character", default=NULL,
              help="interaction to model", metavar='character'), 
  make_option(c("-g", "--genotype_loc"), type="character", default=NULL,
              help="genotype file in plink1 format, without extension", metavar="character")
)


# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
# location of the region-to-gene files
confinement_loc <- NULL
# input directory
in_dir <- NULL
# location of SMF
smf_loc <- NULL
# location of the output
output_loc <- NULL
# expression filename
expression_file <- NULL
# expression filename
accessibility_file <- NULL
# genotype location
genotype_loc <- NULL
# whether to gausnorm the expression data
expression_gausnorm <- T
# whether to gausnorm the accessibility/TF data
accessibility_gausnorm <- T
# fixed effects string
fixed_effects_string <- NULL
# random effects string
random_effects_string <- NULL
# interaction terms string
interaction_terms_string <- NULL
# covariates file
covariates_file <- NULL

if (debug) {
  # set all of the variables hardcoded for a testing debug run
  confinement_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/CD4T/mo_gt_tf_gene_confinement.tsv.gz'
  in_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD4T/chr9-99361653-100377592/'
  smf_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD4T/smf.tsv.gz'
  output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/CD4T/chr9-99361653-100377592/'
  expression_file <- 'expression.tsv.gz'
  accessibility_file <- 'accessibility.tsv.gz'
  genotype_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr9'
  covariates_file <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
  fixed_effects_string <- 'lane,region,genotype'
  random_effects_string <- 'sample_final'
  interaction_terms_string <- 'genotype*region'
  
  
} else {
  # grab from the actual parameters if doing a real run
  # TODO implement this
}

# make the full path to the expression data
full_exp_path <- NULL
# depending on if it is an absolute path, we do things differently
if (startsWith(expression_file, '/')) {
  full_exp_path <- expression_file
} else {
  full_exp_path <- paste(in_dir, expression, sep = '/')
}
# same for the accessibility/TF data
full_acc_path <- NULL
if (startsWith(expression_file, '/')) {
  full_acc_path <- accessibility_file
} else {
  full_acc_path <- paste(in_dir, accessibility_file, sep = '/')
}
# and for covariates
full_covariates_path <- NULL
if (!is.null(covariates_file) & !is.na(covariates_file) & startsWith(covariates_file, '/')) {
  full_covariates_path <- covariates_file
} else if (!is.null(covariates_file) & !is.na(covariates_file)){
  full_covariates_path <- paste(in_dir, covariates_file, sep = '/')
}
# read the confinement file
confinement <- fread(confinement_loc, header = T, sep = '\t')
# set harmonized column names to make it easier for ourselves
colnames(confinement) <- c('variant', 'region', 'gene')

# read the expression file
expression_data <- fread(expression_file, header = T, sep = '\t', check.names = F)
# read the TF/accessibility data
accessibility_data <- fread(accessibility_file, header = T, sep = '\t', check.names = F)
# also set the first column name so we can refer to it later
colnames(expression_data)[[1]] <- 'gene'
colnames(accessibility_data)[[1]] <- 'region'
# subset both sets
expression_data <- expression_data[expression_data[['gene']] %in% confinement[['gene']], ]
accessibility_data <- accessibility_data[accessibility_data[['region']] %in% confinement[['region']], ]
# check if we have any data left
if (nrow(expression_data) > 0) {
  # check if we have any data left
  if (nrow(accessibility_data) > 0) {
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
    # filter the confinement on the variants we have in the genotype data as well
    confinement <- confinement[confinement[['variant']] %in% overlapping_variants, ]
    # and the genes we have
    confinement <- confinement[confinement[['gene']] %in% expression_data[['gene']], ]
    # and regions or TFs we have
    confinement <- confinement[confinement[['region']] %in% accessibility_data[['region']], ]
    # check if we have any data left
    if (nrow(confinement) > 0) {
      # read smf
      smf <- fread(smf_loc, header = T, sep = '\t')
      # harmonize names
      colnames(smf) <- c('cell', 'participant')
      # intersect the smf with the expression data and accessibility/TF data
      intersecting_cells <- intersect(smf[['cell']], colnames(expression_data))
      intersecting_cells <- intersect(intersecting_cells, colnames(accessibility_data))
      
      # check if there is any data left
      if (length(intersecting_cells) > 0) {
        # subset the matrices
        smf <- smf[smf[['cell']] %in% intersecting_cells, ]
        accessibility_data <- accessibility_data[, c('region', intersecting_cells)]
        expression_data <- expression_data[, c('gene', intersecting_cells)]
        
        # split the fixed effects
        fixed_effects <- c()
        if (!is.null(fixed_effects_string) & !is.na(fixed_effects_string) & fixed_effects_string != '') {
          fixed_effects <- strsplit(fixed_effects_string, ',')[[1]]
        }
        # warn if we are not including the genotype
        if ('genotype' %in% fixed_effects) {
          warning(paste('\'genotype\' term not present in fixed effects!'))
        }
        # split random effects
        random_effects <- c()
        if (!is.null(random_effects_string) & !is.na(random_effects_string) & random_effects_string != '') {
          random_effects <- strsplit(random_effects_string, ',')[[1]]
        }
        # create formula
        base_formula <- get_formula(var_of_interest = 'expression', fixed_effects = fixed_effects, random_effects = random_effects)
        message(paste('Using base formula', base_formula))
        # split interactions
        interactions <- strsplit(interaction_terms_string, ',')[[1]]
        # paste together
        interactions_formula_part <- paste(interactions, collapse = '*')
        # create interaction formula
        interaction_formula <- paste(base_formula, interactions_formula_part, sep = '+')
        message(paste('Using interaction formula', interaction_formula))
        # if a covariate matrix was supplied, we'll load it
        if (!is.null(full_covariates_path)) {
          covariates_data <- fread(full_covariates_path, header = T, sep = '\t')
          # the first column should be the cell
          colnames(covariates_data)[[1]] <- 'cell'
          # intersect this
          intersecting_cells <- intersect(intersecting_cells, covariates_data[['cell']])
          # subset the covariates
          covariates_data <- covariates_data[covariates_data[['cell']] %in% intersecting_cells, ]
        }
        # check if we have any data left
        if (length(intersecting_cells) > 0) {
          # do gaussnorm if so requested
          if (expression_gausnorm) {
            message('Yeo-Johnson gausnorm on expression data...')
            expression_data <- gausnorm_independent_variable_matrix(independent_variable_matrix = expression_data, feature_id_column = 'gene')
          }
          if (accessibility_gausnorm) {
            message('Yeo-Johnson gausnorm on accessibility/TF data...')
            accessibility_data <- gausnorm_independent_variable_matrix(independent_variable_matrix = accessibility_data, feature_id_column = 'region')
          }
          
          
        } else {
          message('No cells left after intersecting with covariates matrix. No more work to be done')
        }
      } else {
        message('No cells left after intersecting smf with expression and accessibility/TF. No more work to be done')
      }
    } else {
      message('No triplets left after intersecting confinement with data. No more work to be done')
    }
  } else {
    message('No regions/TFs left after filtering confinement. No more work to be done')
  }
} else {
  message('No genes left after filtering confinement. No more work to be done')
}

