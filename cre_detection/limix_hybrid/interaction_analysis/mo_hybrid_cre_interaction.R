#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_hybrid_cre_interaction.R
# Function: perform interaction-eQTL analysis at single-cell level with TF or ATAC as interaction terms
# Example: 
# ~/start_Rscript.sh \
#   /groups/umcg-franke-scrna/tmp02/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_hybrid_cre_interaction.R \
#   --in /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD4T/chr7-159197098-159254288 \
#   --out /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/CD4T/chr7-159197098-159254288 \
#   --confinement /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement.tsv.gz \
#   --smf_loc /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD4T/smf.tsv.gz \
#   --covariates_file /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz \
#   --fixed_effects region,genotype \
#   --random_effects sample_final,lane \
#   --interaction_terms genotype,region \
#   --barcode_column barcode_lane \
#   --genotype_loc /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr7 \
#   --expression_gausnorm \
#   --accessibility_gausnorm
# 
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)
# load command line parameters
library(optparse)
# use plink files
library(snpStats)
# transformation into gaussian normal distribution
library(bestNormalize)
# for the model
library(lme4)
library(lmerTest)


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
get_formula <- function(var_of_interest, fixed_effects=NULL, random_effects=NULL) {
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


model_to_row <- function(model) {
  # make summary of model
  model_summary <- summary(model)
  # extract coeficients
  model_coefficients <- model_summary[['coefficients']]
  # get the values we have metrics for
  covariates <- rownames(model_coefficients)
  # create rename dictionary
  rename_dict <- list(
    '^Estimate$' = 'beta', 
    '^Std. Error$' = 'se', 
    '^t value$' = 'tval', 
    '^Pr\\(>\\|t\\|\\)$' = 'p'
  )
  # rename all
  for (original in names(rename_dict)) {
    colnames(model_coefficients) <- gsub(original, rename_dict[[original]], colnames(model_coefficients))
  }
  # get the type of values
  stats <- colnames(model_coefficients)
  # create a df of one row, with columns that are a combination of covariates and their stats
  row_created <- data.frame(matrix(, nrow = 1, ncol = length(covariates) * length(covariates)))
  # set index
  i <- 1
  # check each variable
  for (covariate in covariates) {
    # and the value
    for (stat in stats) {
      # extract value
      row_created[1, i] <- model_coefficients[covariate, stat]
      # update column name
      colnames(row_created)[[i]] <- paste(covariate, stat, sep = '_')
      # update index
      i <- i + 1
    }
  }
  return(row_created)
}


do_interaction_analysis <- function(expression_data, 
                                    accessibility_data, 
                                    genotype_data, 
                                    smf, 
                                    confinement,
                                    covariates_data=NULL, 
                                    fixed_effects=c('lane','region','genotype'), 
                                    random_effects=c('sample_final'), 
                                    interactions=c('genotype','region'), 
                                    family = 'gaussian') {
  # create formula
  base_formula <- get_formula(var_of_interest = 'expression', fixed_effects = fixed_effects, random_effects = random_effects)
  base_formula_string <- (Reduce(paste, deparse(base_formula)))
  message(paste('Using base formula:', base_formula_string, '\n'))
  # paste together
  interactions_formula_part <- paste(interactions, collapse = '*')
  # create interaction formula
  interaction_formula_string <- paste(base_formula_string, interactions_formula_part, sep = ' + ')
  message(paste('Using interaction formula:', interaction_formula_string, '\n'))
  interaction_formula <- as.formula(interaction_formula_string)
  # conver to dataframes
  expression_data <- data.frame(expression_data)
  accessibility_data <- data.frame(accessibility_data)
  smf <- data.frame(smf)
  confinement <- data.frame(confinement)
  if (!is.null(covariates_data)) {
    covariates_data <- data.frame(covariates_data)
    # get which things we need from that dataframe
    covariate_columns <- setdiff(c(fixed_effects, random_effects, interactions), c('expression', 'region', 'genotype'))
    # and subset the covariates data to that
    covariates_data <- covariates_data[, c('cell', covariate_columns)]
  }
  # put all results in a list
  res_per_comparison <- list()
  # check each region
  for (region in unique(accessibility_data[['region']])) {
    # extract the genes and variants
    confinement_region <- confinement[confinement[['region']] == region, ]
    # the specific genes then
    genes_region <- unique(confinement_region[['gene']])
    # subset to these genes
    expression_data_region <- expression_data[
      expression_data[['gene']] %in% genes_region, 
    ]
    # extract the region values
    region_values <- as.vector(unlist(accessibility_data[accessibility_data[['region']] == region, 2:ncol(accessibility_data)]))
    # check each gene
    for (gene in unique(genes_region)) {
      # extract the gene
      gene_values <- as.vector(unlist(expression_data_region[expression_data_region[['gene']] == gene, 2:ncol(expression_data_region)]))
      # merge the metadata with the gene and the region
      covariates_data[['region']] <- region_values
      covariates_data[['expression']] <- gene_values
      # get the variants for this region-gene combination
      variants_region_gene <- unique(confinement_region[confinement_region[['gene']] == gene, ][['variant']])
      # check each variant
      for (variant in variants_region_gene) {
        # extract genotypes
        genotype <- genotype_data$genotypes[smf[['participant']], variant]
        # then to numeric
        genotype_numeric <- as.vector(as(genotype, 'numeric'))
        # add the genotype
        covariates_data[['genotype']] <- genotype_numeric
        # keep only complete cases
        covariates_data_complete <- covariates_data[complete.cases(covariates_data), ]
        # initialize variables
        base_model <- NULL
        interaction_model <- NULL
        ftest_res <- NULL
        anova_test_used <- NULL
        # depending on the family, the calls and anovas are different
        if (family == 'poisson') {
          # model without interaction
          base_model <- lme4::glmer(formula = base_formula, data = covariates_data_complete, family = poisson)
          # model with interaction
          interaction_model <- lme4::glmer(formula = interaction_formula, data = covariates_data_complete, family = poisson)
          # check if they are different
          ftest_res <- anova(base_model, interaction_model, refit = FALSE, test = 'Chisq')
          anova_test_used <- 'LRT'
        }
        else if (family == 'gaussian') {
          # base model
          base_model <- lmerTest::lmer(formula = base_formula, data = covariates_data_complete)
          # interaction model
          interaction_model <- lmerTest::lmer(formula = interaction_formula, data = covariates_data_complete)
          # check if they are different
          ftest_res <- anova(base_model, interaction_model, refit = FALSE, test = 'F')
          anova_test_used <- 'F'
        }
        else {
          stop(paste0('unknown family ', family, ', only valid families are gaussian and poisson'))
        }
        # convert the model to a df
        interaction_model_df <- model_to_row(interaction_model)
        # add the family used
        interaction_model_df[['family']] <- family
        # add the interaction test used
        interaction_model_df[['anova']] <- anova_test_used
        # and finally the value (which is the last column, but the name will differ based on the type of test used)
        interaction_model_df[['anova_p']] <- ftest_res['interaction_model', ncol(ftest_res)]
        # keep the number of cells we have
        interaction_model_df[['ncell']] <- nrow(covariates_data_complete)
        # and participants
        interaction_model_df[['nparticipant']] <- length(unique(smf[smf[['cell']] %in% covariates_data_complete[['cell']], ][['participant']]))
        # store result
        res_per_comparison[[paste(region, gene, variant)]] <- interaction_model_df
      }
    }
  }
  # merge all results
  res_all <- do.call('rbind', res_per_comparison)
  return(res_all)
}


write_empty_result <- function(output_loc) {
  # gz file ends with .gz
  if (grepl('.gz$', output_loc)) {
    # gzip if ends with .gz
    con <- gzfile(output_loc, 'w')
    close(con)
  }
  else {
    file.create(output_loc)
  }
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
  make_option(c("-v", "--covariates_file"), type="character", default=NULL, 
             help="accessibility filename for chunk", metavar="character"), 
  make_option(c("-f", "--fixed_effects"), type="character", default=NULL,
              help="comman separated list of fixed effects to correct for [default= %default]", metavar="character"),
  make_option(c("-r", "--random_effects"), type="character", default=NULL,
              help="comma separated list of random effects to correct for [default= %default]", metavar="character"),
  make_option(c("-n", "--accessibility_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on accessibility before modelling [default: %default]"), 
  make_option(c("-y", "--expression_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on expression before modelling [default: %default]"), 
  make_option(c("-t", "--interaction_terms"), type="character", default=NULL,
              help="interaction to model", metavar='character'), 
  make_option(c("-g", "--genotype_loc"), type="character", default=NULL,
              help="genotype file in plink1 format, without extension", metavar="character"), 
  make_option(c("-b", "--barcode_column"), type="character", default='barcode_lane', 
              help="barcode column for metadata", metavar="character")
)


# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
# input directory
in_dir <- NULL
# location of the output
output_loc <- NULL
# location of the region-to-gene files
confinement_loc <- NULL
# location of SMF
smf_loc <- NULL
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
# barcode column in the covariates data
barcode_column <- NULL

if (debug) {
  # set all of the variables hardcoded for a testing debug run
  confinement_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement.tsv.gz'
  in_dir <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD4T/chr12-8899578-9674043/'
  smf_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/CD4T/smf.tsv.gz'
  output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/CD4T/chr12-8899578-9674043/'
  expression_file <- 'expression.tsv.gz'
  accessibility_file <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_CD4T_nonsparse_transposed.tsv.gz'
  genotype_loc <- '/groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr12'
  covariates_file <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
  fixed_effects_string <- 'region,genotype'
  random_effects_string <- 'sample_final,lane'
  interaction_terms_string <- 'genotype,region'
  barcode_column <- 'barcode_lane'
  
  
} else {
  # obligatory parameters without a default
  if (is.null(opt[['in']])) {
    error("i/--in is an obligatory parameter")
  } else {
    in_dir <- opt[['in']]
  }
  if (is.null(opt[['out']])) {
    error("o/--out is an obligatory parameter")
  } else {
    output_loc <- opt[['out']]
  }
  if (is.null(opt[['confinement']])) {
    error("c/--confinement is an obligatory parameter")
  } else {
    confinement_loc <- opt[['confinement']]
  }
  if (is.null(opt[['smf_loc']])) {
    error("s/--smf_loc is an obligatory parameter")
  } else {
    smf_loc <- opt[['smf_loc']]
  }
  # if (is.null(opt[['covariates_file']])) {
  #   error("v/--covariates_file is an obligatory parameter")
  # } else {
  #   covariates_file <- opt[['covariates_file']]
  # }
  if (is.null(opt[['genotype_loc']])) {
    error("g/--genotype_loc is an obligatory parameter")
  } else {
    genotype_loc <- opt[['genotype_loc']]
  }
  if (is.null(opt[['fixed_effects']])) {
    error("f/--fixed_effects is an obligatory parameter")
  } else {
    fixed_effects_string <- opt[['fixed_effects']]
  }
  if (is.null(opt[['interaction_terms']])) {
    error("t/--interaction_terms is an obligatory parameter")
  } else {
    interaction_terms_string <- opt[['interaction_terms']]
  }
  # parameters that have a sane default
  random_effects_string <- opt[['random_effects']]
  barcode_column <- opt[['barcode_column']]
  expression_file <- opt[['expression_file']]
  accessibility_file <- opt[['accessibility_file']]
  covariates_file <- opt[['covariates_file']]
  accessibility_gausnorm <- opt[['accessibility_gausnorm']]
  expression_gausnorm <- opt[['expression_gausnorm']]
}

# make the full path to the expression data
full_exp_path <- NULL
# depending on if it is an absolute path, we do things differently
if (startsWith(expression_file, '/')) {
  full_exp_path <- expression_file
} else {
  full_exp_path <- paste(in_dir, expression_file, sep = '/')
}
# same for the accessibility/TF data
full_acc_path <- NULL
if (startsWith(accessibility_file, '/')) {
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
confinement <- fread(confinement_loc, header = T, sep = '\t', )
# set harmonized column names to make it easier for ourselves
colnames(confinement) <- c('variant', 'region', 'gene')

# initialize variables
expression_data <- NULL
accessibility_data <- NULL
# check if there is expression data
if (length(count.fields(full_exp_path)) > 1) {
  # read the expression file
  # expression_data <- read.table(full_exp_path, header = T, sep = '\t', check.names = F, row.names = 1)
  expression_data <- fread(full_exp_path, header = T, sep = '\t', check.names = F, skip = 1)
  # read the header
  expression_data_header_line <- readLines(full_exp_path, n = 1)
  # split by sep
  expression_data_header <- strsplit(expression_data_header_line, '\t')[[1]]
  # add this header
  if (length(expression_data_header) == ncol(expression_data)) {
    colnames(expression_data) <- expression_data_header
  } else {
    # otherwise we need an extra column
    colnames(expression_data) <- c('gene', expression_data_header)
  }
  # also set the first column name so we can refer to it later
  colnames(expression_data)[[1]] <- 'gene'
} else {
  warning('no data fields for expression data, will do no further work')
  # to avoid further nesting, we'll make a dummy entry that makes it so that we dont continue further
  expression_data <- data.table('gene' = c())
}
# check if there is TF/accessibility data
if (length(count.fields(full_exp_path)) > 1) {
  # read the TF/accessibility data
  # accessibility_data <- read.table(full_acc_path, header = T, sep = '\t', check.names = F, row.names = 1)
  accessibility_data <- fread(full_acc_path, header = T, sep = '\t', check.names = F, skip = 1)
  # read the header
  accessibility_data_header_line <- readLines(full_acc_path, n = 1)
  # split by sep
  accessibility_data_header <- strsplit(accessibility_data_header_line, '\t')[[1]]
  # add this header
  if (length(accessibility_data_header) == ncol(accessibility_data)) {
    colnames(accessibility_data) <- accessibility_data_header
  } else {
    # otherwise we need an extra column
    colnames(accessibility_data) <- c('region', accessibility_data_header)
  }
  # set same colnames always
  colnames(accessibility_data)[[1]] <- 'region'
} else {
  warning('no data fields for accessibility data, will do no further work')
  # to avoid further nesting, we'll make a dummy entry that makes it so that we dont continue further
  accessibility_data <- data.table('region' = c())
}
  


# expression_data <- cbind(data.frame('gene' = rownames(expression_data)), expression_data)
# accessibility_data <- cbind(data.frame('region' = rownames(accessibility_data)), accessibility_data)
# subset both sets
expression_data <- expression_data[!is.na(expression_data[['gene']]) & expression_data[['gene']] %in% confinement[['gene']], ]
accessibility_data <- accessibility_data[!is.na(accessibility_data[['region']]) & accessibility_data[['region']] %in% confinement[['region']], ]

# format output loc
tsv_output_loc_full <- paste(output_loc, 'result.tsv.gz', sep = '/')
# set output loc as the tsv
output_loc_full <- tsv_output_loc_full
# gz file ends with .gz
if (grepl('.gz$', tsv_output_loc_full)) {
  # gzip if ends with .gz
  output_loc_full <- gzfile(tsv_output_loc_full)
}
# initialize the result
interaction_result <- NULL

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
      colnames(smf) <- c('participant', 'cell')
      # intersect the smf with the expression data and accessibility/TF data
      intersecting_cells <- intersect(smf[['cell']], colnames(expression_data))
      intersecting_cells <- intersect(intersecting_cells, colnames(accessibility_data))
      
      # check if there is any data left
      if (length(intersecting_cells) > 0) {
        # subset the matrices
        smf <- smf[smf[['cell']] %in% intersecting_cells, ]
        # get the columns
        accessibility_columns <- c('region', intersecting_cells)
        expression_columns <- c('gene', intersecting_cells)
        accessibility_data <- accessibility_data[, ..accessibility_columns]
        expression_data <- expression_data[, ..expression_columns]
        
        # split the fixed effects
        fixed_effects <- c()
        if (!is.null(fixed_effects_string) & !is.na(fixed_effects_string) & fixed_effects_string != '') {
          fixed_effects <- strsplit(fixed_effects_string, ',')[[1]]
        }
        # warn if we are not including the genotype
        if (!('genotype' %in% fixed_effects)) {
          warning(paste('\'genotype\' term not present in fixed effects!\n'))
        }
        # split random effects
        random_effects <- c()
        if (!is.null(random_effects_string) & !is.na(random_effects_string) & random_effects_string != '') {
          random_effects <- strsplit(random_effects_string, ',')[[1]]
        }
        # split random effects
        interactions <- c()
        if (!is.null(interaction_terms_string) & !is.na(interaction_terms_string) & interaction_terms_string != '') {
          interactions <- strsplit(interaction_terms_string, ',')[[1]]
        }
        # if a covariate matrix was supplied, we'll load it
        if (!is.null(full_covariates_path)) {
          covariates_data <- fread(full_covariates_path, header = T, sep = '\t')
          # the first column should be the cell
          covariates_data <- cbind(data.frame('cell' = covariates_data[[barcode_column]]), covariates_data)
          # intersect this
          intersecting_cells <- intersect(intersecting_cells, covariates_data[['cell']])
          # subset the covariates
          covariates_data <- covariates_data[covariates_data[['cell']] %in% intersecting_cells, ]
          # get the expression and accessibility/TF columns again
          accessibility_columns <- c('region', intersecting_cells)
          expression_columns <- c('gene', intersecting_cells)
          # and subset
          accessibility_data <- accessibility_data[, ..accessibility_columns]
          expression_data <- expression_data[, ..expression_columns]
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
          # order cells
          intersecting_cells <- intersecting_cells[order(intersecting_cells)]
          # make columns again
          accessibility_columns <- c('region', intersecting_cells)
          expression_columns <- c('gene', intersecting_cells)
          # and use this order
          accessibility_data <- accessibility_data[, ..accessibility_columns]
          expression_data <- expression_data[, ..expression_columns]
          covariates_data <- covariates_data[match(intersecting_cells, covariates_data[['cell']]), ]
          smf <- smf[match(intersecting_cells, smf[['cell']])]
          # perform the analysis
          message('Starting analysis..\n')
          # into a variable
          interaction_result <- do_interaction_analysis(
            expression_data = expression_data, 
            accessibility_data = accessibility_data, 
            genotype_data = genotypes, 
            smf = smf, 
            confinement = confinement,
            covariates_data = covariates_data, 
            fixed_effects = fixed_effects, 
            random_effects = random_effects, 
            interactions = interactions
          )
          # write result
          write.table(interaction_result, output_loc_full, sep = '\t', row.names = F, col.names = T, quote = F)
          # make a checksum
          mdfiver::create_sha256_for_file(tsv_output_loc_full)
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

# if the interaction result is still null, we didn't end up doing anything
if (is.null(interaction_result)) {
  # so we'll store an empty file
  write_empty_result(tsv_output_loc_full)
}
