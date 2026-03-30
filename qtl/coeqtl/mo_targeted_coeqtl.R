#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_targeted_coeqtl.R
# Function: perform targeted co-eQTL analysis
# 
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)
# use plink files
library(snpStats)
# transformation into gaussian normal distribution
library(bestNormalize)
# for the model
library(lme4)
library(lmerTest)
# for seurat object
library(Seurat)
# for plots
library(ggplot2)
library(cowplot)


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


gausnorm_independent_variable <- function(x, boxcox=F, min_value=1e-6) {
  # initialize value 
  y <- NULL
  # only if numeric we can convert
  if (is.numeric(x)) {
    if (!is.null(min_value) && boxcox) {
      x[x < min_value] <- min_value
    }
    if (boxcox) {
      y <- boxcox(x)$x.t
    } else {
      y <- yeojohnson(x)$x.t
    }
  }
  else {
    y <- x
  }
  return(y)
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
  row_created <- data.frame(matrix(, nrow = 1, ncol = length(covariates) * length(stats)))
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
                                    genotype_data, 
                                    smf, 
                                    confinement,
                                    covariates_data=NULL, 
                                    fixed_effects=c('lane','region','genotype'), 
                                    random_effects=c('sample_final'), 
                                    interactions=c('genotype','expression2'), 
                                    family = 'gaussian', 
                                    expression_gausnorm=T, 
                                    expression_boxcox=F) {
  # create formula
  base_formula <- get_formula(var_of_interest = 'expression', fixed_effects = fixed_effects, random_effects = random_effects)
  base_formula_string <- (Reduce(paste, deparse(base_formula)))
  message(paste('Using base formula:', base_formula_string, ''))
  # paste together
  interactions_formula_part <- paste(interactions, collapse = '*')
  # create interaction formula
  interaction_formula_string <- paste(base_formula_string, interactions_formula_part, sep = ' + ')
  message(paste('Using interaction formula:', interaction_formula_string, ''))
  interaction_formula <- as.formula(interaction_formula_string)
  # conver to dataframes
  expression_data <- data.frame(expression_data)
  smf <- data.frame(smf)
  confinement <- data.frame(confinement)
  if (!is.null(covariates_data)) {
    covariates_data <- data.frame(covariates_data)
    # get which things we need from that dataframe
    covariate_columns <- setdiff(c(fixed_effects, random_effects, interactions), c('expression', 'expression2', 'genotype'))
    # and subset the covariates data to that
    covariates_data <- covariates_data[, c('cell', covariate_columns)]
  }
  # put all results in a list
  res_per_comparison <- list()
  # check each region
  for (genea in unique(confinement[['genea']])) {
    # extract the genes and variants
    confinement_genea <- confinement[confinement[['genea']] == genea, ]
    # the specific genes then
    genesb_genea <- unique(confinement_genea[['geneb']])
    # subset to these genes
    expression_data_genea_and_genesb <- expression_data[
      rownames(expression_data) %in% c(genea, genesb_genea), 
    ]
    # check each gene
    for (geneb in unique(genesb_genea)) {
      # extract the gene
      genea_values <- as.vector(unlist(expression_data_genea_and_genesb[genea, ]))
      geneb_values <- as.vector(unlist(expression_data_genea_and_genesb[geneb, ]))
      # merge the metadata with the gene and the region
      covariates_data[['expression']] <- genea_values
      covariates_data[['expression2']] <- geneb_values
      # gausnorm them if requested
      if(expression_gausnorm) {
        covariates_data[['expression']] <- gausnorm_independent_variable(covariates_data[['expression']], expression_boxcox)
        covariates_data[['expression2']] <- gausnorm_independent_variable(covariates_data[['expression2']], expression_boxcox)
      }
      # get the variants for this region-gene combination
      variants_genea_geneb <- unique(confinement_genea[confinement_genea[['genea']] == genea &
                                                         confinement_genea[['geneb']] == geneb, ][['variant']])
      # check each variant
      for (variant in variants_genea_geneb) {
        # extract genotypes
        genotype <- genotype_data$genotypes[smf[['participant']], variant]
        # then to numeric
        genotype_numeric <- as.vector(as(genotype, 'numeric'))
        # add the genotype
        covariates_data[['genotype']] <- genotype_numeric
        # keep only complete cases
        covariates_data_complete <- covariates_data[complete.cases(covariates_data), ]
        # and only finite values
        covariates_data_complete <- covariates_data_complete[is.finite(covariates_data_complete[['expression']]) & is.finite(covariates_data_complete[['expression2']]) & is.finite(covariates_data_complete[['genotype']]), ]
        # check if there is any data left
        if (nrow(covariates_data_complete) > 0) {
          # initialize variables
          base_model <- NULL
          interaction_model <- NULL
          ftest_res <- NULL
          anova_test_used <- NULL
          # try to do 
          tryCatch({
            # depending on the family, the calls and anovas are different
            if (family == 'poisson') {
              # model without interaction
              base_model <- lme4::glmer(formula = base_formula, data = covariates_data_complete, family = poisson)
              # model with interaction
              interaction_model <- lme4::glmer(formula = interaction_formula, data = covariates_data_complete, family = poisson)
              # try to do the Chi-squared test
              tryCatch({
                # check if they are different
                ftest_res <- anova(base_model, interaction_model, refit = FALSE, test = 'Chisq')
                # catch error if model did not converge or some other issue
              }, error = function(e) {
                warning(paste('Error in chi-square', genea, 'gene', geneb, 'variant', variant, ':', e$message, '. This can happen if the model fails to converge'))
                ftest_res <- data.frame('y' = c(NA, NA))
                rownames(ftest_res) <- c('base_model', 'interaction_model')
              })
              anova_test_used <- 'LRT'
            }
            else if (family == 'gaussian') {
              # base model
              base_model <- lmerTest::lmer(formula = base_formula, data = covariates_data_complete)
              # interaction model
              interaction_model <- lmerTest::lmer(formula = interaction_formula, data = covariates_data_complete)
              # try to do F test
              tryCatch({
                # check if they are different
                ftest_res <- anova(base_model, interaction_model, refit = FALSE, test = 'F')
                # catch error if model did not converge or some other issue
              }, error = function(e) {
                warning(paste('Error in F-test', genea, 'gene', geneb, 'variant', variant, ':', e$message, '. This can happen if the model fails to converge'))
                ftest_res <- data.frame('y' = c(NA, NA))
                rownames(ftest_res) <- c('base_model', 'interaction_model')
              })
              anova_test_used <- 'F'
            }
            else {
              stop(paste0('unknown family ', family, ', only valid families are gaussian and poisson'))
            }
            
          }, error = function(e) {
            warning(paste('Error in model fitting', genea, geneb, variant, ':', e$message, '. This can happen if the model fails to converge'))
          })
          # initialize the dataframe to add model to
          interaction_model_df_base <- data.frame('variant' = c(variant), 'genea' = c(genea), 'geneb' = c(geneb))
          # initialize the model df
          interaction_model_df <- NULL
          # if we have a model, we can convert to a df
          if(!is.null(interaction_model)) {
            # convert the model to a df
            interaction_model_df <- model_to_row(interaction_model)
            # and what we actually tested
            interaction_model_df <- cbind(interaction_model_df_base, interaction_model_df)
          } else if(!is.null(base_model)) {
            # if we at least have a base model, we can still convert that to a df
            interaction_model_df <- model_to_row(base_model)
            # and what we actually tested
            interaction_model_df <- cbind(interaction_model_df_base, interaction_model_df)
          } else {
            # if the model failed to fit, we still want to have a row for this combination
            interaction_model_df <- interaction_model_df_base
          }
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
          res_per_comparison[[paste(genea, geneb, variant)]] <- data.table(interaction_model_df)
        } else {
          warning(paste('No data left for gene', genea, 'gene', geneb, 'variant', variant, 'after complete cases. Skipping this combination.'))
        }
      }
    }
  }
  # merge all results
  res_all <- rbindlist(res_per_comparison, fill = T)
  return(res_all)
}


get_interaction_inputs <- function(expression_data, 
                                   genotype_data, 
                                   smf, 
                                   confinement,
                                   covariates_data=NULL, 
                                   fixed_effects=c('lane','region','genotype'), 
                                   random_effects=c('sample_final'), 
                                   interactions=c('genotype','region'), 
                                   expression_gausnorm=T, 
                                   expression_boxcox=F
) {

  smf <- data.frame(smf)
  confinement <- data.frame(confinement)
  if (!is.null(covariates_data)) {
    covariates_data <- data.frame(covariates_data)
    # get which things we need from that dataframe
    covariate_columns <- setdiff(c(fixed_effects, random_effects, interactions), c('expression', 'expression2', 'genotype'))
    # and subset the covariates data to that
    covariates_data <- covariates_data[, c('cell', covariate_columns)]
  }
  # put all results in a list
  res_per_comparison <- list()
  # check each region
  for (genea in unique(confinement[['genea']])) {
    # extract the genes and variants
    confinement_genea <- confinement[confinement[['genea']] == genea, ]
    # the specific genes then
    genesb_genea <- unique(confinement_genea[['geneb']])
    # subset to these genes
    expression_data_genea_and_genesb <- expression_data[
      rownames(expression_data) %in% c(genea, genesb_genea), 
    ]
    # check each gene
    for (geneb in unique(genesb_genea)) {
      # extract the gene
      genea_values <- as.vector(unlist(expression_data_genea_and_genesb[genea, ]))
      geneb_values <- as.vector(unlist(expression_data_genea_and_genesb[geneb, ]))
      # merge the metadata with the gene and the region
      covariates_data[['expression']] <- genea_values
      covariates_data[['expression2']] <- geneb_values
      # gausnorm them if requested
      if(expression_gausnorm) {
        covariates_data[['expression_raw']] <- covariates_data[['expression']]
        covariates_data[['expression2_raw']] <- covariates_data[['expression2']]
        covariates_data[['expression']] <- gausnorm_independent_variable(covariates_data[['expression']], expression_boxcox)
        covariates_data[['expression2']] <- gausnorm_independent_variable(covariates_data[['expression2']], expression_boxcox)
      }
      # get the variants for this region-gene combination
      variants_genea_geneb <- unique(confinement_genea[confinement_genea[['genea']] == genea &
                                                         confinement_genea[['geneb']] == geneb, ][['variant']])
      # check each variant
      for (variant in variants_genea_geneb) {
        # extract genotypes
        genotype <- genotype_data$genotypes[smf[['participant']], variant]
        # then to numeric
        genotype_numeric <- as.vector(as(genotype, 'numeric'))
        # add the genotype
        covariates_data[['genotype']] <- genotype_numeric
        # keep only complete cases
        covariates_data_complete <- covariates_data[complete.cases(covariates_data), ]
        # store result
        res_per_comparison[[paste(genea, geneb, variant)]] <- covariates_data_complete
      }
    }
  }
  return(res_per_comparison)
}


calculate_per_sample_correlation <- function(full_variates_table, formula_string='expression~expression2', correlation=T, method='spearman', sample_column='sample_id', beta_variate_column='expression2', family='gaussian') {
  # make into datatable
  full_variates_table <- data.table(full_variates_table)
  # we'll store in a list first
  cor_per_sample <- list()
  # list the samples
  samples_present <- unique(full_variates_table[[sample_column]])
  # filter where we don't know the sample
  samples_present <- samples_present[!is.na(samples_present)]
  # check each sample
  for (sample_present in samples_present) {
    # subset the table
    variates_table_sample <- full_variates_table[
      !is.na(full_variates_table[[sample_column]]) & full_variates_table[[sample_column]] == sample_present, 
    ]
    # get number of cells
    ncell <- nrow(variates_table_sample)
    # now calculate a correlation
    estimate <- NULL
    p <- NULL
    # try to do this, we might error if we have too few observations
    tryCatch({
      if (correlation) {
        # reformat the formula
        formula_string <- gsub(' ', '', formula_string)
        # then split by predictor
        formula_string_split <- strsplit(formula_string, '~')[[1]]
        # check if of correct length
        if (length(formula_string_split) > 2) {
          stop('split contains more than a two values. A correlation can only be made up of two variables')
        }
        else if (length(formula_string_split) < 2) {
          stop('split contains less than two values. A correlation can only be made up of two variables')
        } else {
          # do the correlation test
          cor_test <- cor.test(y = variates_table_sample[[formula_string_split[[1]]]], 
                               x = variates_table_sample[[formula_string_split[[2]]]], 
                               method = method)
          # extract p
          p <- cor_test$p.value
          estimate <- as.vector(cor_test$estimate[1])
        }
      } else {
        # initialize regression model
        regression_model <- NULL
        # depending on the family, the calls are different
        if (family == 'poisson') {
          # use poisson model
          regression_model <- lme4::glmer(formula = as.formula(formula_string), data = variates_table_sample, family = poisson)
        }
        else if (family == 'gaussian') {
          # use gaussian model
          regression_model <- lmerTest::lmer(formula = as.formula(formula_string), data = variates_table_sample)
        } else {
          # error if weird formula is given
          stop(paste0('unknown family ', family, ', only valid families are gaussian and poisson when using a regression model'))
        }
        # convert model
        regression_model_row <- model_to_row(regression_model)
        # extract the values
        p <- as.vector(unlist(regression_model_row[paste(beta_variate_column, 'p', sep = '_')][1]))
        estimate <- as.vector(unlist(regression_model_row[paste(beta_variate_column, 'beta', sep = '_')][1]))
      }
    }, error = function(e) {
      warning(paste('Error in correlation or regression model for sample', ':', e$message, '. This can happen if the model fails to converge'))
    })
    # make into df
    cor_per_sample[[sample_present]] <- data.table('sample' = c(sample_present), 'estimate' = c(estimate), 'p' = c(p), 'ncell' = c(ncell))
  }
  # merge all
  cor_all <- rbindlist(cor_per_sample, fill = T)
  return(cor_all)
}


calculate_per_sample_prediction <- function(full_variates_table, formula_string='expression~expression2', correlation=T, method='spearman', sample_column='sample_id', beta_variate_column='region', family='gaussian', base_lm=F) {
  # make into datatable
  full_variates_table <- data.table(full_variates_table)
  # we'll store in a list first
  cor_per_sample <- list()
  # list the samples
  samples_present <- unique(full_variates_table[[sample_column]])
  # filter where we don't know the sample
  samples_present <- samples_present[!is.na(samples_present)]
  # check each sample
  for (sample_present in samples_present) {
    # subset the table
    variates_table_sample <- full_variates_table[
      !is.na(full_variates_table[[sample_column]]) & full_variates_table[[sample_column]] == sample_present, 
    ]
    # get number of cells
    ncell <- nrow(variates_table_sample)
    # now calculate a correlation
    estimate <- NULL
    p <- NULL
    # try to do this, we might error if we have too few observations
    tryCatch({
      if (correlation) {
        # reformat the formula
        formula_string <- gsub(' ', '', formula_string)
        # then split by predictor
        formula_string_split <- strsplit(formula_string, '~')[[1]]
        # check if of correct length
        if (length(formula_string_split) > 2) {
          stop('split contains more than a two values. A correlation can only be made up of two variables')
        }
        else if (length(formula_string_split) < 2) {
          stop('split contains less than two values. A correlation can only be made up of two variables')
        } else {
          # do the correlation test
          cor_test <- cor.test(y = variates_table_sample[[formula_string_split[[1]]]], 
                               x = variates_table_sample[[formula_string_split[[2]]]], 
                               method = method)
          # extract p
          p <- cor_test$p.value
          estimate <- as.vector(cor_test$estimate[1])
        }
      } else {
        # initialize regression model
        regression_model <- NULL
        # depending on the family, the calls are different
        if (family == 'poisson') {
          # use poisson model
          regression_model <- lme4::glmer(formula = as.formula(formula_string), data = variates_table_sample, family = poisson)
        }
        else if (family == 'gaussian') {
          # use gaussian model
          if (base_lm) {
            regression_model <- lm(formula = as.formula(formula_string), data = variates_table_sample)
          } else {
            regression_model <- lmerTest::lmer(formula = as.formula(formula_string), data = variates_table_sample)
          }
        } else {
          # error if weird formula is given
          stop(paste0('unknown family ', family, ', only valid families are gaussian and poisson when using a regression model'))
        }
        variates_table_sample[['predicted']] <- as.vector(unlist(predict(regression_model, variates_table_sample)))
      }
    }, error = function(e) {
      warning(paste('Error in correlation or regression model for sample', region, 'gene', gene, 'variant', variant, ':', e$message, '. This can happen if the model fails to converge'))
    })
    # make into df
    cor_per_sample[[sample_present]] <- variates_table_sample
  }
  # merge all
  cor_all <- rbindlist(cor_per_sample, fill = T)
  return(cor_all)
}


####################
# Main code        #
####################

# location of the genotype
genotype_loc <- '/groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr8'
# genes to correlate with one another
var_gene_gene_confinement <- data.frame(
  'variant' = c('8:11491452:G:A', '8:11491452:G:A'), 
  'genea' = c('FAM167A', 'BLK'), 
  'geneb' = c('BLK', 'FAM167A')
)
# locatin of the Seurat object
seurat_object_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds'
# location of SMF file
smf_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/all/smf.tsv.gz'
# effects
fixed_effects_string <- 'expression2,genotype,nCount_RNA,condition_final,celltype_imputed_lowerres'
random_effects_string <- 'sample_final,lane'
interaction_terms_string <- 'genotype,expression2'
aggregate_columns_string <- 'sample_final,condition_final'
# whether to gausnorm
expression_gausnorm <- T
expression_boxcox <- F
correlations_gausnorm <- T
ncell_cutoff <- 10

# split the fixed effects
fixed_effects <- c()
if (!is.null(fixed_effects_string) && !is.na(fixed_effects_string) && fixed_effects_string != '') {
  fixed_effects <- strsplit(fixed_effects_string, ',')[[1]]
}
# warn if we are not including the genotype
if (!('genotype' %in% fixed_effects)) {
  warning(paste('\'genotype\' term not present in fixed effects!\n'))
}
# split random effects
random_effects <- c()
if (!is.null(random_effects_string) && !is.na(random_effects_string) && random_effects_string != '') {
  random_effects <- strsplit(random_effects_string, ',')[[1]]
}
# split random effects
interactions <- c()
if (!is.null(interaction_terms_string) && !is.na(interaction_terms_string) && interaction_terms_string != '') {
  interactions <- strsplit(interaction_terms_string, ',')[[1]]
}
# split aggregate columns
aggregate_columns <- c()
if (!is.null(aggregate_columns_string) && !is.na(aggregate_columns_string) && aggregate_columns_string != '') {
  aggregate_columns <- strsplit(aggregate_columns_string, ',')[[1]]
}


# load seurat object
seurat_object <- readRDS(seurat_object_loc)
# SUBSET!!!
# seurat_object <- seurat_object[, 
#                                !is.na(seurat_object@meta.data[['celltype_imputed_lowerres']]) & 
#                                  seurat_object@meta.data[['celltype_imputed_lowerres']] == 'B']
# SUBSET!!!

# extract covariates
covariates_data <- seurat_object@meta.data
# set the cell id as the first column
covariates_data <- cbind(data.table('cell' = rownames(covariates_data)), data.table(covariates_data))
# read smf
smf <- fread(smf_loc, header = T, sep = '\t')
# harmonize names
colnames(smf) <- c('participant', 'cell')
# intersect the smf with the expression data and accessibility/TF data
intersecting_cells <- intersect(smf[['cell']], covariates_data[['cell']])
# subset all
seurat_object <- seurat_object[, intersecting_cells]
smf <- smf[match(intersecting_cells, smf[['cell']]), ]
covariates_data <- covariates_data[match(intersecting_cells, covariates_data[['cell']]), ]
# extract the actual expression
expression_data <- seurat_object@assays$MJ@layers$data
# set the feature names
rownames(expression_data) <- rownames(data.frame(seurat_object@assays$MJ@features))
# and cell names
colnames(expression_data) <- colnames(seurat_object)
# keep only genes we care about
genes_both <- intersect(rownames(expression_data), unique(c(var_gene_gene_confinement[['genea']], var_gene_gene_confinement[['geneb']])))
# then subset the data
expression_data <- expression_data[genes_both, ]
var_gene_gene_confinement <- var_gene_gene_confinement[
  var_gene_gene_confinement[['genea']] %in% genes_both &
  var_gene_gene_confinement[['geneb']] %in% genes_both, 
  
]
# get the variants from the confinement file
variants <- var_gene_gene_confinement[['variant']]
# read the bim
variants_in_gt <- fread(paste(genotype_loc, '.bim', sep = ''), header = F)[[2]]
# get overlapping variants
overlapping_variants <- intersect(variants, variants_in_gt)
# initialize genotype data
genotypes <- NULL
if (length(overlapping_variants) > 0) {
  # read the genotypes, but only those in the file and in the confinement
  genotypes <- read.plink(
    bed = paste(genotype_loc, '.bed', sep = ''),
    bim = paste(genotype_loc, '.bim', sep = ''),
    fam = paste(genotype_loc, '.fam', sep = ''), 
    select.snps = overlapping_variants
  )
  # filter the confinement on the variants we have in the genotype data as well
  var_gene_gene_confinement <- var_gene_gene_confinement[var_gene_gene_confinement[['variant']] %in% overlapping_variants, ]
}
# perform the analysis
message('Starting analysis..')
# into a variable
interaction_result <- do_interaction_analysis(
  expression_data = expression_data, 
  genotype_data = genotypes, 
  smf = smf, 
  confinement = var_gene_gene_confinement,
  covariates_data = covariates_data, 
  fixed_effects = fixed_effects, 
  random_effects = random_effects, 
  interactions = interactions, 
  expression_gausnorm = expression_gausnorm, 
  expression_boxcox = expression_boxcox
)

# do the co-eQTL version
interaction_inputs <- get_interaction_inputs(
  expression_data = expression_data, 
  genotype_data = genotypes, 
  smf = smf, 
  confinement = var_gene_gene_confinement,
  covariates_data = covariates_data, 
  fixed_effects = fixed_effects, 
  random_effects = random_effects, 
  interactions = interactions, 
  expression_gausnorm = expression_gausnorm, 
  expression_boxcox = expression_boxcox
)

# store the plots
cor_gt_plots <- list()
# the inputs
cor_gt_dfs <- list()
# and the results
cor_gt_results <- list()
# get co-eQTL style plots
for (confinement_i in 1:nrow(var_gene_gene_confinement)) {
  # extract the variant, region and gene
  variant <- as.vector(unlist(var_gene_gene_confinement[confinement_i, 'variant']))
  genea <- as.vector(unlist(var_gene_gene_confinement[confinement_i, 'genea']))
  geneb <- as.vector(unlist(var_gene_gene_confinement[confinement_i, 'geneb']))
  # paste together the naming
  confinement_name <- paste(genea, geneb, variant)
  # check if this combination is in the results
  if (confinement_name %in% names(interaction_inputs)) {
    # extract the plottable table
    plot_df <- interaction_inputs[[confinement_name]]
    # add an aggregated column for the sample, by pasting the aggregate columns together
    plot_df[['aggregated_sample']] <- apply(plot_df[, aggregate_columns, drop = F], 1, function(x) paste(x, collapse = '_'))
    # add an aggregated column for the sample, by pasting the aggregate columns together
    covariates_data[['aggregated_sample']] <- apply(covariates_data[, ..aggregate_columns, drop = F], 1, function(x) paste(x, collapse = '_'))
    # get the per-sample plot
    per_sample_df <- calculate_per_sample_correlation(plot_df, sample_column = 'aggregated_sample', formula_string = paste(c('expression', 'expression2'), sep = '~', collapse = '~'))
    # store original estimate
    per_sample_df[['estimate_raw']] <- per_sample_df[['estimate']]
    # gausnorm if requested
    if (correlations_gausnorm) {
      per_sample_df[['estimate']] <- gausnorm_independent_variable(per_sample_df[['estimate_raw']])
    }
    # get the predictions per sample
    per_sample_predictions <- calculate_per_sample_prediction(plot_df, sample_column = 'aggregated_sample', correlation = F, formula_string = 'expression~expression2', base_lm=T)
    # take the unique sets of the covariates from the plot df, to add this to the per sample df
    unique_covariate_columns <- unique(c(fixed_effects, random_effects, interactions, aggregate_columns))
    # but remove region and expression and celltype and lane
    # unique_covariate_columns <- setdiff(unique_covariate_columns, c('expression2', 'expression', 'lane', 'nCount_RNA'))
    unique_covariate_columns <- setdiff(unique_covariate_columns, c('expression2', 'expression', 'nCount_RNA', 'celltype_imputed_lowerres'))
    # subset the plot df to these columns and the sample column, and take unique rows
    plot_df_unique_covariates <- unique(plot_df[, c('aggregated_sample', unique_covariate_columns)])
    # then add this to the plot df
    per_sample_df <- merge(per_sample_df, plot_df_unique_covariates, by.x = 'sample', by.y = 'aggregated_sample', all.x = T)
    # make character string
    per_sample_df[['gt']] <- as.character(per_sample_df[['genotype']])
    # create a formula
    base_formula <- get_formula(var_of_interest = 'estimate', fixed_effects = c(setdiff(fixed_effects, c('expression2', 'nCount_RNA')), 'ncell'), random_effects = setdiff(random_effects, c('sample_final', 'lane')))
    base_formula <- get_formula(var_of_interest = 'estimate', fixed_effects = c(setdiff(fixed_effects, c('expression2', 'nCount_RNA')), 'ncell'), random_effects = random_effects)
    # initialize variable
    lm_gt_to_cor_table_base <- data.table('variant' = c(variant), 'genea' = c(genea), 'geneb' = c(geneb))
    lm_gt_to_cor_table <- NULL
    # try to do modelling
    tryCatch({
      # fit model with the genotype and the correlation
      lm_gt_to_cor <- NULL
      # if (!is.null(setdiff(random_effects, c('sample_final', 'lane'))) && length(setdiff(random_effects, c('sample_final', 'lane'))) > 0) {
      if (!is.null(random_effects) && length(random_effects) > 0) {
        # use glm if we have random effects
        lm_gt_to_cor <- lmerTest::lmer(formula = base_formula, data = per_sample_df[per_sample_df$ncell >= ncell_cutoff, ])
      } else {
        # or simple model if here are not
        message('calling base lm due to no random effects specified')
        lm_gt_to_cor <- lm(data = per_sample_df[per_sample_df$ncell >= ncell_cutoff, ], formula = base_formula)
      }
      # extract p value for the genotype term
      lm_gt_to_cor_summary <- summary(lm_gt_to_cor)
      lm_gt_to_cor_p <- NULL
      # the p values is always the last
      lm_gt_to_cor_p <- lm_gt_to_cor_summary[['coefficients']]['genotype', ncol(lm_gt_to_cor_summary[['coefficients']])]
      # convert the result to a table
      lm_gt_to_cor_table <- model_to_row(lm_gt_to_cor)
      # add the variant, region and gene to the table
      lm_gt_to_cor_table <- cbind(lm_gt_to_cor_table_base, lm_gt_to_cor_table)
    }, error = function(e) {
      warning(paste('Error in model fitting', genea, geneb, variant, ':', e$message))
      # make empty table
      lm_gt_to_cor_table <- lm_gt_to_cor_table_base
    })
    # add the number of samples
    lm_gt_to_cor_table[['nsample']] <- nrow(per_sample_df[per_sample_df[['ncell']] >= ncell_cutoff, ])
    # and the distribution of cells
    lm_gt_to_cor_table[['ncell']] <- paste(as.character(min(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])),
                                           as.character(quantile(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])[['25%']]),
                                           as.character(quantile(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])[['50%']]),
                                           as.character(quantile(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])[['75%']]),
                                           as.character(max(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])),
                                           sep = ';'
    )
    # check how many genotypes we have
    unique_genotypes <- unique(per_sample_df[['gt']])
    # filter genotypes
    unique_genotypes <- unique_genotypes[!is.na(unique_genotypes)]
    # and sort first
    unique_genotypes <- unique_genotypes[order(unique_genotypes)]
    # get colors for the genotypes
    geno_colors <- roycols::get_color_list(unique_genotypes)
    # plot both of them
    p_cor <- ggplot(data = per_sample_df[per_sample_df[['ncell']] >= ncell_cutoff, ], mapping = aes(x = gt, y = estimate_raw, fill = gt)) + 
      geom_boxplot(outlier.shape = NA) + 
      # geom_point() +
      # and add jitter
      geom_jitter(size = 0.5, alpha = 0.5, data = per_sample_df[per_sample_df[['ncell']] >= ncell_cutoff, ], mapping = aes(x = gt, y = estimate_raw, colour = ncell)) +
      # geom_jitter(size = 0.5, alpha = 0.5) + 
      scale_fill_manual(values = geno_colors) + 
      # and colour of ncell
      scale_colour_gradient2(low='blue', mid = 'white', high='red', midpoint = max(per_sample_df[['ncell']]/2)) + 
      xlab(paste('genotype')) + 
      ylab(paste(genea, 'expression ~', geneb, 'expression2')) + 
      labs(fill = 'Genotype', colour = 'Ncell') + 
      ggtitle(paste('cor', variant, genea, geneb, 'p < ', as.character(round(lm_gt_to_cor_p, digits = 5)))) +
      theme(legend.title = element_text(size=14), 
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            axis.text.y = element_text(size=12),
            axis.text.x = element_text(size=12),
            strip.text.x = element_text(size=12)) + 
      theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
    # get the covariate columns
    per_sample_predictions_columns <- c('cell', setdiff(colnames(covariates_data), colnames(per_sample_predictions)))
    # next, add the same information to our predicted expression values
    per_sample_predictions <- merge(per_sample_predictions, covariates_data[, ..per_sample_predictions_columns], by.x = 'cell', by.y = 'cell', all.x = T)
    # make character string
    per_sample_predictions[['gt']] <- as.character(per_sample_predictions[['genotype']])
    # put each genotype in a list
    geno_p_list <- list()
    # get the x range
    x_min <- min(per_sample_predictions[['expression2']])
    x_max <- max(per_sample_predictions[['expression2']])
    # get the y range
    y_min <- min(per_sample_predictions[['predicted']])
    y_max <- max(per_sample_predictions[['predicted']])
    # check each genotype
    for (unique_genotype in unique_genotypes) {
      # extract taht genotype
      per_sample_predictions_geno <- per_sample_predictions[per_sample_predictions[['gt']] == unique_genotype, ]
      # we'll colour all these by the genotype, which is the same
      color_list <- list()
      for (s in unique(per_sample_predictions_geno[['aggregated_sample']])) {
        color_list[[s]] <- geno_colors[[unique_genotype]]
      }
      # plot each sample
      # p_gt_group <- ggplot(data = per_sample_predictions_geno, mapping = aes(x = region, y = predicted, colour = aggregated_sample)) + 
      p_gt_group <- ggplot(data = per_sample_predictions_geno, mapping = aes(x = expression2, y = predicted, group = aggregated_sample)) + 
        geom_line(colour = geno_colors[[unique_genotype]]) + 
        xlim(c(x_min * 1.1, x_max * 1.1)) + 
        ylim(c(y_min * 1.1, y_max * 1.1)) +
        # and colour of ncell
        # scale_colour_manual(values = color_list) + 
        xlab(paste(geneb, 'expression')) + 
        ylab(paste(genea, 'expression')) + 
        labs(fill = 'sample', colour = 'sample') + 
        ggtitle(paste(unique_genotype)) +
        theme(legend.title = element_text(size=14), 
              legend.text = element_text(size=12),
              axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=12),
              strip.text.x = element_text(size=12)) + 
        theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
        theme(legend.position = 'none')
      # put in the list
      geno_p_list[[unique_genotype]] <- p_gt_group
    }
    # merge the ps
    geno_p_all <- plot_grid(plotlist = geno_p_list, nrow = 1)
    # place in table
    cor_gt_plots[[paste(variant, genea, geneb, 'cor', sep = '_')]] <- plot_grid(geno_p_all, p_cor, nrow = 2, rel_heights = c(1, 2)) 
    cor_gt_dfs[[paste(variant, genea, geneb, 'cor', sep = '_')]] <- per_sample_df
    cor_gt_results[[paste(variant, genea, geneb, 'cor', sep = '_')]] <- lm_gt_to_cor_table
  }
}
