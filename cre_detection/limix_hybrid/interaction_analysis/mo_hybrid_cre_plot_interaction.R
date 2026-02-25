#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_hybrid_cre_plot_interaction.R
# Function: plot interaction-eQTL at single-cell level with TF or ATAC as interaction terms
# Example: 
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
# to plot
library(ggplot2)
library(roycols)


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


get_interaction_inputs <- function(expression_data, 
                                    accessibility_data, 
                                    genotype_data, 
                                    smf, 
                                    confinement,
                                    covariates_data=NULL, 
                                    fixed_effects=c('lane','region','genotype'), 
                                    random_effects=c('sample_final'), 
                                    interactions=c('genotype','region')
                                   ) {
  # # create formula
  # base_formula <- get_formula(var_of_interest = 'expression', fixed_effects = fixed_effects, random_effects = random_effects)
  # base_formula_string <- (Reduce(paste, deparse(base_formula)))
  # message(paste('Using base formula:', base_formula_string, ''))
  # # paste together
  # interactions_formula_part <- paste(interactions, collapse = '*')
  # # create interaction formula
  # interaction_formula_string <- paste(base_formula_string, interactions_formula_part, sep = ' + ')
  # message(paste('Using interaction formula:', interaction_formula_string, ''))
  # interaction_formula <- as.formula(interaction_formula_string)
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
        # store result
        res_per_comparison[[paste(region, gene, variant)]] <- covariates_data_complete
      }
    }
  }
  return(res_per_comparison)
}


calculate_per_sample_correlation <- function(full_variates_table, formula_string='expression~region', correlation=T, sample_column='sample_id', beta_variate_column='region', family='gaussian') {
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
                               x = variates_table_sample[[formula_string_split[[2]]]])
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
        p <- regression_model_row[paste(beta_variate_column, 'p', sep = '_')][1]
        estimate <- regression_model_row[paste(beta_variate_column, 'beta', sep = '_')][1]
      }
    }, error = function(e) {
      warning(paste('Error in correlation or regression model for sample', sample_present, e$message))
      estimate <- NA
      p <- NA
    })
    # make into df
    cor_per_sample[[sample_present]] <- data.table('sample' = c(sample_present), 'estimate' = c(estimate), 'p' = c(p), 'ncell' = c(ncell))
  }
  # merge all
  cor_all <- rbindlist(cor_per_sample, fill = T)
  return(cor_all)
}


####################
# Settings         #
####################

# luck seed
set.seed(7777)


####################
# Main code        #
####################

# the chunk the data is in
chunk <- 'chr6-6385763-6874450'
# where the files are placed
in_dir_base <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/'
# which genotype to use
genotype_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr6'
# the variant to look at
variants <- c('6:6588959:G:GGA', '6:6569080:G:C', '6:6588959:G:GGA', '6:6569080:G:C', '6:6588959:G:GGA', '6:6569080:G:C')
# the TF or region
regions <- c('BCL11A_direct_+/+_(146g)', 'BCL11A_direct_+/+_(146g)', 'EBF1_direct_+/+_(142g)', 'EBF1_direct_+/+_(142g)', 'PAX5_direct_+/+_(115g)', 'PAX5_direct_+/+_(115g)')
# the gene
genes <- c('LY86', 'LY86', 'LY86', 'LY86', 'LY86', 'LY86')
# finally the cell type
cell_type <- 'monocyte'
# whether to gausnorm first
accessibility_gausnorm <- T
expression_gausnorm <- T
# model family
family <- 'gaussian'

# other parameters like in the actual mapping
in_dir <- paste(in_dir_base, cell_type, chunk, sep = '/')
smf_loc <- paste(in_dir_base, cell_type, 'smf.tsv.gz', sep = '/')
expression_file <- 'expression.tsv.gz'
accessibility_file <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_monocyte_nonsparse_transposed.tsv.gz'
covariates_file <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
fixed_effects_string <- 'region,genotype'
random_effects_string <- 'sample_final,lane'
interaction_terms_string <- 'genotype,region'
barcode_column <- 'barcode_lane'

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

# create confinement with just those parameters
confinement <- data.table('variant' = c(variants), 'region' = c(regions), 'gene' = c(genes))

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
if (length(count.fields(full_acc_path)) > 1) {
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

# subset both sets
expression_data_confined <- expression_data[!is.na(expression_data[['gene']]) & expression_data[['gene']] %in% confinement[['gene']], ]
accessibility_data_confined <- accessibility_data[!is.na(accessibility_data[['region']]) & accessibility_data[['region']] %in% confinement[['region']], ]

# initialize table with plot information
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
    confinement <- confinement[confinement[['gene']] %in% expression_data_confined[['gene']], ]
    # and regions or TFs we have
    confinement <- confinement[confinement[['region']] %in% accessibility_data_confined[['region']], ]
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
          message('Starting analysis..')
          # into a variable
          interaction_result <- get_interaction_inputs(
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

# create formula
base_formula <- get_formula(var_of_interest = 'expression', fixed_effects = fixed_effects, random_effects = random_effects)
base_formula_string <- (Reduce(paste, deparse(base_formula)))
# paste together
interactions_formula_part <- paste(interactions, collapse = '*')
# create interaction formula
interaction_formula_string <- paste(base_formula_string, interactions_formula_part, sep = ' + ')
interaction_formula <- as.formula(interaction_formula_string)
# also make a simplified formula
simplified_formula <- get_formula(var_of_interest = 'expression', fixed_effects = c(fixed_effects, random_effects), random_effects = NULL)
# whether to simplify
simplify <- T

# now plot each of these
for (confinement_i in 1:nrow(confinement)) {
  # extract the variant, region and gene
  variant <- confinement[confinement_i, 'variant']
  region <- confinement[confinement_i, 'region']
  gene <- confinement[confinement_i, 'gene']
  # extract the plottable table
  plot_df <- interaction_result[[paste(region, gene, variant)]]
  # model the interaction
  interaction_model <- NULL
  # depending on the family, the calls and anovas are different
  if (family == 'poisson') {
    if (simplify) {
      interaction_model <- lme4::glmer(formula = simplified_formula, data = plot_df, family = poisson)
    } else {
      # model with interaction
      interaction_model <- lme4::glmer(formula = interaction_formula, data = plot_df, family = poisson)
    }
  }
  else if (family == 'gaussian') {
    if (simplify) {
      interaction_model <- lm(formula = simplified_formula, data = plot_df)
    } else {
      # interaction model
      interaction_model <- lmerTest::lmer(formula = interaction_formula, data = plot_df)
    }
  }
  # predict using the model
  plot_df[['expression_predicted']] <- as.vector(unlist(predict(interaction_model, plot_df)))
  # make the variant a character string
  plot_df[['gt']] <- as.character(plot_df[['genotype']])
  # plot this
  p <- ggplot(data = plot_df, mapping = aes(x = region, y = expression, colour = gt)) +
    geom_point() + 
    # geom_line(mapping = aes(x = region, y = expression_predicted, colour = gt))  +
    geom_smooth(method = 'lm', formula = y~x) + 
    # geom_smooth(method = 'lm', formula = simplified_formula) + 
    scale_colour_manual(values = roycols::get_color_list(unique(plot_df[['gt']]))) + 
    xlab(paste(region, 'accessibility')) + 
    ylab(paste(gene, 'expression')) + 
    labs(colour = 'Genotype') + 
    theme(legend.title = element_text(size=14), 
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(size=12)) + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  p
}
