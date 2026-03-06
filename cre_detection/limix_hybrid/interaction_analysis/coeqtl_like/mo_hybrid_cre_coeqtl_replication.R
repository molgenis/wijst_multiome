#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_hybrid_cre_coeqtl_replication.R
# Function: perform interaction-eQTL at single-cell level with TF or ATAC as interaction terms, using the co-eQTL framework
# Example: 
# ~/start_Rscript.sh \
#   /groups/umcg-franke-scrna/tmp04/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_hybrid_cre_coeqtl_replication.R \
#   --in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/chr6-6385763-6874450/ \
#   --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/coeqtl/monocyte/lane_donor/chr6-6385763-6874450// \
#   --confinement /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement.tsv.gz \
#   --smf_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/smf.tsv.gz \
#   --covariates_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz \
#   --accessibility_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_monocyte_nonsparse_transposed.tsv.gz \
#   --fixed_effects genotype \
#   --random_effects sample_final,lane \
#   --interaction_terms expression,region \
#   --aggregate_columns sample_final,lane \
#   --barcode_column barcode_lane \
#   --genotype_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr6 \
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
# to plot
library(ggplot2)
library(roycols)


####################
# Functions        #
####################

#' get a label dict that replaces the posix safe names into printable versions
#' 
#' @returns a label dict that replaces the posix safe names into printable versions
#' label_dict_names <- get_label_dict()
get_label_dict <- function() {
  label_dict <- list()
  label_dict[['CD4_T_cells']] <- 'CD4+ T cells'
  label_dict[['CD8_T_cells']] <- 'CD8+ T cells'
  label_dict[['CD4T']] <- 'CD4+ T'
  label_dict[['CD8T']] <- 'CD8+ T'
  label_dict[['CD4_T']] <- 'CD4+ T'
  label_dict[['CD8_T']] <- 'CD8+ T'
  label_dict[['Dendritic_cells']] <- 'Dendritic cells'
  label_dict[['Endothelial_cells']] <- 'Endothelial cells'
  label_dict[['Fibroblasts']] <- 'Fibroblasts'
  label_dict[['Glia_cells']] <- 'Glia cells'
  label_dict[['Mast_cells']] <- 'MAST cells'
  label_dict[['Mature_absorptive_enterocytes']] <- 'Mature absorptive enterocytes'
  label_dict[['Mature_secretory_enterocytes']] <- 'Mature secretory enterocytes'
  label_dict[['Memory_B']] <- 'Memory B cells'
  label_dict[['Monocytes']] <- 'Monocytes'
  label_dict[['monocyte']] <- 'Monocyte'
  label_dict[['Mono']] <- 'Monocyte'
  label_dict[['Plasma_cells']] <- 'Plasma cells'
  label_dict[['Stem_cells']] <- 'Stem cells'
  label_dict[['Stromal_cells']] <- 'Stromal cells'
  label_dict[['T_others']] <- 'other T cells'
  label_dict[['Transit_amplifying_cells']] <- 'Transit amplifying cells'
  label_dict[['AI']] <- 'Actively Inflamed'
  label_dict[['NI']] <- 'Non-Inflamed'
  return(label_dict)
}


rename_labels <- function(vector_to_rename) {
  # get the labels that are present
  label_renaming <- get_label_dict()
  # now check which labels we are missing
  missing_renames <- setdiff(unique(vector_to_rename), names(label_renaming))
  # add those renames as not being renames
  for (missing_rename in missing_renames) {
    label_renaming[[missing_rename]] <- missing_rename
  }
  # now replace each value with the rename
  renamed_vector <- as.vector(unlist(label_renaming[vector_to_rename]))
  # and return that
  return(renamed_vector)
}


remap_with_label_dict <- function(vector_of_names) {
  # get the label dict
  relabels <- get_label_dict()
  # get the labels available for renaming
  labels_available <- names(relabels)
  # get the ones we cant remap
  unmappable <- setdiff(unique(vector_of_names), labels_available)
  # report on those
  if (length(unmappable) > 0) {
    print(paste('cannot remap the following names, they will be returned unchanged:', paste(unmappable, collapse = ',')))
    # and put those in our remapping list as their originals
    relabels[unmappable] <- unmappable
  }
  # now actually do the remapping
  remapped <- as.vector(unlist(relabels[vector_of_names]))
  return(remapped)
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
  color_coding_dict[['plasmablast']] <- '#DB8E00'
  color_coding_dict[['Stem cells']] <- '#66A61E'
  color_coding_dict[['Stromal cells']] <- '#8DD3C7'
  #color_coding_dict[['other T cells']] <- '#A6761D'
  color_coding_dict[['other T cells']] <- '#FF63B6'
  color_coding_dict[['T_other']] <- '#FF63B6'
  color_coding_dict[['T_other']] <- '#FF63B6'
  color_coding_dict[['Transit amplifying cells']] <- '#FF7F00'
  color_coding_dict[['Transit\namplifying cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  color_coding_dict[['unannotated']] <- 'gray'
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
gausnorm_independent_variable_matrix <- function(independent_variable_matrix, feature_id_column='feature', boxcox=F, min_value=1e-6) {
  # take the features
  features <- independent_variable_matrix[[feature_id_column]]
  # remove the feature ID
  independent_variable_matrix[[feature_id_column]] <- NULL
  # take the donor names, as they are the columns
  colnames_original <- colnames(independent_variable_matrix)
  # transpose the matrix, as we'll do this on a per-column basis
  independent_variable_matrix_t <- t(as.matrix(independent_variable_matrix))
  # transform the data
  transformed_data <- apply(independent_variable_matrix_t, 2, function(x) {
    if (is.numeric(x)) {
      if (!is.null(min_value)) {
        x[x < min_value] <- min_value
      }
      if (boxcox) {
        boxcox(x)$x.t
      } else {
        yeojohnson(x)$x.t
      }
    }
    else {
      x
    }
  })
  # make into table
  transformed_data <- do.call('cbind', transformed_data)
  # transform back and make datatable
  transformed_data <- as.data.table(t(transformed_data))
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

gausnorm_independent_variable <- function(x, boxcox=F, min_value=1e-6) {
  # initialize value 
  y <- NULL
  # only if numeric we can convert
  if (is.numeric(x) && boxcox) {
    if (!is.null(min_value)) {
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
                                   interactions=c('genotype','region'), 
                                   aggregate_columns=c('sample_final', 'lane'),
                                   accessibility_gausnorm=T, 
                                   expression_gausnorm=T, 
                                   accessibility_boxcox=F, 
                                   expression_boxcox=F
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
    covariate_columns <- setdiff(c(fixed_effects, random_effects, interactions, aggregate_columns), c('expression', 'region', 'genotype'))
    
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
      # gausnorm them if requested
      if (accessibility_gausnorm) {
        covariates_data[['region']] <- gausnorm_independent_variable(covariates_data[['region']], accessibility_boxcox)
      }
      if(expression_gausnorm) {
        covariates_data[['expression']] <- gausnorm_independent_variable(covariates_data[['expression']], expression_boxcox)
      }
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


calculate_per_sample_correlation <- function(full_variates_table, formula_string='expression~region', correlation=T, method='spearman', sample_column='sample_id', beta_variate_column='region', family='gaussian') {
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
      warning(paste('Error in correlation or regression model for sample', region, 'gene', gene, 'variant', variant, ':', e$message, '. This can happen if the model fails to converge'))
    })
    # make into df
    cor_per_sample[[sample_present]] <- data.table('sample' = c(sample_present), 'estimate' = c(estimate), 'p' = c(p), 'ncell' = c(ncell))
  }
  # merge all
  cor_all <- rbindlist(cor_per_sample, fill = T)
  return(cor_all)
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
# set debug flag
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
  make_option(c("-m", "--aggregate_columns"), type="character", default=NULL,
              help="comma separated list of columns to calculate the per-sample correlation for [default= %default]", metavar="character"),
  make_option(c("-n", "--accessibility_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on accessibility before modelling [default: %default]"), 
  make_option(c("-y", "--expression_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on expression before modelling [default: %default]"), 
  make_option(c("-z", "--correlations_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on correlations before modelling [default: %default]"), 
  make_option(c("-t", "--interaction_terms"), type="character", default=NULL,
              help="interaction to model", metavar='character'), 
  make_option(c("-g", "--genotype_loc"), type="character", default=NULL,
              help="genotype file in plink1 format, without extension", metavar="character"), 
  make_option(c("-b", "--barcode_column"), type="character", default='barcode_lane', 
              help="barcode column for metadata", metavar="character"), 
  make_option(c("-l", "--ncell_cutoff"), type="character", default=10, 
              help="number of cells required to plot and model a sample", metavar='numeric')
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
# whether to gausnorm the expression data
expression_boxcox <- F
# whether to gausnorm the accessibility/TF data
accessibility_boxcox <- F
# fixed effects string
fixed_effects_string <- NULL
# random effects string
random_effects_string <- ''
# interaction terms string
interaction_terms_string <- NULL
# aggregate columns string
aggregate_columns_string <- NULL
# covariates file
covariates_file <- NULL
# barcode column in the covariates data
barcode_column <- NULL
# number of cells required to keep a correlation for a sample
ncell_cutoff <- 10
# gausnorm the correlations
correlations_gausnorm <- T

if (debug) {
  # set all of the variables hardcoded for a testing debug run
  confinement_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement.tsv.gz'
  in_dir <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/chr6-6385763-6874450/'
  smf_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/smf.tsv.gz'
  output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/coeqtl/monocyte/lane_donor/chr6-6385763-6874450/'
  expression_file <- 'expression.tsv.gz'
  accessibility_file <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_monocyte_nonsparse_transposed.tsv.gz'
  genotype_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr6'
  covariates_file <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
  fixed_effects_string <- 'genotype'
  random_effects_string <- 'sample_final,lane'
  interaction_terms_string <- 'region,expression'
  barcode_column <- 'barcode_lane'
  aggregate_columns_string <- 'sample_final,lane'
  ncell_cutoff <- 10
  correlations_gausnorm <- T
  
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
  if (is.null(opt[['aggregate_columns']])) {
    error("m/--aggregate_columns is an obligatory parameter")
  } else {
    aggregate_columns_string <- opt[['aggregate_columns']]
  }
  if (is.null(opt[['random_effects']])) {
    random_effects_string <- ''
  } else {
    random_effects_string <- opt[['random_effects']]
  }
  # parameters that have a sane default
  random_effects_string <- opt[['random_effects']]
  barcode_column <- opt[['barcode_column']]
  expression_file <- opt[['expression_file']]
  accessibility_file <- opt[['accessibility_file']]
  covariates_file <- opt[['covariates_file']]
  accessibility_gausnorm <- opt[['accessibility_gausnorm']]
  expression_gausnorm <- opt[['expression_gausnorm']]
  correlations_gausnorm <- opt[['correlations_gausnorm']]
  ncell_cutoff <- opt[['ncell_cutoff']]
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
        # split random effects
        aggregate_columns <- c()
        if (!is.null(aggregate_columns_string) && !is.na(aggregate_columns_string) && aggregate_columns_string != '') {
          aggregate_columns <- strsplit(aggregate_columns_string, ',')[[1]]
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
            if (expression_boxcox) {
              message('Yeo-Johnson gausnorm on expression data...')
              # expression_data <- gausnorm_independent_variable_matrix(independent_variable_matrix = expression_data, feature_id_column = 'gene', boxcox = T)
            } else {
              message('Yeo-Johnson gausnorm on expression data...')
              # expression_data <- gausnorm_independent_variable_matrix(independent_variable_matrix = expression_data, feature_id_column = 'gene')
            }
          }
          if (accessibility_gausnorm) {
            if (accessibility_boxcox) {
              message('Yeo-Johnson gausnorm on accessibility/TF data...')
              # accessibility_data <- gausnorm_independent_variable_matrix(independent_variable_matrix = accessibility_data, feature_id_column = 'region', boxcox = T)
            } else {
              message('Yeo-Johnson gausnorm on accessibility/TF data...')
              # accessibility_data <- gausnorm_independent_variable_matrix(independent_variable_matrix = accessibility_data, feature_id_column = 'region')
            }
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
            interactions = interactions, 
            aggregate_columns = aggregate_columns,
            accessibility_gausnorm = accessibility_gausnorm, 
            expression_gausnorm = expression_gausnorm, 
            accessibility_boxcox = accessibility_boxcox, 
            expression_boxcox = expression_boxcox
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

# store the plots
cor_gt_plots <- list()
# the inputs
cor_gt_dfs <- list()
# and the results
cor_gt_results <- list()

if (!is.null(interaction_result)) {
  # get co-eQTL style plots
  for (confinement_i in 1:nrow(confinement)) {
    # extract the variant, region and gene
    variant <- as.vector(unlist(confinement[confinement_i, 'variant']))
    region <- as.vector(unlist(confinement[confinement_i, 'region']))
    gene <- as.vector(unlist(confinement[confinement_i, 'gene']))
    # paste together the naming
    confinement_name <- paste(region, gene, variant)
    # check if this combination is in the results
    if (confinement_name %in% names(interaction_result)) {
      # extract the plottable table
      plot_df <- interaction_result[[confinement_name]]
      # add an aggregated column for the sample, by pasting the aggregate columns together
      plot_df[['aggregated_sample']] <- apply(plot_df[, aggregate_columns, drop = F], 1, function(x) paste(x, collapse = '_'))
      # add an aggregated column for the sample, by pasting the aggregate columns together
      covariates_data[['aggregated_sample']] <- apply(covariates_data[, aggregate_columns, drop = F], 1, function(x) paste(x, collapse = '_'))
      # get the per-sample plot
      per_sample_df <- calculate_per_sample_correlation(plot_df, sample_column = 'aggregated_sample', formula_string = paste(interactions, sep = '~', collapse = '~'))
      # store original estimate
      per_sample_df[['estimate_raw']] <- per_sample_df[['estimate']]
      # gausnorm if requested
      if (correlations_gausnorm) {
        per_sample_df[['estimate']] <- gausnorm_independent_variable(per_sample_df[['estimate']])
      }
      # take the unique sets of the covariates from the plot df, to add this to the per sample df
      unique_covariate_columns <- unique(c(fixed_effects, random_effects, interactions, aggregate_columns))
      # but remove region and expression
      unique_covariate_columns <- setdiff(unique_covariate_columns, c('region', 'expression'))
      # subset the plot df to these columns and the sample column, and take unique rows
      plot_df_unique_covariates <- unique(plot_df[, c('aggregated_sample', unique_covariate_columns)])
      # then add this to the plot df
      per_sample_df <- merge(per_sample_df, plot_df_unique_covariates, by.x = 'sample', by.y = 'aggregated_sample', all.x = T)
      # make character string
      per_sample_df[['gt']] <- as.character(per_sample_df[['genotype']])
      # create a formula
      base_formula <- get_formula(var_of_interest = 'estimate', fixed_effects = c(fixed_effects, 'ncell'), random_effects = random_effects)
      # initialize variable
      lm_gt_to_cor_table_base <- data.table('variant' = c(variant), 'region' = c(region), 'gene' = c(gene))
      lm_gt_to_cor_table <- NULL
      # try to do modelling
      tryCatch({
        # fit model with the genotype and the correlation
        lm_gt_to_cor <- NULL
        if (!is.null(random_effects) && length(random_effects) > 0) {
          # use glm if we have random effects
          lm_gt_to_cor <- lmerTest::lmer(formula = base_formula, data = per_sample_df[per_sample_df$ncell >= ncell_cutoff, ])
        } else {
          # or simple model if here are not
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
        warning(paste('Error in model fitting', region, gene, variant, ':', e$message))
        # make empty table
        lm_gt_to_cor_table <- lm_gt_to_cor_table_base
      })
      # add the number of samples
      lm_gt_to_cor_table[['nsample']] <- nrow(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ])
      # and the distribution of cells
      lm_gt_to_cor_table[['ncell']] <- paste(as.character(min(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])),
                                             as.character(quantile(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])[['25%']]),
                                             as.character(quantile(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])[['50%']]),
                                             as.character(quantile(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])[['75%']]),
                                             as.character(max(per_sample_df[per_sample_df$ncell >= ncell_cutoff, ][['ncell']])),
                                             sep = ';'
      )
      # get rounded significance
      significance <- 'NA'
      if (is.numeric(lm_gt_to_cor_p)) {
        significance <- as.character(round(lm_gt_to_cor_p, digits = 5))
      }
      # plot both of them
      p_cor <- ggplot(data = per_sample_df[per_sample_df$ncell >= ncell_cutoff, ], mapping = aes(x = gt, y = estimate_raw, fill = gt)) + 
        geom_boxplot(outlier.shape = NA) + 
        # geom_point() +
        # and add jitter
        geom_jitter(size = 0.5, alpha = 0.5, data = per_sample_df[per_sample_df$ncell >= ncell_cutoff, ], mapping = aes(x = gt, y = estimate_raw, colour = ncell)) +
        # geom_jitter(size = 0.5, alpha = 0.5) + 
        scale_fill_manual(values = roycols::get_color_list(unique(per_sample_df[['gt']]))) + 
        # and colour of ncell
        scale_colour_gradient2(low='blue', mid = 'white', high='red', midpoint = max(per_sample_df$ncell/2)) + 
        xlab(paste('genotype')) + 
        ylab(paste(gene, 'region ~ expression')) + 
        labs(fill = 'Genotype', colour = 'Ncell') + 
        ggtitle(paste('cor', variant, region, gene, 'p < ', significance)) +
        theme(legend.title = element_text(size=14), 
              legend.text = element_text(size=12),
              axis.title.x = element_text(size=14),
              axis.title.y = element_text(size=14),
              axis.text.y = element_text(size=12),
              axis.text.x = element_text(size=12),
              strip.text.x = element_text(size=12)) + 
        theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
      # place in table
      cor_gt_plots[[paste(variant, region, gene, 'cor', sep = '_')]] <- p_cor
      cor_gt_dfs[[paste(variant, region, gene, 'cor', sep = '_')]] <- per_sample_df
      cor_gt_results[[paste(variant, region, gene, 'cor', sep = '_')]] <- lm_gt_to_cor_table
    }
  }
}

# make output directory if it doesn't exist
if (!dir.exists(output_loc)) {
  dir.create(output_loc, recursive = T)
}

# now do each combination
for (comb in names(cor_gt_results)) {
  # extract the result table
  result_table <- cor_gt_results[[comb]]
  # check if the result table has more than 0 rows
  if ((!is.null(result_table)) && (!is.null(dim(result_table))) && (length(dim(result_table)) == 2) && (nrow(result_table) > 0)) {
    # check if the genotype association was significant
    if (result_table[['genotype_p']][1] < 0.05) {
      # if so, save the plot
      p_to_plot <- cor_gt_plots[[comb]]
      # make a safer name
      p_name <- gsub('\\:|/', '.', comb)
      # save the plot
      ggsave(paste0(output_loc, '/', p_name, '.pdf'), plot = p_to_plot, width = 6, height = 6)
    }
  }
}

# check if there were any results
if (!is.null(cor_gt_results) && length(cor_gt_results) > 0) {
  # save the results table as well
  cor_gt_results_table <- rbindlist(cor_gt_results, fill = T)
  # extract the last part of the folder
  chunk_name <- basename(in_dir)
  # save this
  cor_gt_results_table[['chunk']] <- rep(chunk_name, times = nrow(cor_gt_results_table))
  # write result
  write.table(cor_gt_results_table, output_loc_full, sep = '\t', row.names = F, col.names = T, quote = F)
  # make a checksum
  mdfiver::create_sha256_for_file(tsv_output_loc_full)
} else {
  # write an empty result
  write_empty_result(tsv_output_loc_full)
}

