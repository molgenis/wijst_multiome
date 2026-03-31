#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_sc_tfqtl.R
# Function: perform transcription-factor QTL analysis at single-cell level
# Example: 
# ~/start_Rscript.sh \
#   /groups/umcg-franke-scrna/tmp02/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_sc_tfqtl.R \
#   --in /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/all/chr12/ \
#   --out /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/tfqtl/limix_sc/output/tf/all/chr12/ \
#   --confinement /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/tfqtl/limix_sc/confinement/mo_var_tf.tsv.gz \
#   --smf_loc /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/all/smf.tsv.gz \
#   --covariates_file /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz \
#   --tf_file /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_all_nonsparse_transposed.tsv.gz \
#   --fixed_effects genotype,nCount_RNA \
#   --random_effects sample_final,lane \
#   --barcode_column barcode_lane \
#   --genotype_loc /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr12 \
#   --tf_gausnorm
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


do_association_analysis <- function(tf_data, 
                                    genotype_data, 
                                    smf, 
                                    confinement,
                                    covariates_data=NULL, 
                                    fixed_effects=c('lane','region','genotype'), 
                                    random_effects=c('sample_final'), 
                                    family = 'gaussian', 
                                    tf_gausnorm = T, 
                                    tf_boxcox = F) {
  # create formula
  base_formula <- get_formula(var_of_interest = 'tf', fixed_effects = fixed_effects, random_effects = random_effects)
  base_formula_string <- (Reduce(paste, deparse(base_formula)))
  message(paste('Using base formula:', base_formula_string, ''))
  # conver to dataframes
  tf_data <- data.frame(tf_data)
  smf <- data.frame(smf)
  confinement <- data.frame(confinement)
  if (!is.null(covariates_data)) {
    covariates_data <- data.frame(covariates_data)
    # get which things we need from that dataframe
    covariate_columns <- setdiff(c(fixed_effects, random_effects), c('tf', 'genotype'))
    # and subset the covariates data to that
    covariates_data <- covariates_data[, c('cell', covariate_columns)]
  } else {
    # otherwise make a dummy where we'll add enverything
    covariates_data <- data.frame(
      'cell' = setdiff(colnames(tf_data, c('cell'))),
      'dummy' = rep(T, times = (ncol(tf_data) - 1))
    )
  }
  # put all results in a list
  res_per_comparison <- list()
  # get unique TFs
  unique_tfs <- unique(confinement[['tf']])
  # count how many
  n_unique_tfs <- length(unique_tfs)
  # store warnings
  warnings_encountered <- list()
  # make progress bar
  pb = txtProgressBar(min = 0, max = n_unique_tfs, initial = 0) 
  # check each gene
  for (tf_i in 1 : n_unique_tfs) {
    # update progress bar
    setTxtProgressBar(pb,tf_i)
    # get the tf
    tf <- unique_tfs[tf_i]
    # extract the gene
    tf_values <- as.vector(unlist(tf_data[tf_data[['tf']] == tf, 2:ncol(tf_data)]))
    covariates_data[['tf']] <- tf_values
    if(tf_gausnorm) {
      covariates_data[['tf']] <- gausnorm_independent_variable(covariates_data[['tf']], tf_boxcox)
    }
    # get the variants for this region-gene combination
    variants_tf <- unique(confinement[confinement[['tf']] == tf, ][['variant']])
    # check each variant
    for (variant in variants_tf) {
      # extract genotypes
      genotype <- genotype_data$genotypes[smf[['participant']], variant]
      # then to numeric
      genotype_numeric <- as.vector(as(genotype, 'numeric'))
      # add the genotype
      covariates_data[['genotype']] <- genotype_numeric
      # keep only complete cases
      covariates_data_complete <- covariates_data[complete.cases(covariates_data), ]
      # and only finite values
      covariates_data_complete <- covariates_data_complete[is.finite(covariates_data_complete[['tf']]) &  is.finite(covariates_data_complete[['genotype']]), ]
      # check if there is any data left
      if (nrow(covariates_data_complete) > 0) {
        # initialize variables
        base_model_df <- data.frame(
          'variant' = c(variant), 
          'tf' = c(tf)
        )
        # try to do 
        tryCatch({
          # depending on the family, the calls and anovas are different
          if (family == 'poisson') {
            # model without interaction
            base_model <- lme4::glmer(formula = base_formula, data = covariates_data_complete, family = poisson)
          }
          else if (family == 'gaussian') {
            # base model
            base_model <- lmerTest::lmer(formula = base_formula, data = covariates_data_complete)
          }
          else {
            stop(paste0('unknown family ', family, ', only valid families are gaussian and poisson'))
          }
          
        }, error = function(e) {
           # warning(paste('Error in model fitting', tf, variant, ':', e$message, '. This can happen if the model fails to converge'))
          # put warning in list
          warnings_encountered[[paste(variant, tf, 'model')]] <- paste('Error in model fitting', tf, variant, ':', e$message, '. This can happen if the model fails to converge')
        })
        if(!is.null(base_model)) {
          # if we at least have a base model, we can still convert that to a df
          base_model_df <- cbind(base_model_df, model_to_row(base_model))
        } else {
          # if the model failed to fit, we still don't have anything
          
        }
        # add the family used
        base_model_df[['family']] <- family
        # keep the number of cells we have
        base_model_df[['ncell']] <- nrow(covariates_data_complete)
        # and participants
        base_model_df[['nparticipant']] <- length(unique(smf[smf[['cell']] %in% covariates_data_complete[['cell']], ][['participant']]))
        # store result
        res_per_comparison[[paste(tf, variant)]] <- data.table(base_model_df)
      } else {
        # warning(paste('No data left for tf', tf, 'variant', variant, 'after complete cases. Skipping this combination.'))
        # put warning in list
        warnings_encountered[[paste(variant, tf, 'data')]] <- paste('No data left for tf', tf, 'variant', variant, 'after complete cases. Skipping this combination.')
      }
    }
  }
  # talk about which warnings were encountered
  if (length(names(warnings_encountered)) > 0) {
    warning('the following warnings were encountered:')
    for (warning_name in names(warnings_encountered)) {
      print(warnings_encountered[[warning_name]])
    }
  }
  # close progress bar
  close(pb)
  # merge all results
  res_all <- rbindlist(res_per_comparison, fill = T)
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
              help="confinement file of variant-TF combinations to test", metavar="character"),
  make_option(c("-s", "--smf_loc"), type="character", default=NULL, 
              help="sample mapping file", metavar="character"), 
  make_option(c("-t", "--tf_file"), type="character", default='expression.tsv.gz', 
              help="transcription factor filename for chunk", metavar="character"), 
  make_option(c("-v", "--covariates_file"), type="character", default=NULL, 
              help="accessibility filename for chunk", metavar="character"), 
  make_option(c("-f", "--fixed_effects"), type="character", default=NULL,
              help="comman separated list of fixed effects to correct for [default= %default]", metavar="character"),
  make_option(c("-r", "--random_effects"), type="character", default=NULL,
              help="comma separated list of random effects to correct for [default= %default]", metavar="character"),
  make_option(c("-n", "--tf_gausnorm"), action="store_true", default=FALSE,
              help="Apply the yeo-johnson transformation on TF activity before modelling [default: %default]"), 
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
tf_file <- NULL
# genotype location
genotype_loc <- NULL
# whether to gausnorm the expression data
tf_gausnorm <- T
# fixed effects string
fixed_effects_string <- NULL
# random effects string
random_effects_string <- NULL
# covariates file
covariates_file <- NULL
# barcode column in the covariates data
barcode_column <- NULL

# we'll keep this the same
tf_boxcox <- F

if (debug) {
  # set all of the variables hardcoded for a testing debug run
  confinement_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/tfqtl/limix_sc/confinement/mo_var_tf.tsv.gz'
  in_dir <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/all/chr22/'
  smf_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/all/smf.tsv.gz'
  output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/tfqtl/limix_sc/output/tf/all/chr22/'
  tf_file <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_all_nonsparse_transposed.tsv.gz'
  genotype_loc <- '/groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr22'
  covariates_file <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz'
  fixed_effects_string <- 'nCount_SCT,genotype'
  random_effects_string <- 'sample_final,lane'
  barcode_column <- 'barcode_lane'
  tf_gausnorm <- T
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
  # parameters that have a sane default
  random_effects_string <- opt[['random_effects']]
  barcode_column <- opt[['barcode_column']]
  tf_file <- opt[['tf_file']]
  covariates_file <- opt[['covariates_file']]
  tf_gausnorm <- opt[['tf_gausnorm']]
}
# make the full path to the expression data
full_tf_path <- NULL
# depending on if it is an absolute path, we do things differently
if (startsWith(tf_file, '/')) {
  full_tf_path <- tf_file
} else {
  full_tf_path <- paste(in_dir, tf_file, sep = '/')
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
colnames(confinement) <- c('variant', 'tf')

# initialize variables
tf_data <- NULL
# check if there is expression data
if (length(count.fields(full_tf_path)) > 1) {
  # read the expression file
  # expression_data <- read.table(full_exp_path, header = T, sep = '\t', check.names = F, row.names = 1)
  tf_data <- fread(full_tf_path, header = T, sep = '\t', check.names = F, skip = 1)
  # read the header
  tf_data_header_line <- readLines(full_tf_path, n = 1)
  # split by sep
  tf_data_header <- strsplit(tf_data_header_line, '\t')[[1]]
  # add this header
  if (length(tf_data_header) == ncol(tf_data)) {
    colnames(tf_data) <- tf_data_header
  } else {
    # otherwise we need an extra column
    colnames(tf_data) <- c('tf', tf_data_header)
  }
  # also set the first column name so we can refer to it later
  colnames(tf_data)[[1]] <- 'tf'
} else {
  warning('no data fields for TF data, will do no further work')
  # to avoid further nesting, we'll make a dummy entry that makes it so that we dont continue further
  tf_data <- data.table('tf' = c())
}
tf_data_confined <- tf_data[!is.na(tf_data[['tf']]) & tf_data[['tf']] %in% confinement[['tf']], ]

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

# check if we have any data left
if (nrow(tf_data) > 0) {
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
  confinement <- confinement[confinement[['tf']] %in% tf_data_confined[['tf']], ]
  # check if we have any data left
  if (nrow(confinement) > 0) {
    # read covariates data
    covariates_data <- NULL
    if (!is.null(full_covariates_path)) {
      covariates_data <- fread(full_covariates_path, header = T, sep = '\t')
      # the first column should be the cell
      covariates_data <- cbind(data.frame('cell' = covariates_data[[barcode_column]]), covariates_data)
      # intersect this
      intersecting_cells <- intersect(colnames(tf_data), covariates_data[['cell']])
      # subset the covariates
      covariates_data <- covariates_data[covariates_data[['cell']] %in% intersecting_cells, ]
      # get the expression and accessibility/TF columns again
      tf_columns <- c('tf', intersecting_cells)
      # and subset
      tf_data <- tf_data[, ..tf_columns]
    } else {
      # if we have no covariates, make a dummy table
      covariates_data <- data.table(
        'cell' = setdiff(colnames(tf_data, c('cell'))),
        'dummy' = rep(T, times = (ncol(tf_data) - 1))
      )
    }
    # read smf
    smf <- fread(smf_loc, header = T, sep = '\t')
    # harmonize names
    colnames(smf) <- c('participant', 'cell')
    # intersect the smf with the expression data and accessibility/TF data
    intersecting_cells <- intersect(smf[['cell']], colnames(tf_data))
    # make columns again
    tf_columns <- c('tf', intersecting_cells)
    # and use this order
    tf_data <- tf_data[, ..tf_columns]
    covariates_data <- covariates_data[match(intersecting_cells, covariates_data[['cell']]), ]
    smf <- smf[match(intersecting_cells, smf[['cell']])]
    # perform the analysis
    message('Starting analysis..')
    # into a variable
    interaction_result <- do_association_analysis(
      tf_data = tf_data, 
      genotype_data = genotypes, 
      smf = smf, 
      confinement = confinement,
      covariates_data = covariates_data, 
      fixed_effects = fixed_effects, 
      random_effects = random_effects, 
      tf_gausnorm = tf_gausnorm, 
      tf_boxcox = tf_boxcox
    )
    # extract the last part of the folder
    chunk_name <- basename(in_dir)
    # if there were any results, we'll write one
    if ((!is.null(interaction_result)) && (!is.null(nrow(interaction_result))) && (nrow(interaction_result) > 0)) {
      # add the chunk as a column
      interaction_result[['chunk']] <- rep(chunk_name, times = nrow(interaction_result))
      # write result
      write.table(interaction_result, output_loc_full, sep = '\t', row.names = F, col.names = T, quote = F)
      # make a checksum
      mdfiver::create_sha256_for_file(tsv_output_loc_full)
    } else {
      # just make the result null again
      interaction_result <- NULL
    }
  } else {
    message('No varian-TF pairs left after intersecting confinement with data. No more work to be done')
  }
} else {
  message('No genes left after filtering confinement. No more work to be done')
}

# if the interaction result is still null, we didn't end up doing anything
if (is.null(interaction_result)) {
  # so we'll store an empty file
  write_empty_result(tsv_output_loc_full)
}
