#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_perform_qtl_mediation_analysis.R
# Function: perform mediation analysis of eQTLs through caQTLs
############################################################################################################################

####################
# libraries        #
####################


library(data.table)
library(snpStats)
library(lme4)
library(MEPS)
library(bda)
library(mediation)
library(mdfiver)
library(parallel)
library(doParallel)
library(pbapply)

####################
# Functions        #
####################


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
  }, cl = cl)
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
  samples_genotypes <- intersect(rownames(genotypes$genotypes), sample_to_donor[['donor']])
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
  genotypes$genotypes <- genotypes$genotypes[sample_to_donor[['donor']], ]
  genotypes$fam <-  genotypes$fam[rownames(genotypes$genotypes), ]
  
  # filter expression
  expression <- expression[, .SD, .SDcols=c('feature', sample_to_donor[['sample']])]
  # filter accessibility
  accessibility <- accessibility[, .SD, .SDcols=c('feature', sample_to_donor[['sample']])]
  # and finally metadata
  metadata <- metadata[metadata[['sample']] %in% sample_to_donor[['sample']], ]
  
  # put back into the list
  inputs[['expression']] <- expression
  inputs[['accessibility']] <- accessibility
  inputs[['genotypes']] <- genotypes_data
  inputs[['metadata']] <- metadata
  
  return(inputs)
}


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
    'out' = options[['out']]
  )
  return(inputs)
}


run_full_analysis <- function(options, verbose=T) {
  # start multicore
  cl <- makeCluster(options[['threads']])
  registerDoParallel(cl)
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
    cluster = cl)
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
  make_option(c("-t", "--threads"), type="numeric", default=1,
              help="number of threads to use [default= %default]", metavar="numeric")
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
