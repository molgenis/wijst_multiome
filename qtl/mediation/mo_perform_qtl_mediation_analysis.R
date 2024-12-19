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
  # set model
  indirect_model1 <- NULL
  indirect_model2 <- NULL
  # build random
  if (sum(grepl('\\|', as.character(formula_indirect1))) > 0) {
    ## First-part indirect effect model:
    indirect_model1 <- do.call(what = 'glmer', list(formula = formula_indirect1, data = mediation_table, family = gaussian(link = "identity")))
    ## Second-part direct + indirect effect model:
    indirect_model2 <- do.call(what = 'glmer', list(formula = formula_mediation, data = mediation_table, family = gaussian(link = "identity")))
  }
  # or only linear
  else {
    ## First-part indirect effect model:
    indirect_model1 <- do.call(what = 'glm', list(formula = formula_indirect1, data = mediation_table, family = gaussian(link = "identity")))
    ## Second-part direct + indirect effect model:
    indirect_model2 <- do.call(what = 'glm', list(formula = formula_mediation, data = mediation_table, family = gaussian(link = "identity")))
  }

  ## Mediation analysis with 1000 simulations
  mediation_results <- mediate(indirect_model1, indirect_model2, treat = 'genotype', mediator = 'accessibility', boot = TRUE, sims = sims)

  return(mediation_results)
}


mediate_all_effects <- function(accessibility, expression, metadata, genotypes, confinement, form_indirect1, form_indirect2) {
  # save result per set
  res_per_set <- apply(confinement, 1, function(x) {
    # get from the row
    variant <- x[1]
    region <- x[2]
    gene <- x[3]
    # get the mediation table
    medation_table <- create_mediation_table(accessibility_table = accessibility, expression_table = expression, metadata = metadata, genotypes = genotypes, variant = variant, region = region, gene = gene)
    # do the mediation analysis
    med_result <- perform_mediation_analysis(mediation_table = medation_table, formula_indirect1 = form_indirect1, formula_mediation = form_indirect2)
    return(med_result)
  })
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
  form_indirect1 <- get_formula('accessibility', c('genotype'), c())
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


do_debug <- function() {
  options_debug <- list()
  options_debug[['eqtl_file']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/L1/UT/B.qtlInput.txt.gz'
  options_debug[['caqtl_file']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/B.qtlInput.txt.gz'
  options_debug[['genotype_file']] <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr7'
  options_debug[['confinement_list']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/confinements/B.confinement.tsv.gz'
  options_debug[['metadata']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/metadata/B.metadata.tsv.gz'
  options_debug[['fixed_effects']] <- 'RNA_UT_PC1,RNA_UT_PC2,RNA_UT_PC3,RNA_UT_PC4,RNA_UT_PC5,RNA_UT_PC6,RNA_UT_PC7,RNA_UT_PC8,RNA_UT_PC9,RNA_UT_PC10'
  #options_debug[['random_effects']] <- 'lane,donor
  options_debug[['random_effects']] <- ''
  options_debug[['out']] <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/output/atac_to_expression/'
  # get the inputs
  inputs <- load_inputs(options_debug)
  # filter the inputs
  inputs <- filter_inputs(inputs)
  return(inputs)
}


####################
# Settings        #
####################

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)

# size of chunks to normalize
chunk_size <- 5000000
registerDoParallel(cores = 8)


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
              help="output location of analysis [default= %default]", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
