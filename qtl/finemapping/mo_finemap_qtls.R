#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_finemap_qtls.R
# Function: perform statistical finemapping on the QTL results
# Example: 
# Rscript ~/mo_finemap_qtls.R \
#   --qtl_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/nominal_condition/L1/monocyte/inflammation_final/iqtl_results_all.txt.gz \
#   --genotype_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/genotype/imputed_hg38_all_anc \
#   --significance_column empirical_feature_p_value \
#   --significance_cutoff 0 \
#   --output_rds /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/nominal_condition/L1/monocyte/inflammation_final/iqtl_results_all_finemapped_genotype.rds \
#   --num_threads 2 \
#   --variant_column snp_id \
#   --feature_column feature_id \
#   --slope_column beta_SNP \
#   --se_column beta_se_SNP \
#   --n_sample 318
#
############################################################################################################################


####################
# libraries        #
####################

library(r2r)
library(ReigenMT)
library(data.table)
library(snpStats)
library(corpcor)
library(susieR)
library(Rfast) # not in container
library(parallel)
library(optparse)


####################
# Functions        #
####################

#' Get Genotype Correlations
#'
#' This function calculates the correlation matrix for a given set of genotypes.
#'
#' @param genotypes A matrix, data.table, data.frame, or list containing genotype data. If the data is in PLINK format, it should be a list with a 'genotypes' element of class 'SnpMatrix'.
#' @param use_covshrink Boolean, whether or not to use a shrank correlation matrix
#' @return A correlation matrix.
#' @details The function transposes the genotype data if it is in a numeric format. If the data is in PLINK format, it converts the 'SnpMatrix' to a numeric matrix. It then replaces any missing values with the mean of the respective columns and calculates the correlation matrix using the `corpcor` package.
#' @examples
#' \dontrun{
#' # Example with a matrix
#' genotypes_matrix <- matrix(rnorm(100), nrow=10)
#' get_genotype_correlations(genotypes_matrix)
#'
#' # Example with a data.table
#' library(data.table)
#' genotypes_dt <- data.table(matrix(rnorm(100), nrow=10))
#' get_genotype_correlations(genotypes_dt)
#'
#' # Example with a PLINK format list
#' library(snpStats)
#' genotypes_plink <- list(genotypes = as(SnpMatrix(matrix(rbinom(100, 2, 0.5), nrow=10)), "SnpMatrix"))
#' get_genotype_correlations(genotypes_plink)
#' }
#' @import data.table
#' @import corpcor
#' @importFrom stats cov2cor
get_genotype_correlations <- function(genotypes, use_covshrink=F) {
  # we need a transposed matrix to calculate correlations
  genotypes_t <- NULL
  # Transpose the genotype data if in a numeric format
  if (is.matrix(genotypes) | is.data.table(genotypes) | is.data.frame(genotypes)) {
    genotypes_t <- data.table::data.table(t(genotypes))
  }
  # or if in plink format
  else if(is.list(genotypes) & class(genotypes$genotypes) == 'SnpMatrix') {
    genotypes_t <- data.table::data.table(as(genotypes$genotypes, "numeric"))
  }
  # we'll assume is a matrix-compatible format
  else {
    warning(paste('not recognizing genotype format, assuming matrix-compatible'))
    genotypes_t <- data.table::data.table(t(genotypes))
  }
  # make into double
  genotypes_t[, (names(genotypes_t)) := lapply(.SD, function(x){as.double(x)})]
  # get the correlation matrix
  cor_shrink <- NULL
  # use covariance shrinkage method
  if (use_covshrink) {
    # get the means of each column
    var_means <- colMeans(genotypes_t, na.rm = T)
    # make var means into list
    var_means_list <- as.list(var_means)
    names(var_means_list) <- colnames(genotypes_t)
    # replace all the NAs
    for (col in names(var_means_list)) {
      setnafill(genotypes_t, type=c("const","locf","nocb"), fill=var_means_list[[col]], cols=col)
    }
    # Perform the fit using the corpcor package
    fitted <- corpcor::cov.shrink(genotypes_t)
    # Extract the alpha (shrinkage intensity)
    alpha <- attributes(fitted)$lambda
    cor_shrink <- stats::cov2cor(fitted)
  }
  # or canonical pearson correlation that is present in R
  else {
    # calculate correlations
    cor_shrink <- cor(genotypes_t)
    # and set the dimension names
    rownames(cor_shrink) <- colnames(genotypes_t)
    colnames(cor_shrink) <- colnames(genotypes_t)
  }
  return(cor_shrink)
}


#' Finemap Feature
#'
#' This function performs fine-mapping for a given feature using summary statistics and genotype data.
#'
#' @param results_feature A data.frame or data.table containing the summary statistics for the feature.
#' @param genotypes A matrix, data.table, data.frame, or list containing genotype data.
#' @param variant_column A character string specifying the column name for variant IDs in `results_feature`. Default is 'variant_id'.
#' @param slope_column A character string specifying the column name for effect sizes in `results_feature`. Default is 'slope'.
#' @param se_column A character string specifying the column name for standard errors in `results_feature`. Default is 'slope_se'.
#' @param L An integer specifying the maximum number of causal variants to consider. Default is 10.
#' @param estimate_residual_variance A logical value indicating whether to estimate residual variance. Default is TRUE.
#' @param n_sample_column A character string specifying the column name for the number of samples in `results_feature`. Default is 'n_sample'.
#' @param n_sample A number specifying the number of sample. Default is 'NULL'.
#' @param initial_iter the initial number of iterations to try for finemapping
#' @param n_retries the max number of retries to try before giving up on finemapping
#' @return A list containing the fine-mapped results.
#' @details The function subsets the genotype data to the variants of interest, calculates the correlation matrix for these variants, and performs fine-mapping using the `susie_rss` function.
#' @examples
#' \dontrun{
#' # Example with summary statistics and genotype data
#' results_feature <- data.table(variant_id = c("rs1", "rs2"), slope = c(0.1, 0.2), slope_se = c(0.01, 0.02))
#' genotypes <- matrix(rnorm(200), nrow=20)
#' finemap_feature(results_feature, genotypes)
#' }
#' @import data.table
#' @import ReigenMT
#' @import susieR
finemap_feature <- function(results_feature, genotypes, variant_column='variant_id', slope_column='slope', se_column='slope_se', L=10, estimate_residual_variance=T, n_sample_column=NULL, n_sample=NULL, initial_iter=100, n_retries=3) {
  # get the n
  n <- NA
  if (!is.null(n_sample_column)) {
    # get that n
    n = round(max(results_feature[[n_sample_column]]))
  }
  else if(!is.null(n_sample)){
    # grab the max n
    n <- n_sample
  }
  # get the variants for this feature
  variants_feature <- results_feature[[variant_column]]
  # subset the genotype data to those features
  genotypes_variants <- ReigenMT::subset_genotypes(genotypes, variants_feature)
  # get the correlation matrix for those variants
  genotype_correlations <- get_genotype_correlations(genotypes_variants)
  # if there is only one variant, we need to not order, because it will be made into a single value
  if (nrow(genotype_correlations) > 1) {
    # make sure it is in the same order as the variants
    genotype_correlations <- genotype_correlations[variants_feature, variants_feature]
  }
  # perform finemapping
  finemapped_feature <- susie_rss(bhat = results_feature[[slope_column]], shat = results_feature[[se_column]], n = n, R = genotype_correlations, L = L, estimate_residual_variance = estimate_residual_variance)
  # let's see how often we retried
  n_retried <- 0
  # and how many iterations we used
  n_iterations_used <- initial_iter
  # if we didnt converge, let's keep trying
  if (!finemapped_feature[['susie_rss']]$converged) {
    # set up our converge parameter
    converged <- F
    # and keep trying until we run out of retries or we converge
    while(!(converged) & n_retried < n_retries) {
      # get the number of times we tried
      n_tried <- n_retried + 1
      # and use that to get the new number of iterations
      n_iterations_used <- (2 ^ (n_tried)) * initial_iter
      # then rerun
      finemapped_feature <- susie_rss(bhat = results_feature[[slope_column]], shat = results_feature[[se_column]], n = n, R = genotype_correlations, L = L, estimate_residual_variance = estimate_residual_variance, max_iter = n_iterations_used)
      # get whether we converged
      converged <- finemapped_feature[['susie_rss']]$converged
      # increase the number of times we retried
      n_retried <- n_retried + 1
    }
  }
  # put into a list
  finemapped_feature <- list('susie_rss' = finemapped_feature)
  # get the lambda
  finemapped_lambda <- NA
  if (finemapped_feature[['susie_rss']]$converged & !is.na(n)) {
    finemapped_lambda <- estimate_s_rss(z = results_feature[[slope_column]] / results_feature[[se_column]],
                                        R = genotype_correlations, 
                                        n = n)
  }
  # add to the finemapping
  finemapped_feature[['lambda']] <- finemapped_lambda
  # and the number of retries
  finemapped_feature[['n_retries']] <- n_retried
  # and the eventual number of iterations
  finemapped_feature[['n_iterations']] <- n_iterations_used
  return(finemapped_feature)
}


#' Finemap Features in Parallel
#'
#' This function performs fine-mapping for multiple features in parallel using summary statistics and genotype data.
#'
#' @param qtl_results A data.frame or data.table containing the QTL results.
#' @param genotypes A matrix, data.table, data.frame, or list containing genotype data.
#' @param feature_column A character string specifying the column name for feature IDs in `qtl_results`. Default is 'phenotype_id'.
#' @param significance_columns A named list specifying the significance thresholds for filtering features. Default is list('qval' = 0.05).
#' @param variant_column A character string specifying the column name for variant IDs in `qtl_results`. Default is 'variant_id'.
#' @param slope_column A character string specifying the column name for effect sizes in `qtl_results`. Default is 'slope'.
#' @param se_column A character string specifying the column name for standard errors in `qtl_results`. Default is 'slope_se'.
#' @param L An integer specifying the maximum number of causal variants to consider. Default is 10.
#' @param estimate_residual_variance A logical value indicating whether to estimate residual variance. Default is TRUE.
#' @param nthreads An integer specifying the number of threads to use for parallel processing. Default is the number of cores detected by `detectCores()`.
#' @return A list of fine-mapped results for each feature.
#' @details The function filters the QTL results to retain only significant features, then performs fine-mapping for each feature in parallel using the `mclapply` function.
#' @examples
#' \dontrun{
#' # Example with QTL results and genotype data
#' qtl_results <- data.table(phenotype_id = c("gene1", "gene2"), variant_id = c("rs1", "rs2"), slope = c(0.1, 0.2), slope_se = c(0.01, 0.02), qval = c(0.01, 0.04))
#' genotypes <- matrix(rnorm(200), nrow=20)
#' finemap_features_parallel(qtl_results, genotypes)
#' }
#' @import data.table
#' @import parallel
#' @import susieR
finemap_features_parallel <- function(qtl_results, genotypes, feature_column='phenotype_id', significance_columns=list('qval' = 0.05), variant_column='variant_id', slope_column='slope', se_column='slope_se', L=10, estimate_residual_variance=T, nthreads=detectCores()) {
  # subset the data to only contain significant featurs
  for (sig_feature in names(significance_columns)) {
    qtl_results <- qtl_results[qtl_results[[sig_feature]] < significance_columns[[sig_feature]], ]
  }
  # check which features are left
  features_present <- unique(qtl_results[[feature_column]])
  # do parallel processing of all the features
  finemapped_per_feature <- 
    mclapply(
      features_present, 
      FUN = function(feature) {
        # subset to this feature
        qtl_results_feature <- qtl_results[qtl_results[[feature_column]] == feature, ]
        # get finemapping result
        finemapped_feature <- finemap_feature(qtl_results_feature, 
                                              genotypes, 
                                              variant_column = variant_column, 
                                              slope_column = slope_column, 
                                              se_column = se_column, 
                                              L = L,
                                              estimate_residual_variance = estimate_residual_variance)
        
        return(list('feature' = feature, 'result' = finemapped_feature))
      }, 
      mc.cores = nthreads)
  return(finemapped_per_feature)
}


#' Finemap Features Serially
#'
#' This function performs fine-mapping for multiple features serially using summary statistics and genotype data.
#'
#' @param qtl_results A data.frame or data.table containing the QTL results.
#' @param genotypes A matrix, data.table, data.frame, or list containing genotype data.
#' @param feature_column A character string specifying the column name for feature IDs in `qtl_results`. Default is 'phenotype_id'.
#' @param significance_columns A named list specifying the significance thresholds for filtering features. Default is list('qval' = 0.05).
#' @param variant_column A character string specifying the column name for variant IDs in `qtl_results`. Default is 'variant_id'.
#' @param slope_column A character string specifying the column name for effect sizes in `qtl_results`. Default is 'slope'.
#' @param se_column A character string specifying the column name for standard errors in `qtl_results`. Default is 'slope_se'.
#' @param L An integer specifying the maximum number of causal variants to consider. Default is 10.
#' @param estimate_residual_variance A logical value indicating whether to estimate residual variance. Default is TRUE.
#' @param nthreads An integer specifying the number of threads to use for parallel processing. Default is the number of cores detected by `detectCores()`.
#' @return A list of fine-mapped results for each feature.
#' @details The function filters the QTL results to retain only significant features, then performs fine-mapping for each feature serially. It prints the feature being processed and updates the user every 100 features.
#' @examples
#' \dontrun{
#' # Example with QTL results and genotype data
#' qtl_results <- data.table(phenotype_id = c("gene1", "gene2"), variant_id = c("rs1", "rs2"), slope = c(0.1, 0.2), slope_se = c(0.01, 0.02), qval = c(0.01, 0.04))
#' genotypes <- matrix(rnorm(200), nrow=20)
#' finemap_features_serial(qtl_results, genotypes)
#' }
#' @import data.table
#' @import susieR
finemap_features_serial <- function(qtl_results, genotypes, feature_column='phenotype_id', significance_columns=list('qval' = 0.05), variant_column='variant_id', slope_column='slope', se_column='slope_se', L=10, estimate_residual_variance=T, nthreads=detectCores(), n_sample=NULL, n_sample_column=NULL) {
  # subset the data to only contain significant featurs
  for (sig_feature in names(significance_columns)) {
    qtl_results <- qtl_results[qtl_results[[sig_feature]] < significance_columns[[sig_feature]], ]
  }
  # check which features are left
  features_present <- unique(qtl_results[[feature_column]])
  # keep a counter
  i <- 0
  # do parallel processing of all the features
  finemapped_per_feature <- list()
  for (feature in features_present) {
    print(feature)
    # subset to this feature
    qtl_results_feature <- qtl_results[qtl_results[[feature_column]] == feature, ]
    # get finemapping result
    finemapped_feature <- finemap_feature(qtl_results_feature, 
                                          genotypes, 
                                          variant_column = variant_column, 
                                          slope_column = slope_column, 
                                          se_column = se_column, 
                                          L = L,
                                          estimate_residual_variance = estimate_residual_variance, 
                                          n_sample_column = n_sample_column, 
                                          n_sample = n_sample)
    finemapped_per_feature[[feature]] <- list('feature' = feature, 'result' = finemapped_feature)
    # update counter
    i <- i + 1
    if (i %% 100 == 0) {
      message(paste('processed', as.character(i), 'features'))
    }
  }
  return(finemapped_per_feature)
}


#' Debugging Function for QTL Fine-Mapping
#'
#' This function performs fine-mapping on QTL results using genotype data and saves the results to an RDS file.
#'
#' @return None. The function saves the fine-mapped results to an RDS file.
#' @details The function reads QTL results and genotype data, filters the QTL results based on significance, subsets the genotype data, and performs fine-mapping using the `finemap_features_serial` function. The results are then saved to an RDS file.
#' @examples
#' \dontrun{
#' do_debug()
#' }
#' @import data.table
#' @import snpStats
#' @import ReigenMT
#' @import susieR
do_debug <- function() {
  # initialize variables
  qtl_file <- '/groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/eqtl/tensorqtl/output/elmentaite_adult_martin_immune/cell_type_medhigh_inflammationsplit_mincor/unconfined/pcs/CD4_T_cells_merged.tsv.gz'
  genotype_file <- '/groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/genotype/combined/imputed/all/sc-eqtlgen-pipeline-wg1/lpmcv2_imputed_hg38_info_filled_rsid'
  output_rds <- '/groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/eqtl/finemapping/output/elmentaite_adult_martin_immune/cell_type_medhigh_inflammationsplit_mincor/unconfined/pcs/CD4_T_cells_merged.rds'
  # these always have a value
  significance_cutoff <- 0.05
  significance_column <- 'qval'
  num_threads <- 4
  variant_column <- 'variant_id'
  feature_column <- 'phenotype_id'
  slope_column <- 'slope'
  se_column <- 'slope_se'
  # this is easier flipped
  in_sample_id <- !(F)
  # read the QTL data
  qtl_result <- read.table(qtl_file, header = T, sep = '\t')
  # read genotype data
  genotype <- snpStats::read.plink(
    bed = paste(genotype_file, '.bed', sep = ''),
    bim = paste(genotype_file, '.bim', sep = ''),
    fam = paste(genotype_file, '.fam', sep = '')
  )
  # only look at significant results
  qtl_result <- qtl_result[qtl_result[[significance_column]] < significance_cutoff, ]
  # only use genotype data that we used in the QTLs
  genotype <- ReigenMT::subset_genotypes(genotype, unique(qtl_result[[variant_column]]))
  # create significance filtering list
  significance_columns <- list()
  significance_columns[[significance_column]] <- significance_cutoff
  # do finemapping
  finemapped <- finemap_features_serial(qtl_result, 
                                        genotypes=genotype, 
                                        feature_column=feature_column, 
                                        significance_columns=significance_columns, 
                                        variant_column=variant_column, 
                                        slope_column=slope_column, 
                                        se_column=se_column, 
                                        L=10, 
                                        estimate_residual_variance=in_sample_id, 
                                        nthreads=num_threads)
  # finally save the result
  saveRDS(finemapped, output_rds)
}


####################
# Settings        #
####################

set.seed(7777)


####################
# Debug            #
####################

#do_debug()


####################
# Main Code        #
####################


# make command line options
option_list <- list(
  make_option(c("-q", "--qtl_file"), type="character", default=NULL, 
              help="eqtl file to read", metavar="character"),
  make_option(c("-g", "--genotype_file"), type="character", default=NULL, 
              help="output file to write", metavar="character"),
  make_option(c("-s", "--significance_column"), type="character", default='qval', 
              help="column in the QTL data that has the significance [default]", metavar="character"), 
  make_option(c("-c", "--significance_cutoff"), type="numeric", default=0.05, 
              help="significance column to filter the QTL output on [default]", metavar="numeric"),
  make_option(c("-o", "--output_rds"), type="character", default=NULL, 
              help="rds output location of the finemapping", metavar="character"),
  make_option(c("-t", "--num_threads"), type="numeric", default=2, 
              help="number of parallel threads to use [default]", metavar="numeric"), 
  make_option(c("-v", "--variant_column"), type="character", default='variant_id', 
              help="column in the QTL data that has the variant identifier [default]", metavar="character"), 
  make_option(c("-f", "--feature_column"), type="character", default='phenotype_id', 
              help="column in the QTL data that has the feature identifier [default]", metavar="character"), 
  make_option(c("-b", "--slope_column"), type="character", default='slope', 
              help="column in the QTL data that has the slop [default]", metavar="character"), 
  make_option(c("-e", "--se_column"), type="character", default='slope_se', 
              help="column in the QTL data that has the standard error of the slope [default]", metavar="character"), 
  make_option(c("-n", "--not_in_sample_ld"), action="store_false", default=T,
              help="genotypes are not the same as used in QTL mapping [default]"), 
  make_option(c("-m", "--n_sample"), type="numeric", default=NULL, 
              help="number of samples", metavar="numeric"),
  make_option(c("-a", "--n_sample_column"), type="character", default=NULL, 
              help="the column in the QTL file that has the number of samples for that test", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
qtl_file <- NULL
genotype_file <- NULL
output_rds <- NULL
n_sample <- NULL
n_sample_column <- NULL
# these always have a value
significance_cutoff <- opt[['significance_cutoff']]
significance_column <- opt[['significance_column']]
num_threads <- opt[['num_threads']]
variant_column <- opt[['variant_column']]
feature_column <- opt[['feature_column']]
slope_column <- opt[['slope_column']]
se_column <- opt[['se_column']]
# this is easier flipped
in_sample_id <- !(opt[['not_in_sample_ld']])
# check variables
if (is.null(opt[['qtl_file']])) {
  stop(paste('-q/--qtl_file is an obligatory parameter\n'))
} else {
  qtl_file <- opt[['qtl_file']]
}
if (is.null(opt[['genotype_file']])) {
  stop(paste('-g/--genotype_file is an obligatory parameter\n'))
} else {
  genotype_file <- opt[['genotype_file']]
}
if (is.null(opt[['output_rds']])) {
  stop(paste('-o/--output_rds is an obligatory parameter\n'))
} else {
  output_rds <- opt[['output_rds']]
}

# these are optional
if (!is.null(opt[['n_sample']])) {
  n_sample <- opt[['n_sample']]
}
if (!is.null(opt[['n_sample_column']])) {
  n_sample_column <- opt[['n_sample_column']]
}

# check if the directory for the output exists, otherwise we would fail at the very last step
output_rds_dir <- dirname(output_rds)
if (!(dir.exists(output_rds_dir))) {
  stop(paste('directory for output RDS', output_rds, 'does not exist. This would mean the final step writing to file would fail, please create the directory first\n'))
}

# read the QTL data
qtl_result <- read.table(qtl_file, header = T, sep = '\t')

# check columns
if (!(variant_column %in% colnames(qtl_result))) {
  stop(paste('variant column', variant_column, 'not in column names of qtl data\n'))
}
if (!(feature_column %in% colnames(qtl_result))) {
  stop(paste('feature column', feature_column, 'not in column names of qtl data\n'))
}
if (!(slope_column %in% colnames(qtl_result))) {
  stop(paste('slope column', slope_column, 'not in column names of qtl data\n'))
}
if (!(se_column %in% colnames(qtl_result))) {
  stop(paste('standard error column', se_column, 'not in column names of qtl data\n'))
}

# read genotype data
genotype <- snpStats::read.plink(
  bed = paste(genotype_file, '.bed', sep = ''),
  bim = paste(genotype_file, '.bim', sep = ''),
  fam = paste(genotype_file, '.fam', sep = '')
)
# only look at significant results
qtl_result <- qtl_result[qtl_result[[significance_column]] < significance_cutoff, ]
# only use genotype data that we used in the QTLs
genotype <- ReigenMT::subset_genotypes(genotype, unique(qtl_result[[variant_column]]))
# create significance filtering list
significance_columns <- list()
significance_columns[[significance_column]] <- significance_cutoff
# do finemapping
finemapped <- finemap_features_serial(qtl_result, 
                                      genotypes=genotype, 
                                      feature_column=feature_column, 
                                      significance_columns=significance_columns, 
                                      variant_column=variant_column, 
                                      slope_column=slope_column, 
                                      se_column=se_column, 
                                      L=10, 
                                      estimate_residual_variance=in_sample_id, 
                                      nthreads=num_threads, 
                                      n_sample_column = n_sample_column, 
                                      n_sample = n_sample)
# finally save the result
saveRDS(finemapped, output_rds)