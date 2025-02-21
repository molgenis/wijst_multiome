#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Marc-Jan Bonder, Roy Oelen
# Name: mo_meta_analyse_creqtl_cres.R
# Function: use CRE-QTL input files to meta-analyse normal CREs
# example: 
# ~/start_Rscript.sh /groups/umcg-franke-scrna/tmp02/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_meta_analyse_creqtl_cres.R \
# --beta_file /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/betas_ps/monocyte/betas.tsv.gz \
# --se_file /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/betas_ps/monocyte/ses.tsv.gz \
# --output_file /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/betas_ps/monocyte/meta_result.tsv.gz \
# --cell_proportions /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/matrices/monocyte/ncells.tsv.gz \
# --smf_file /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/matrices/monocyte/sample_mapping.tsv.gz \
# --donor_column sample_final \
# --meta_sample_dups

############################################################################################################################


####################
# libraries        #
####################

library(qvalue)
library(optparse)
library(mdfiver)


####################
# Functions        #
####################

#' Convert Z-Score to P-Value
#'
#' This function converts a Z-score to a p-value. It can handle both one-sided and two-sided tests.
#'
#' @param z Numeric. The Z-score to be converted.
#' @param one.sided Character. Specifies the type of test. Use `NULL` for a two-sided test, `"-"` for a one-sided test where the alternative hypothesis is that the true mean is less than the hypothesized mean, and any other value for a one-sided test where the alternative hypothesis is that the true mean is greater than the hypothesized mean.
#'
#' @return Numeric. The p-value corresponding to the given Z-score.
#'
#' @examples
#' convert_z_score(1.96) # Two-sided test
#' convert_z_score(1.96, one.sided="-") # One-sided test (less than)
#' convert_z_score(1.96, one.sided="+") # One-sided test (greater than)
convert_z_score <- function(z, one.sided=NULL) {
  if(is.null(one.sided)) {
    pval = pnorm(-abs(z));
    pval = 2 * pval
  } else if(one.sided=="-") {
    pval = pnorm(z);
  } else {
    pval = pnorm(-z);
  }
  return(pval);
}


#' Meta-Analyse Effect Sizes
#'
#' This function performs a meta-analysis on effect sizes and their standard errors, returning the combined effect size, standard error, and total weight.
#'
#' @param betas Numeric vector. The effect sizes to be meta-analysed.
#' @param beta_ses Numeric vector. The standard errors of the effect sizes.
#' @param weights Numeric vector. The weights for each effect size.
#'
#' @return A list containing the meta-analysed effect size (`beta`), the meta-analysed standard error (`se`), and the total weight (`weight`).
#'
#' @examples
#' betas <- c(0.2, 0.3, 0.4)
#' beta_ses <- c(0.1, 0.1, 0.1)
#' weights <- c(1, 1, 1)
#' meta_analyse(betas, beta_ses, weights)
meta_analyse <- function(betas, beta_ses, weights) {
  # calculate variance from standard errors
  variance = beta_ses^2
  # get beta divided by variance
  effect_size_divided_by_variance = betas / variance
  # do for each row, so each region-gene pair
  effect_size_divided_by_variance_total = rowSums(effect_size_divided_by_variance,na.rm=T)
  # do the same for one
  one_divided_by_variance = 1 / variance
  one_divided_by_variance_total = rowSums(one_divided_by_variance, na.rm = T)
  # use this to meta beta and meta se
  meta_analysed_effect_size = effect_size_divided_by_variance_total / one_divided_by_variance_total
  meta_analysed_standard_error = sqrt(1 / one_divided_by_variance_total)
  # get the sum of weights
  meta_weight <- sum(weights)
  # put in a list
  meta_stats <- list('beta' = meta_analysed_effect_size, 'se' = meta_analysed_standard_error, 'weight' = meta_weight)
  return(meta_stats)
}


####################
# Settings        #
####################


####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-b", "--beta_file"), type="character", default=NULL, 
              help="file of betas to read", metavar="character"),
  make_option(c("-s", "--se_file"), type="character", default=NULL, 
              help="file of standard errors to read", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="output file to save meta-analysis results", metavar="character"), 
  make_option(c("-p", "--cell_proportions"), type="character", default=NULL, 
              help="tab separated file, linking number of cells to samples", metavar="character"),
  make_option(c("-m", "--smf_file"), type="character", default=NULL, 
              help="sample mapping file, linking sample to specific donor [default]", metavar="character"), 
  make_option(c("-a", "--sample_column"), type="character", default='sample', 
              help="column in the sample mapping file that has the sample name [default]", metavar="character"), 
  make_option(c("-d", "--donor_column"), type="character", default='donor', 
              help="column in the sample mapping file that has the donor [default]", metavar="character"), 
  
  make_option(c("-f", "--meta_sample_dups"), action="store_true", default=F,
              help="meta-analyse the samples from the same donor first, instead of meta-analysing all in one go [default]"), 
  make_option(c("-u", "--dedeup_sample_dups"), action="store_true", default=F,
              help="remove second sample of same donor, keeping one with most cells, instead of meta-analysing all in one go [default]")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize variables
beta_loc <- NULL
se_loc <- NULL
ncell_loc <- NULL
smf_loc <- NULL
out_loc <- NULL
# and read the values we definitely have
meta_analyse_sample_duplicates <- opt[['meta_sample_dups']]
deduplicate_sample_duplicates <- opt[['dedeup_sample_dups']]
sample_column <- opt[['sample_column']]
donor_column <- opt[['donor_column']]

# we need to check if these parameters were given, otherwise we'll immediatly error out
if (is.null(opt[['beta_file']])) {
  stop(paste('-b/--beta_file is an obligatory parameter\n'))
} else {
  beta_loc <- opt[['beta_file']]
}
if (is.null(opt[['se_file']])) {
  stop(paste('-s/--se_file is an obligatory parameter\n'))
} else {
  se_loc <- opt[['se_file']]
}
if (is.null(opt[['output_file']])) {
  stop(paste('-o/--output_file is an obligatory parameter\n'))
} else {
  out_loc <- opt[['output_file']]
}
if (is.null(opt[['cell_proportions']])) {
  stop(paste('-p/--cell_proportions is an obligatory parameter\n'))
} else {
  ncell_loc <- opt[['cell_proportions']]
}
if (is.null(opt[['smf_file']])) {
  stop(paste('-m/--smf_file is an obligatory parameter\n'))
} else {
  smf_loc <- opt[['smf_file']]
}
if (meta_analyse_sample_duplicates & deduplicate_sample_duplicates) {
  stop('paste, only -f/--meta_sample_dups or -u/-dedeup_sample_dups can be set. They cannot be set at the same time\n')
}

# read the input files
beta <- read.delim(beta_loc)
beta_se <- read.delim(se_loc)
smf <- read.delim(smf_loc)
nCells <- read.delim(ncell_loc)

# keep only the things where we have all data
stat_columns <- c('region', 'gene')
beta_samples <- setdiff(colnames(beta), stat_columns)
beta_se_samples <- setdiff(colnames(beta_se), stat_columns)
smf_samples <- smf[[sample_column]]
ncell_samples <- nCells[['sample']]
# check which we have
all_samples <- intersect(beta_samples, beta_se_samples)
all_samples <- intersect(all_samples, smf_samples)
all_samples <- intersect(all_samples, ncell_samples)
# and use the same order for all of them
beta <- beta[, c(stat_columns, all_samples)]
beta_se <- beta_se[, c(stat_columns, all_samples)]
smf <- smf[match(all_samples, smf[[sample_column]]), ]
nCells <- nCells[match(all_samples, nCells[['sample']]), ]

# keep track of the duplicate sample problem
dupsample_solve <- 'none'

# check if we need to keep only the unique donors
if (deduplicate_sample_duplicates) {
  # get the order of the cell numbers
  cell_nr_order <- order(nCells[['ncell']], decreasing = T)
  # now order all the tables
  nCells <- nCells[cell_nr_order, ]
  beta <- beta[, c(1, 2 ,(cell_nr_order+2))] # 1 and 2 are the region and gene ones
  beta_se <- beta_se[, c(1, 2 ,(cell_nr_order+2))] # 1 and 2 are the region and gene ones
  smf <- smf[cell_nr_order, ]
  # now deduplicate
  indices_deduplicated <- !duplicated(smf[[donor_column]])
  # and use the deduplicated indices
  smf <- smf[indices_deduplicated, ]
  beta <- beta[, c(T, T ,(indices_deduplicated))] # 1 and 2 are the region and gene ones
  beta_se <- beta_se[, c(T, T ,(indices_deduplicated))] # 1 and 2 are the region and gene ones
  nCells <- nCells[indices_deduplicated, ]
  # update the duplicate sample solve
  dupsample_solve <- 'keep_biggest'
}


# get distribution of cells
ncell_dist <- paste(as.character(min(nCells[['ncell']])),
                    as.character(quantile(nCells[['ncell']])[['25%']]),
                    as.character(quantile(nCells[['ncell']])[['50%']]),
                    as.character(quantile(nCells[['ncell']])[['75%']]),
                    as.character(max(nCells[['ncell']])),
                    sep = ';'
)
# and total
ncell_tot <- sum(nCells[['ncell']])
# and number of samples
n_sample <- nrow(smf)
# and number of donors
n_donor <- length(unique(smf[[donor_column]]))

# check if we need to meta-analyse duplicate samples first
if (meta_analyse_sample_duplicates) {
  # check each unique donor
  unique_donor <- unique(smf[[donor_column]])
  # make deduplicated matrices
  dedup_betas <- data.frame(matrix(NA, nrow = nrow(beta), ncol = length(unique_donor), dimnames = list(NULL, unique_donor)))
  dedup_beta_ses <- data.frame(matrix(NA, nrow = nrow(beta_se), ncol = length(unique_donor), dimnames = list(NULL, unique_donor)))
  dedup_nCells <- data.frame(matrix(NA, nrow = length(unique_donor), ncol = 2, dimnames = list(unique_donor, c('sample', 'ncell'))))
  # check each of these donors
  for (donor in unique_donor) {
    # get the samples for that donor
    donor_samples <- smf[smf[[donor_column]] == donor, ][[sample_column]]
    # check if there is more than one sample
    if (length(donor_samples) > 1) {
      # now subset the three table for those samples
      donor_betas <- beta[, match(donor_samples, colnames(beta)), drop = F]
      donor_ses <- beta_se[, match(donor_samples, colnames(beta_se)), drop = F]
      donor_nCells <- nCells[match(donor_samples, nCells[['sample']]), , drop = F]
      # get meta stats
      donor_stats <- meta_analyse(donor_betas[ , setdiff(colnames(donor_betas), stat_columns)], donor_ses[ , setdiff(colnames(donor_ses), stat_columns)], donor_nCells[['ncell']])
      # put them into the tables
      dedup_betas[[donor]] <- donor_stats[['beta']]
      dedup_beta_ses[[donor]] <- donor_stats[['se']]
      dedup_nCells[donor, 'sample'] <- donor
      dedup_nCells[donor, 'ncell'] <- donor_stats[['weight']]
    }
    # otherwise it is just the stats of the one sample
    else {
      dedup_betas[[donor]] <- beta[[donor_samples[1]]]
      dedup_beta_ses[[donor]] <- beta_se[[donor_samples[1]]]
      dedup_nCells[donor, 'sample'] <- donor
      dedup_nCells[donor, 'ncell'] <- nCells[nCells[['sample']] == donor_samples[1], 'ncell']
    }
  }
  # the smf is simplified to just be donor-to-donor
  dedup_smf <- data.frame(x = unique_donor, y = unique_donor)
  colnames(dedup_smf) <- c(sample_column, donor_column)
  # make these the new matrices
  beta <- cbind(beta[, stat_columns], dedup_betas)
  beta_se <- cbind(beta_se[, stat_columns], dedup_beta_ses)
  nCells <- dedup_nCells
  smf <- dedup_smf
  # and clear memory of old variables
  rm(dedup_betas)
  rm(dedup_beta_ses)
  rm(dedup_nCells)
  rm(dedup_smf)
  # update the duplicate sample solve
  dupsample_solve <- 'meta_sep'
}

# grab the descriptive columns
beta_stats <- beta[, stat_columns]
# now do the meta-analysis with only the numeric columns
meta_stats <- meta_analyse(beta[ , setdiff(colnames(beta), stat_columns)], beta_se[ , setdiff(colnames(beta_se), stat_columns)], nCells[['ncell']])
# get the Z score
meta_z_score = meta_stats[['beta']] / meta_stats[['se']]
# make Z score into p value
meta_p_value = convert_z_score(meta_z_score)
# calculate q value from p
meta_q_value = qvalue::qvalue(meta_p_value,pi0 = 1)$qvalue

# put all the results in a table
res_tbl <- data.frame(
  'meta_beta' = meta_stats[['beta']], 
  'meta_se' = meta_stats[['se']], 
  'meta_z' = meta_z_score, 
  'meta_p' = meta_p_value, 
  'meta_q' = meta_q_value, 
  'n_sample' = rep(n_sample, times = length(meta_z_score)), 
  'n_donor' = rep(n_donor, times = length(meta_z_score)), 
  'cell_total' = rep(ncell_tot, times = length(meta_z_score)), 
  'cell_dist' = rep(ncell_dist, times = length(meta_z_score)), 
  'dupsample_solve' = rep(dupsample_solve, times = length(meta_z_score))
)
# add decriptive columns back
res_tbl <- cbind(beta_stats, res_tbl)

# get zipped output loc
output_loc_gz <- out_loc
# gz file ends with .gz
if (grepl('.gz$', out_loc)) {
  # gzip if ends with .gz
  output_loc_gz <- gzfile(out_loc)
}
write.table(res_tbl, output_loc_gz, sep = '\t', row.names = F, col.names = T, quote = F)
# also make a checksum
mdfiver::create_md5_for_file(out_loc)
