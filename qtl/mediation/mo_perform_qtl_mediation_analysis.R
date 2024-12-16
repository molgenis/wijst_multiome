#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_perform_qtl_mediation_analysis.R
# Function: perform mediation analysis of eQTLs through caQTLs
############################################################################################################################

####################
# libraries        #
####################

library(mediate)
library(data.table)


####################
# Functions        #
####################


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
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output location of analysis [default= %default]", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
