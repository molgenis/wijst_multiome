#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_calc_lds.R
# Function: 
# Example: 
# Rscript mo_calc_lds.R \
#   --genotypes /groups/umcg-franke-scrna/tmp04/external_datasets/sc-eqtlgen-imputation-ref-hg38/ref_panel_QC/30x-GRCh38-norsid \
#   --output_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full-moldpairs \
#   --variant_list_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz \
#   --second_variant_list_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full-chromposrefalt.tsv.gz \ 
#   --variant_list_column snp_id \
#   --second_variant_list_column chromposaltref
#
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)
# use plink files
library(snpStats)
# use sparse matrices
library(Matrix)
# load command line parameters
library(optparse)
# for zipping the matrix
library(R.utils)


####################
# Functions        #
####################



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
  make_option(c("-g", "--genotypes"), type="character", default=NULL, 
              help="plink based path genotype file", metavar="character"),
  make_option(c("-o", "--output_file"), type="character", default=NULL, 
              help="base output file", metavar="character"), 
  make_option(c("-v", "--variant_list_file"), type="character", default=NULL, 
              help="variant list file (default without header)", metavar="character"),
  make_option(c("-a", "--second_variant_list_file"), type="character", default=NULL, 
              help="(optional) second variant file list to calculate LD against (default without header)", metavar="character"), 
  make_option(c("-c", "--variant_list_column"), type="character", default=NULL, 
              help="column in variant list file containing variant ID", metavar="character"),
  make_option(c("-s", "--second_variant_list_column"), type="character", default=NULL, 
              help="column in second variant list file containing variant ID", metavar="character"),
  make_option(c("-d", "--depth"), type="numeric", default=1000, 
              help="number of consecutive variants to check for LD", metavar="numeric"), 
  make_option(c("-l", "--ld_cutoff"), type="numeric", default=0.1, 
              help="make values smaller than the cutoff into 0, so the matrix is more sparse", metavar="numeric")
)


# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize the options in debug mode
if (debug) {
  opt <- list(
    genotypes = '/groups/umcg-franke-scrna/tmp04/external_datasets/sc-eqtlgen-imputation-ref-hg38/ref_panel_QC/30x-GRCh38-norsid',
    # output_file = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_ld/mo_qtl_variants_ld_chr',
    output_file = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full-moldpairs',
    variant_list_file = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz',
    variant_list_column = 'snp_id',
    second_variant_list_file = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full-chromposrefalt.tsv.gz', 
    second_variant_list_column = 'chromposaltref',
    depth = 1000, 
    ld_cutoff = .1
  )
}
# set non-optional parameters
genotypes_loc <- NULL
output_file <- NULL
variant_list_file <- NULL
# read them
if (!is.null(opt[['genotypes']])) {
  genotypes_loc <- opt[['genotypes']]
} else {
  print_help(opt_parser)
  stop("Please provide the path to the genotype files (without .bed/.bim/.fam)", call.=FALSE)
}
if(!is.null(opt[['output_file']])) {
  output_file <- opt[['output_file']]
} else {
  print_help(opt_parser)
  stop("Please provide the path to the output file", call.=FALSE)
}
if (!is.null(opt[['variant_list_file']])) {
  variant_list_file <- opt[['variant_list_file']]
} else {
  print_help(opt_parser)
  stop("Please provide the path to the variant list file", call.=FALSE)
}
# extract the optional parameters
depth <- as.numeric(opt[['depth']])
# get ld cutoff
ld_cutoff <- as.numeric(opt[['ld_cutoff']])

# load first set of variants
first_variants <- NULL
if (is.null(opt[['variant_list_column']])) {
  first_variants <- fread(variant_list_file, header = F)[[1]]
} else {
  first_variants <- fread(variant_list_file, header = T)[[opt[['variant_list_column']]]]
}
# make sure they are unique
first_variants <- unique(first_variants)

# load the second set of variants if provided
second_variants <- NULL
if (!is.null(opt[['second_variant_list_file']])) {
  if (is.null(opt[['second_variant_list_column']])) {
    second_variants <- fread(opt[['second_variant_list_file']], header = F)[[1]]
  } else {
    second_variants <- fread(opt[['second_variant_list_file']], header = T)[[opt[['second_variant_list_column']]]]
  }
  # make sure they are unique
  second_variants <- unique(second_variants)
}
# make the outer join of variants
variants <- unique(c(first_variants, second_variants))

# read the bim
variants_in_gt <- fread(paste(genotypes_loc, '.bim', sep = ''), header = F)[[2]]
# get overlapping variants
overlapping_variants <- intersect(variants, variants_in_gt)
# read the genotypes, but only those in the file and in the confinement
genotypes <- read.plink(
  bed = paste(genotypes_loc, '.bed', sep = ''),
  bim = paste(genotypes_loc, '.bim', sep = ''),
  fam = paste(genotypes_loc, '.fam', sep = ''), 
  select.snps = overlapping_variants
)
# subset the two sets of variants as well
first_variants <- first_variants[first_variants %in% overlapping_variants]
if (!is.null(second_variants)) {
  second_variants <- second_variants[second_variants %in% overlapping_variants]
}
overlapping_variants <- unique(c(first_variants, second_variants))
# get the chromosomes present
chromosomes_gt <- unique(genotypes$map$chromosome)

# check each chromosome
for (chromosome in chromosomes_gt) {
  message(paste('doing chromosome', chromosome))
  # get variants on this chromosome
  variants_chr <- rownames(genotypes$map[genotypes$map$chromosome == chromosome, ])
  # subset the snpstats object to these variants
  genotypes_chr <- list(
    genotypes = genotypes$genotypes[, variants_chr],
    map = genotypes$map[variants_chr, ],
    fam = genotypes$fam
  )
  # order the variants by position
  variants_order <- rownames(genotypes_chr$map[order(genotypes_chr$map$position), ])
  # remove duplicates
  variants_order <- gsub('\\.[0-9]+$', '', variants_order)
  variants_order <- unique(variants_order)
  # and order them
  genotypes_chr <- list(
    genotypes = genotypes_chr$genotypes[, variants_order],
    map = genotypes_chr$map[variants_order, ],
    fam = genotypes_chr$fam
  )
  # set attributes
  attr(genotypes_chr$genotypes, "snp.support") <- genotypes_chr$map
  # initialize
  ld_mat_chr <- NULL
  # depending on whether we have a second list, use a different approach
  if (is.null(second_variants)) {
    # calculate all the LD, because we only have the variants in the first list anyway
    ld_mat_chr <- snpStats::ld(genotypes_chr$genotypes, depth = depth, stats = "R.squared", symmetric = T)
    # make more sparse
    if (ld_cutoff > 0) {
      ld_mat_chr[abs(ld_mat_chr) < ld_cutoff] <- 0
    }
    # write variants as well
    ld_mat_chr_rows <- rownames(ld_mat_chr)
    ld_mat_chr_cols <- colnames(ld_mat_chr)
    # write these
    output_file_chrom_full_rows_loc <- paste(output_file, chromosome, '.y.txt.gz', sep = '')
    output_file_chrom_full_cols_loc <- paste(output_file, chromosome, '.x.txt.gz', sep = '')
    write.table(ld_mat_chr_rows, gzfile(output_file_chrom_full_rows_loc), row.names = F, col.names = F, quote = F)
    write.table(ld_mat_chr_cols, gzfile(output_file_chrom_full_cols_loc), row.names = F, col.names = F, quote = F)
    # make checksums
    mdfiver::create_sha256_for_file(output_file_chrom_full_rows_loc)
    mdfiver::create_sha256_for_file(output_file_chrom_full_cols_loc)
    # set the output file loc
    output_file_chrom_full <- paste(output_file, chromosome, '.mtx', sep = '')
    Matrix::writeMM(ld_mat_chr, output_file_chrom_full)
    # zip this file
    gzip(output_file_chrom_full)
    # then make a checksum
    mdfiver::create_sha256_for_file(paste0(output_file_chrom_full, '.gz'))
  } else {
    # now subset to list only the variants in list 1
    genotypes_chr_l1 <- list(
      genotypes = genotypes_chr$genotypes[, colnames(genotypes_chr$genotypes) %in% first_variants],
      map = genotypes_chr$map[rownames(genotypes_chr$map) %in% first_variants, ],
      fam = genotypes_chr$fam
    )
    # subset to list 2 as well
    genotypes_chr_l2 <- list(
      genotypes = genotypes_chr$genotypes[, colnames(genotypes_chr$genotypes) %in% second_variants],
      map = genotypes_chr$map[rownames(genotypes_chr$map) %in% second_variants, ],
      fam = genotypes_chr$fam
    )
    # calculate the LD between the two lists
    ld_mat_chr <- snpStats::ld(genotypes_chr_l1$genotypes, genotypes_chr_l2$genotypes, stats = "R.squared")
    # make more sparse
    if (ld_cutoff > 0) {
      ld_mat_chr[abs(ld_mat_chr) < ld_cutoff] <- 0
    }
    # write variants as well
    ld_mat_chr_rows <- rownames(ld_mat_chr)
    ld_mat_chr_cols <- colnames(ld_mat_chr)
    # write these
    output_file_chrom_full_rows_loc <- paste(output_file, chromosome, '.y.txt.gz', sep = '')
    output_file_chrom_full_cols_loc <- paste(output_file, chromosome, '.x.txt.gz', sep = '')
    write.table(ld_mat_chr_rows, gzfile(output_file_chrom_full_rows_loc), row.names = F, col.names = F, quote = F)
    write.table(ld_mat_chr_cols, gzfile(output_file_chrom_full_cols_loc), row.names = F, col.names = F, quote = F)
    # make checksums
    mdfiver::create_sha256_for_file(output_file_chrom_full_rows_loc)
    mdfiver::create_sha256_for_file(output_file_chrom_full_cols_loc)
    # now make the matrix sparse
    ld_mat_chr <- Matrix(ld_mat_chr, sparse = TRUE)
    # set the output file loc
    output_file_chrom_full <- paste(output_file, chromosome, '.mtx', sep = '')
    Matrix::writeMM(ld_mat_chr, output_file_chrom_full)
    # zip this file
    gzip(output_file_chrom_full)
    # then make a checksum
    mdfiver::create_sha256_for_file(paste0(output_file_chrom_full, '.gz'))
  }
}
