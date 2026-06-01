#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_calc_lds.R
# Function: calulate pair-wise LD between SNPs in single variant list or LD between variants in list 1 and list 2. NOTE! NA ARE CODED AS HIGH NUMBER! 1+308
# Example: 
# Rscript mo_calc_lds.R \
#   --genotypes /groups/umcg-franke-scrna/tmp04/external_datasets/sc-eqtlgen-imputation-ref-hg38/ref_panel_QC/30x-GRCh38-norsid \
#   --output_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/moldpairs/1000G_HC/allpop/immune-gwas-catalog-download-associations-alt-full-moldpairs \
#   --variant_list_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz \
#   --second_variant_list_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full-chromposrefalt.tsv.gz \ 
#   --variant_list_column snp_id \
#   --second_variant_list_column chromposaltref \
#   --ld_cutoff 0.1
# Example 2: 
# Rscript mo_calc_lds.R \
#   --genotypes /groups/umcg-franke-scrna/tmp04/external_datasets/sc-eqtlgen-imputation-ref-hg38/ref_panel_QC/30x-GRCh38-EUR-norsid \
#   --output_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_ld/1000G_HC/eurpop/moresparse/mo_qtl_variants_tested_cpeaks_overlap_chr \
#   --variant_list_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz \
#   --variant_list_column snp_id \
#   --ld_cutoff 0.1
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
  make_option(c("-l", "--ld_cutoff"), type="numeric", default=0.0, 
              help="make values smaller than the cutoff into 0, so the matrix is more sparse", metavar="numeric"), 
  make_option(c("-n", "--chromosomes"), type="character", default=NULL, 
              help="comma separacted string of specific chromosomes to calculate LD for", metavar="character"), 
  make_option(c("-k", "--n_chunks"), type="numeric", default=NULL, 
              help="chunks to divide work into for l1 vs l2 type of LDs", metavar="numeric")
)


# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# initialize the options in debug mode
if (debug) {
  opt <- list(
    genotypes = '/groups/umcg-franke-scrna/tmp04/external_datasets/sc-eqtlgen-imputation-ref-hg38/ref_panel_QC/30x-GRCh38-EUR-norsid',
    output_file = '/groups/umcg-fg/tmp04/projects/mpra/generate_panel/mo_qtl_1000g_ld/ld_matrices/mpra_top_qtl_1000g_eur_chr',
    variant_list_file = '/groups/umcg-fg/tmp04/projects/mpra/generate_panel/mo_qtl_1000g_ld/variants_mo/mo_qtl_top_variants_chr1.txt.gz',
    # variant_list_column = 'snp_id',
    second_variant_list_file = '/groups/umcg-fg/tmp04/projects/mpra/generate_panel/mo_qtl_1000g_ld/variants_1000g/1000g_all_eur_variants_chr1.txt.gz', 
    # second_variant_list_column = 'chromposaltref',
    depth = 1000, 
    ld_cutoff = 0.29, 
    chromosomes = '1', 
    n_chunks = 3
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
# set chromosomes string
chromosomes_string <- opt[['chromosomes']]
# as vector
chromosomes_to_do <- NULL
# then into list
if (!is.null(chromosomes_string)) {
  chromosomes_to_do <- strsplit(chromosomes_string, ',')[[1]]
}
# number of chunks
n_chunks <- opt[['n_chunks']]

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
# if we defined the chromosome, filter
if (!is.null(chromosomes_to_do)) {
  chromosomes_gt <- intersect(chromosomes_gt, chromosomes_to_do)
}

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
    message(paste('calculating LD between', ncol(genotypes_chr$genotypes), 'variants in list 1'))
    # let the user know chunking doesn't work
    if (!is.null(n_chunks)) {
      warning('chunking is not supported (or required) for full pairwise LD calculation')
    }
    # calculate all the LD, because we only have the variants in the first list anyway
    ld_mat_chr <- snpStats::ld(genotypes_chr$genotypes, depth = depth, stats = "R.squared", symmetric = T)
    # make more sparse
    if (ld_cutoff > 0) {
      message(paste('making values smaller than', ld_cutoff, 'into 0'))
      ld_mat_chr@x[ld_mat_chr@x < ld_cutoff] <- 0
      ld_mat_chr <- drop0(ld_mat_chr)
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
    # set the LD low for where the value does not make sense, remove R2 much larger than 1, as these are the NA values
    ld_mat_chr@x[ld_mat_chr@x > 1.001] <- -0.0001
    # and explicit NAs as well
    ld_mat_chr@x[is.na(ld_mat_chr@x)] <- -0.0001
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
    message(paste('calculating LD between', ncol(genotypes_chr_l1$genotypes), 'variants in list 1 and', ncol(genotypes_chr_l2$genotypes), 'variants in list 2'))
    # calculate the LD between the two lists
    if (is.null(n_chunks)) {
      ld_mat_chr <- snpStats::ld(genotypes_chr_l1$genotypes, genotypes_chr_l2$genotypes, stats = "R.squared")
      # now make the matrix sparse
      ld_mat_chr <- Matrix(ld_mat_chr, sparse = TRUE)
      # make more sparse
      if (ld_cutoff > 0) {
        message(paste('making values smaller than', ld_cutoff, 'into 0'))
        # ld_mat_chr[ld_mat_chr < ld_cutoff] <- 0
        ld_mat_chr@x[ld_mat_chr@x < ld_cutoff] <- 0
        ld_mat_chr <- drop0(ld_mat_chr)
      }
    } else {
      # in chunks if required
      ld_per_chunk <- list()
      # get variants in l2
      l2_variants <- colnames(genotypes_chr$genotypes)
      # get the size per chunk
      size_per_chunk <- floor(length(l2_variants) / n_chunks)
      # start index
      l2_variant_i <- 1
      # and check in chunks
      while(l2_variant_i < length(l2_variants)) {
        # get the end index
        end_chunk_i <- l2_variant_i + size_per_chunk - 1
        # if we go over the max size, just go to max size
        if (end_chunk_i > length(l2_variants)) {
          end_chunk_i <- length(l2_variants)
        }
        # let the user where we are
        message(paste('doing chunk from', as.character(l2_variant_i), 'to', as.character(end_chunk_i), 'out of total', length(l2_variants), 'variants'))
        # extract the variants
        l2_variants_chunk <- l2_variants[l2_variant_i : end_chunk_i]
        # subset set number two
        genotypes_chr_l2_chunk <- list(
          genotypes = genotypes_chr_l2$genotypes[, colnames(genotypes_chr_l2$genotypes) %in% l2_variants_chunk],
          map = genotypes_chr_l2$map[rownames(genotypes_chr_l2$map) %in% l2_variants_chunk, ],
          fam = genotypes_chr_l2$fam
        )
        # calculate LD
        ld_mat_chr_chunk <- snpStats::ld(genotypes_chr_l1$genotypes, genotypes_chr_l2_chunk$genotypes, stats = "R.squared")
        # then make sparse
        ld_mat_chr_chunk <- Matrix(ld_mat_chr_chunk, sparse = TRUE)
        # make more sparse
        if (ld_cutoff > 0) {
          message(paste('making values smaller than', ld_cutoff, 'into 0'))
          ld_mat_chr_chunk@x[ld_mat_chr_chunk@x < ld_cutoff] <- 0
          ld_mat_chr_chunk <- drop0(ld_mat_chr_chunk)
        }
        # put in the list
        ld_per_chunk[[paste(l2_variant_i, end_chunk_i, sep = '-')]] <- ld_mat_chr_chunk
        # increase index
        l2_variant_i <- l2_variant_i + size_per_chunk
      }
      # merge all by column
      ld_mat_chr <- do.call('cbind', ld_per_chunk)
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

    # set the LD low for where the value does not make sense, remove R2 much larger than 1, as these are the NA values
    # ld_mat_chr[ld_mat_chr > 1.001] <- -0.0001
    ld_mat_chr@x[ld_mat_chr@x > 1.001] <- -0.0001
    # and explicit NAs as well
    ld_mat_chr@x[is.na(ld_mat_chr@x)] <- -0.0001

    # set the output file loc
    output_file_chrom_full <- paste(output_file, chromosome, '.mtx', sep = '')
    Matrix::writeMM(ld_mat_chr, output_file_chrom_full)
    # zip this file
    gzip(output_file_chrom_full)
    # then make a checksum
    mdfiver::create_sha256_for_file(paste0(output_file_chrom_full, '.gz'))
  }
}
