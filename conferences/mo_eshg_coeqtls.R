#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_eshg_coeqtls.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(foreach)
library(doParallel)
library(vcfR)
library(data.table)
library(ggplot2)


####################
# objects.         #
####################
# object to hold the genotype information
Genotypes <- setRefClass('Genotypes',
                         fields = list(
                           # the snps
                           calls = 'data.table',
                           # the metadata
                           meta = 'data.table'
                         ),
                         methods = list(
                           get_genotypes = function(.self, snp_id=NULL, chrom_pos=NULL, snp_pos=NULL, convert_alleles=F, to_numeric=F) {
                             # we need the index of the SNP
                             snp_index <- NULL
                             if (!is.null(snp_id)) {
                               # extract snp index
                               snp_index <- grep(snp_id, .self$meta[['ID']])
                               
                             } else if (!is.null(chrom_pos) & !is.null(snp_pos)) {
                               # extract snp index
                               snp_index <- which(.self$meta[['CHROM']] == chrom_pos & .self$meta[['POS']] == snp_pos)
                               
                             } else {
                               stop("Need to supply the variant ID, or both the CHROM and SNP positions")
                             }
                             # grab the genotypes
                             genotypes_raw <- as.vector(unlist(.self$calls[snp_index, ]))
                             # convert from phased
                             #genotypes_raw <- gsub('', '|', genotypes_raw)
                             # convert to alleles if requested
                             if (convert_alleles) {
                               # get the ref
                               ref <- .self$meta[snp_index, 'REF'][['REF']][1]
                               # and alt
                               alt <- .self$meta[snp_index, 'ALT'][['ALT']][1]
                               # harmonize placement of slash
                               genotypes_raw <- gsub('1\\|0', '0\\|1', genotypes_raw)
                               #genotypes_raw <- gsub('10', '01', genotypes_raw)
                               # replace letters
                               genotypes_raw <- gsub('0', ref, genotypes_raw)
                               genotypes_raw <- gsub('1', alt, genotypes_raw)
                             }
                             if (to_numeric) {
                               genotypes_raw <- as.numeric(as.factor(genotypes_raw))
                             }
                             # turn into list
                             genotypes_raw <- as.list(genotypes_raw)
                             # add the participants as the keys
                             names(genotypes_raw) <- colnames(.self$calls)
                             return(genotypes_raw)
                           },
                           get_available_variants = function(.self){
                             return(.self$meta[['ID']])
                           }
                         )
)


####################
# Functions        #
####################


get_correlation_subset <- function(seurat_object, gene_a, gene_b, subset_value, subset_column='assignment', method='spearman') {
  # subset to the relevant subset
  seurat_subset <- seurat_object[, !is.na(seurat_object@meta.data[[subset_column]]) & seurat_object@meta.data[[subset_column]] == subset_value]
  # get the correlation
  cor_result <- cor.test(
    x = as.vector(unlist(seurat_subset@assays$SCT@counts[gene_a, ])),
    y = as.vector(unlist(seurat_subset@assays$SCT@counts[gene_b, ])),
    method = method
  )
  cor_val <- cor_result$estimate[[1]]
  return(cor_val)
}

get_correlations_subsets <- function(seurat_object, gene_a, gene_b, subset_column='assignment', method='spearman') {
  # only works if the values are strings
  seurat_object@meta.data[[subset_column]] <- as.character(seurat_object@meta.data[[subset_column]])
  # get the unique values in the subset column
  subset_values <- unique(seurat_object@meta.data[[subset_column]])
  # remove any empty values there might be
  subset_values <- subset_values[!is.na(subset_values)]
  # loop through all values
  result_per_subset <- foreach(i = 1:length(subset_values)) %do% {
    cor_val <- get_correlation_subset(seurat_object, gene_a, gene_b, subset_values[i], subset_column, method)
    # put in small dataframe
    res <- data.frame(val = c(subset_values[i]), cor = cor_val)
    return(res)
  }
  # merge all
  cors <- do.call('rbind', result_per_subset)
  return(cors)
}


get_correlations_per_condition <- function(seurat_object, gene_a, gene_b, condition_column='timepoint', subset_column='assignment', method='spearman') {
  # make condition string
  seurat_object@meta.data[[condition_column]] <- as.character(seurat_object@meta.data[[condition_column]])
  # do each condition
  conditions_to_do <- unique(seurat_object@meta.data[[condition_column]])
  # remove any empty
  conditions_to_do <- conditions_to_do[!is.na(conditions_to_do)]
  # do each condition
  conditions_list <- list()
  for (condition in conditions_to_do) {
    res_condition <- get_correlations_subsets(seurat_object = seurat_object[, !is.na(seurat_object@meta.data[[condition_column]]) & seurat_object@meta.data[[condition_column]] == condition], 
                                              gene_a = gene_a, 
                                              gene_b = gene_b, 
                                              subset_column = subset_column, 
                                              method = method)
    # add condition
    res_condition[['condition']] <- condition
    # put in list
    conditions_list[[condition]] <- res_condition
  }
  # merge conditions
  conditions_all <- do.call('rbind', conditions_list)
  return(conditions_all)
}


#' turn vcfR object into Genotype object
#' 
#' @param genotypes_vcf the vcfR input file to convert to Genotype object
#' @returns Genotype object
create_genotypes_object <- function(genotypes_vcf){
  # get the donors in the genotype data
  donors_vcf <- colnames(genotypes_vcf@gt)
  # get the intersection
  donors <- setdiff(donors_vcf, c('FORMAT'))
  # and the genotype data itself
  gt_calls <- data.table(extract.gt(genotypes_vcf)[, donors])
  # grab the genotypes
  gt_metadata <- data.table(genotypes_vcf@fix[, c('CHROM', 'POS', 'ID', 'REF', 'ALT')])
  # create object
  genos <- Genotypes$new(calls = gt_calls, meta = gt_metadata)
  return(genos)
}


####################
# Settings.        #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 512 * 1000 * 1024^2)

# set seed
set.seed(7777)

# set cores to use
parallel::mcaffinity(1:8)

####################
# Main Code        #
####################


# location of the objects
objects_loc <- '/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/seurat/'
m1_v2_loc <- paste(objects_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')
m1_v3_loc <- paste(objects_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')

# read the objects
m1_v2 <- readRDS(m1_v2_loc)
m1_v3 <- readRDS(m1_v3_loc)

# subset to monocytes only
m1_v2_mono <- m1_v2[, !is.na(m1_v2@meta.data[['cell_type_lowerres']]) & m1_v2@meta.data[['cell_type_lowerres']] == 'monocyte']
m1_v3_mono <- m1_v3[, !is.na(m1_v3@meta.data[['cell_type_lowerres']]) & m1_v3@meta.data[['cell_type_lowerres']] == 'monocyte']

# get correlations
m1_v2_mono_cors <- get_correlations_per_condition(m1_v2_mono, 'TMEM176B', 'ELF1')
m1_v3_mono_cors <- get_correlations_per_condition(m1_v3_mono, 'TMEM176B', 'ELF1')
# add chemistry
m1_v2_mono_cors[['chem']] <- 'V2'
m1_v3_mono_cors[['chem']] <- 'V3'
# merge
m1_mono_cors <- rbind(m1_v2_mono_cors, m1_v3_mono_cors)

# get the TYK2
genotypes_loc <- '/groups/umcg-franke-scrna/tmp03/releases/wijst-2020-hg19/v1/genotype/LL_trityper_plink_converted_rs10255907.vcf.gz'
genotypes_vcf = read.vcfR(genotypes_loc)
vars <- as.vector(unlist(genotypes_vcf@gt))
# remove first entry
vars <- vars[2:length(vars)]
gts <- colnames(genotypes_vcf@gt)
gts <- gts[2:length(gts)]
# add the genotype data
m1_mono_cors[['genotype']] <- vars[match(m1_mono_cors[['val']], gts)]
# rename the conditions
m1_mono_cors[['condition']] <- as.vector(unlist(list('UT' = 'UT', 'X3hCA' = '3hCA', 'X24hCA' = '24hCA', 'X3hMTB' = '3hMTB', 'X24hMTB' = '24hMTB', 'X3hPA' = '3hPA', 'X24hPA' = '24hPA')[m1_mono_cors[['condition']]]))
# order them
m1_mono_cors[['condition']] <- factor(m1_mono_cors[['condition']], levels = c('UT', '3hCA', '3hMTB', '3hPA', '24hCA', '24hMTB', '24hPA'))

# turn into plot
p_v2 <- ggplot(data = m1_mono_cors[m1_mono_cors$chem == 'V2' & m1_mono_cors[['condition']] %in% c('UT', '3hCA', '24hCA'), ], mapping = aes(x = genotype, y = cor, fill = genotype)) + 
  geom_boxplot(outlier.shape = NA) + 
  # split by interacion term
  facet_grid(. ~ condition) +
  # and add jitter
  geom_jitter(size = 0.5, alpha = 0.5) +
  # and labels
  xlab('rs1108577') + 
  ylab('TMEM176B-ELF1 coexpression') +
  ggtitle(paste('rs1108577', 'affecting', 'TMEM176B-ELF1', 'in', 'monocytes', 'in v2')) + 
  theme(axis.text.x=element_blank(), axis.ticks = element_blank()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = list('0/0' = 'darkred', '0/1' = 'darkorange', '1/1' = 'darkblue'))

p_v3 <- ggplot(data = m1_mono_cors[m1_mono_cors$chem == 'V3' & m1_mono_cors[['condition']] %in% c('UT', '3hCA', '24hCA'), ], mapping = aes(x = genotype, y = cor, fill = genotype)) + 
  geom_boxplot(outlier.shape = NA) + 
  # split by interacion term
  facet_grid(. ~ condition) +
  # and add jitter
  geom_jitter(size = 0.5, alpha = 0.5) +
  # and labels
  xlab('rs1108577') + 
  ylab('TMEM176B-ELF1 coexpression') +
  ggtitle(paste('rs1108577', 'affecting', 'TMEM176B-ELF1', 'in', 'monocytes', 'in v3')) + 
  theme(axis.text.x=element_blank(), axis.ticks = element_blank()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  scale_fill_manual(values = list('0/0' = 'darkred', '0/1' = 'darkorange', '1/1' = 'darkblue'))

