#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_monique_conference_plots.R
# Function: plots for Monique her plant conference
############################################################################################################################


####################
# libraries        #
####################

library(ggplot2)
# for plink file loading
library(snpStats)
# for regression models
library(lme4)
# Seurat and Signac
library(Seurat)
library(Signac)
# data table
library(data.table)


####################
# Functions        #
####################

#' Get PLINK Genotypes for a Chromosome
#'
#' This function reads PLINK genotype files for a specified chromosome.
#'
#' @param genotype_loc Character. The location of the genotype files.
#' @param chromosome Integer. The chromosome number.
#' @param genotype_prepend Character. The prefix for the genotype files. Default is 'EUR_imputed_hg38_varFiltered_chr'.
#' @param genotype_append Character. The suffix for the genotype files. Default is ''.
#' @return A list containing the genotype data.
#' @export
#' @examples
#' genotypes <- get_plink_genotypes_chromosome("path/to/genotypes", 1)
get_plink_genotypes_chromosome <- function(genotype_loc, chromosome, genotype_prepend='EUR_imputed_hg38_varFiltered_chr', genotype_append='') {
  genotypes_loc <- paste(genotype_loc, genotype_prepend, chromosome, genotype_append, sep = '')
  genotypes <- read.plink(
    bed = paste(genotypes_loc, '.bed', sep = ''),
    bim = paste(genotypes_loc, '.bim', sep = ''),
    fam = paste(genotypes_loc, '.fam', sep = '')
  )
  return(genotypes)
}


get_peak_info <- function(peaks_object) {
  # calculate average expression
  avg_peaks <- AverageExpression(peaks_object)[['peaks']]
  # calculate average expression
  sum_peaks <- AggregateExpression(peaks_object)[['peaks']]
  # calculate pct exp
  counts_as_row_sparse <- as(peaks_object@assays$peaks@counts, "RsparseMatrix") # turn into row-wise sparse matrix
  non_zero_cell_nr <- as.vector(unlist(rowSums(counts_as_row_sparse != 0))) # do a rowsum on the T/F value you get from the zero-or-not comparison
  pct_peaks <- non_zero_cell_nr / ncol(peaks_object@assays$peaks@counts) # the percentage expressed is that number of cells divided by the total number
  # turn the exp per group into a bed
  peaks_bed <- get_peaks_to_bed(sum_peaks)[[1]]
  # add the other info
  peaks_bed[['avg']] <-  unlist(as.vector(avg_peaks[, 1]))
  peaks_bed[['ncell']] <-  ncol(peaks_object)
  peaks_bed[['pct_exp']] <- pct_peaks
  return(peaks_bed)
}


get_peaks_to_bed <- function(exp_per_group) {
  bed_per_group <- list()
  # check each of the columns
  for (identity in colnames(exp_per_group)) {
    print(identity)
    # subset to that identity
    peaks_ident <- exp_per_group[, c(identity), drop = F]
    # get the locations
    positions <- data.frame(do.call('rbind', strsplit(rownames(peaks_ident), '-')))
    # set the colnames properly to bed format
    colnames(positions) <- c('#chrom', 'start', 'end')
    # add the rownames themselves as the name
    positions[['name']] <- rownames(peaks_ident)
    # and the counts as score
    positions[['exp']] <- log10(peaks_ident[, c(identity)])
    # make the start and stop numeric
    positions[['start']] <- as.numeric(positions[['start']])
    positions[['end']] <- as.numeric(positions[['end']])
    # add bed to list
    bed_per_group[[identity]] <- positions
  }
  return(bed_per_group)
}


get_coexpression_per_sample <- function(multimodal_object, gene, region, sample_column) {
  # get the unique samples
  samples_unique <- unique(multimodal_object@meta.data[[sample_column]])
  # remove NA
  samples_unique <- samples_unique[!is.na(samples_unique)]
  # make a table of the correlations
  cor_table <- matrix(data = NA, nrow = length(samples_unique), ncol = 2)
  # set column names
  colnames(cor_table) <- c('sample', 'correlation')
  # set the data types
  cor_table[[1]] <- as(cor_table[[1]], 'character')
  cor_table[[2]] <- as(cor_table[[2]], 'numeric')
  # check each sample by index
  for (i in 1 : length(samples_unique)) {
    # extract the sample
    sample_at_i <- samples_unique[i]
    # subset the data to that sample
    object_sample <- multimodal_object[, !is.na(multimodal_object@meta.data[[sample_column]]) & multimodal_object@meta.data[[sample_column]] == sample_at_i]
    # calculate the correlation
    cor_sample <- cor(x = multimodal_object@assays$RNA@scale.data[gene, ], y = multimodal_object@assays$peaks@counts[region, ], method = 'spearman')
    # put into the matrix
    cor_table[i, 'sample'] <- sample_at_i
    cor_table[i, 'correlation'] <- cor_sample
  }
  # return the matrix
  return(cor_table)
}

####################
# Main Code        #
####################

# cre-QTL example
tbl <- read.table('/groups/umcg-franke-scrna/tmp04/users/umcg-roelen/singularity/rstudio-server/simulated_home/monique_creqtl_example.tsv.gz', header = T, sep = '\t')
ggplot(data = tbl[tbl$omics2 > -.971, ], 
       mapping = aes(
         x = omics2, 
         y = omics1, 
         col = alleles)
       ) + 
  geom_point() + 
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", 
              se = F) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_colour_manual(values = list('AA' = 'darkblue', 'AT' = 'darkorange', 'TT' = 'darkgreen'), name = 'genotype') +
  xlab('chr12-8951819-8952176 accessibility') +
  ylab('PZP expression')

# eQTL example
ggplot(data = tbl, 
       mapping = aes(
         x = alleles, 
         y = omics1, 
         fill = alleles)
) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 0.5, alpha = 0.5) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = list('AA' = 'darkblue', 'AT' = 'darkorange', 'TT' = 'darkgreen'), name = 'genotype') +
  xlab('genotype') +
  ylab('PZP expression')

# get interaction p-value of CRE-QTL
nogeno <- do.call(what = 'lmer', list(formula = 'omics1~omics2+(1|genotype_id)', data = tbl))
yesgeno <- do.call(what = 'lmer', list(formula = 'omics1~omics2+(1|genotype_id)+genotype*omics2', data = tbl))
anova((nogeno), (yesgeno), alternative = "less")

# load genotype data
genotypes_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/'
genotypes <- get_plink_genotypes_chromosome(genotypes_loc, '12')
  
# get co-expression-QTL datat
cor_tbl <- read.table('/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/calculate_correlation_metrics/output/wg3_multiome/wg3_multiome-CD4_T-pearson-weighted-chr-12-final.tsv.gz', header = T, sep = '\t', row.names = 1)
cor_tbl_plottable <- data.frame('donor' = colnames(cor_tbl), 'correlation' = as.vector(unlist(cor_tbl['BACH2_PZP', ])))
cor_tbl_plottable[['alleles']] <- tbl[match(cor_tbl_plottable$donor, tbl$genotype_id), 'alleles']
cor_tbl_plottable[['genotype']] <- tbl[match(cor_tbl_plottable$donor, tbl$genotype_id), 'genotype']
# get p-value of co-eQTL
summary(lm(as.formula('correlation~genotype'), data = cor_tbl_plottable))
# get spearman R of co-eQTL
cor(
  x = cor_tbl_plottable[!is.na(cor_tbl_plottable[['correlation']]) & 
                          !is.na(cor_tbl_plottable[['genotype']]), ][['correlation']], 
  y = cor_tbl_plottable[!is.na(cor_tbl_plottable[['correlation']]) & 
                          !is.na(cor_tbl_plottable[['genotype']]), ][['genotype']], 
  method = 'spearman'
)
# make co-eQTL plot
ggplot(data = cor_tbl_plottable[!is.na(cor_tbl_plottable$alleles), ], 
       mapping = aes(
         x = alleles, 
         y = correlation, 
         fill = alleles)
) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 0.5, alpha = 0.5) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = list('AA' = 'darkblue', 'AT' = 'darkorange', 'TT' = 'darkgreen'), name = 'genotype') +
  xlab('Genotype') +
  ylab('PZP-BACH2 Spearman correlation') + 
  ggtitle('rs7299653 affecting PZP-BACH2 co-expression') +
  theme(legend.position = 'none')

# load metadata of expression data
expression_metadata <- fread('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_celllevel_metadata.tsv.gz', header = T, sep = '\t')
# get the CD4 data
cd4_counts <- data.frame(table(expression_metadata[expression_metadata[['predicted.mo_10x_cell_type']] %in% c('CD4 TCM', 'CD4 Naive', 'CD4 TEM', 'Treg'), c('sample_final', 'lane', 'condition_final', 'predicted.mo_10x_cell_type')]))
# remove empty entries
cd4_counts <- cd4_counts[cd4_counts$Freq > 0, ]
# sum the memory cells for everyone
cd4_mem_to_naive_counts_list <- list()
# check each donor
for (donor in unique(cd4_counts[['sample_final']])) {
  # check each lane
  for (lane in unique(cd4_counts[cd4_counts[['sample_final']] == donor, 'lane'])) {
    # extract that entry
    lane_donor_entry <- cd4_counts[
      cd4_counts[['sample_final']] == donor &
        cd4_counts[['lane']] == lane,
    ]
    # get frequencies of CD4 naive
    lane_donor_cd4_naive_freq <- sum(lane_donor_entry[lane_donor_entry[['predicted.mo_10x_cell_type']] == 'CD4 Naive', 'Freq'])
    # and of memory
    lane_donor_cd4_memory_freq <- sum(lane_donor_entry[lane_donor_entry[['predicted.mo_10x_cell_type']] %in% c('CD4 TCM', 'CD4 TEM'), 'Freq'])
    # calculate the lfc
    lfc <- log2(lane_donor_cd4_naive_freq/lane_donor_cd4_memory_freq)
    # make a row
    lane_donor_row <- data.frame('sample_final' = c(donor), 'lane' = c(lane), 'condition_final' = unique(lane_donor_entry[['condition_final']])[1], 'CD4_naive_vs_memory_lfc' = c(lfc))
    # put in the list
    cd4_mem_to_naive_counts_list[[paste(donor, lane, sep = ';;')]] <- lane_donor_row
  }
}
# merge all
cd4_mem_to_naive_counts <- do.call('rbind', cd4_mem_to_naive_counts_list)
# add genotype
cd4_mem_to_naive_counts[['genotype_coded']] <- as.character(genotypes$genotypes[cd4_mem_to_naive_counts[['sample_final']], '12:8906397:A:T'])
cd4_mem_to_naive_counts[['genotype']] <- as.vector(unlist(list('00' = NA, '01' = 0, '02' = 1, '03' = 2)[cd4_mem_to_naive_counts[['genotype_coded']]]))
cd4_mem_to_naive_counts[['alleles']] <- tbl[match(cd4_mem_to_naive_counts$sample_final, tbl$genotype_id), 'alleles']
cd4_mem_to_naive_counts[['alleles']] <- as.vector(unlist(list('00' = NA, '01' = 'AA', '02' = 'AT', '03' = 'TT')[as.character(cd4_mem_to_naive_counts[['genotype_coded']])]))

# make co-eQTL plot
ggplot(data = cd4_mem_to_naive_counts[!is.na(cd4_mem_to_naive_counts$alleles) & cd4_mem_to_naive_counts[['condition_final']] == '24hCA', ], 
       mapping = aes(
         x = alleles, 
         y = CD4_naive_vs_memory_lfc, 
         fill = alleles)
) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size = 0.5, alpha = 0.5) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) +
  scale_fill_manual(values = list('AA' = 'darkblue', 'AT' = 'darkorange', 'TT' = 'darkgreen'), name = 'genotype') +
  xlab('Genotype') +
  ylab('CD4+ T Memory versus Naive LFC') + 
  ggtitle('rs7299653 affecting CD4+ T Memory to Naive fraction in 24hCA') +
  theme(legend.position = 'none')

# location of the CD4T+ cells
cd4t_acc_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_cd4t_wstatus_1_80_20240709.rds'
# read the object
cd4t <- readRDS(cd4t_acc_loc)
# add previous metadata
cd4t@meta.data[['predicted.mo_10x_cell_type']] <- as.vector(unlist(expression_metadata[match(colnames(cd4t), expression_metadata[['barcode_lane']]), 'predicted.mo_10x_cell_type']))
cd4t_naive_acc <- get_peak_info(cd4t[, !is.na(cd4t@meta.data[['predicted.mo_10x_cell_type']]) & cd4t@meta.data[['predicted.mo_10x_cell_type']] == 'CD4 Naive'])
cd4t_mem_acc <- get_peak_info(cd4t[, !is.na(cd4t@meta.data[['predicted.mo_10x_cell_type']]) & (cd4t@meta.data[['predicted.mo_10x_cell_type']] %in% c('CD4 TCM', 'CD4 TEM'))])
# for conditions
cd4t_naive_acc_ut <- get_peak_info(cd4t[, !is.na(cd4t@meta.data[['predicted.mo_10x_cell_type']]) & 
                                       cd4t@meta.data[['predicted.mo_10x_cell_type']] == 'CD4 Naive' &
                                       !is.na(cd4t@meta.data[['condition_final']]) &
                                         cd4t@meta.data[['condition_final']] == 'UT', ])
cd4t_mem_acc_ut <- get_peak_info(cd4t[, !is.na(cd4t@meta.data[['predicted.mo_10x_cell_type']]) & 
                                        (cd4t@meta.data[['predicted.mo_10x_cell_type']] %in% c('CD4 TCM', 'CD4 TEM')) &
                                     !is.na(cd4t@meta.data[['condition_final']]) &
                                     cd4t@meta.data[['condition_final']] == 'UT', ])
cd4t_naive_acc_24hca <- get_peak_info(cd4t[, !is.na(cd4t@meta.data[['predicted.mo_10x_cell_type']]) & 
                                          cd4t@meta.data[['predicted.mo_10x_cell_type']] == 'CD4 Naive' &
                                          !is.na(cd4t@meta.data[['condition_final']]) &
                                          cd4t@meta.data[['condition_final']] == '24hCA', ])
cd4t_mem_acc_24hca <- get_peak_info(cd4t[, !is.na(cd4t@meta.data[['predicted.mo_10x_cell_type']]) & 
                                        (cd4t@meta.data[['predicted.mo_10x_cell_type']] %in% c('CD4 TCM', 'CD4 TEM')) &
                                        !is.na(cd4t@meta.data[['condition_final']]) &
                                        cd4t@meta.data[['condition_final']] == '24hCA', ])

# read the multiome object as well
cd4t_multiome_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_cd4t_1_80_20240521.rds'
cd4t_multiome <- readRDS(cd4t_multiome_loc)
# get base correlation
cor(x = cd4t_multiome@assays$RNA@scale.data['PZP', ], y = cd4t_multiome@assays$peaks@counts['chr12-8951819-8952176', ], method = 'spearman')
# add previous metadata
cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']] <- as.vector(unlist(expression_metadata[match(colnames(cd4t_multiome), expression_metadata[['barcode_lane']]), 'predicted.mo_10x_cell_type']))
# also add the sample with lane
cd4t_multiome@meta.data[['sample_lane']] <- paste(cd4t_multiome@meta.data[['sample_final']], cd4t_multiome@meta.data[['lane']], sep = ';;')
# cell subtype and condition
cd4t_mem_24hca <- cd4t_multiome[, !is.na(cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']]) & 
       (cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']] %in% c('CD4 TCM', 'CD4 TEM')) &
       !is.na(cd4t_multiome@meta.data[['condition_final']]) &
         cd4t_multiome@meta.data[['condition_final']] == '24hCA', ]
cd4t_naive_24hca <- cd4t_multiome[, !is.na(cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']]) & 
                                    cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']] == 'CD4 Naive' &
                           !is.na(cd4t_multiome@meta.data[['condition_final']]) &
                             cd4t_multiome@meta.data[['condition_final']] == '24hCA', ]
cd4t_naive_ut <- cd4t_multiome[, !is.na(cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']]) & 
                                 cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']] == 'CD4 Naive' &
       !is.na(cd4t_multiome@meta.data[['condition_final']]) &
         cd4t_multiome@meta.data[['condition_final']] == 'UT', ]
cd4t_mem_ut <- cd4t_multiome[, !is.na(cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']]) & 
       (cd4t_multiome@meta.data[['predicted.mo_10x_cell_type']] %in% c('CD4 TCM', 'CD4 TEM')) &
       !is.na(cd4t_multiome@meta.data[['condition_final']]) &
         cd4t_multiome@meta.data[['condition_final']] == 'UT', ]
# correlations
cor(x = cd4t_mem_24hca@assays$RNA@scale.data['PZP', ], y = cd4t_mem_24hca@assays$peaks@counts['chr12-8951819-8952176', ], method = 'spearman')
cor(x = cd4t_naive_24hca@assays$RNA@scale.data['PZP', ], y = cd4t_naive_24hca@assays$peaks@counts['chr12-8951819-8952176', ], method = 'spearman')
cor(x = cd4t_naive_ut@assays$RNA@scale.data['PZP', ], y = cd4t_naive_ut@assays$peaks@counts['chr12-8951819-8952176', ], method = 'spearman')
cor(x = cd4t_mem_ut@assays$RNA@scale.data['PZP', ], y = cd4t_mem_ut@assays$peaks@counts['chr12-8951819-8952176', ], method = 'spearman')
# correlations per sample
sample_cors_cd4t <- get_coexpression_per_sample(cd4t_multiome, 'PZP', 'chr12-8951819-8952176', 'sample_lane')
