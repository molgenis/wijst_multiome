#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_finemapped_eqtl_to_caqtl.R
# Function: check the overlap of effects for colocalizing eQTLs with caQTLs
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(stringr)


####################
# Settings        #
####################

#' Extract Colocalizing Variants
#'
#' This function extracts colocalizing variants from a given dataset of finemapped colocalizations.
#'
#' @param finemapped_colocs A data frame containing finemapped colocalizations.
#' @param cell_type_column A string specifying the column name for cell types. Default is 'ct'.
#' @param feature_column A string specifying the column name for features (e.g., genes). Default is 'Gene'.
#' @param variant_column A string specifying the column name for variants. Default is 'hit1'.
#' @param PP_H4_abf_column A string specifying the column name for posterior probability (PP.H4.abf). Default is 'PP.H4.abf'.
#' @param PP_H4_abf_column_cutoff A numeric value specifying the cutoff for posterior probability. Default is 0.8.
#' @param link_split_char A string specifying the character used to split the feature link. Default is '_'.
#'
#' @return A data frame containing the colocalizing variants, genes, regions, and cell types.
get_coloc_variants <- function(finemapped_colocs, cell_type_column='ct', feature_column='Gene', variant_column='hit1', PP_H4_abf_column='PP.H4.abf', PP_H4_abf_column_cutoff=0.8, link_split_char='_') {
  # subset to significant
  finemapped_colocs_colocing <- finemapped_colocs[finemapped_colocs[[PP_H4_abf_column]] >= PP_H4_abf_column_cutoff, ]
  # extract the variant
  finemapped_colocs_colocing_short <- data.frame('variant' = finemapped_colocs_colocing[[variant_column]])
  # extract the feature link
  finemapped_colocs_colocing_short[, c('gene', 'region')] <- str_split_fixed(finemapped_colocs_colocing[[feature_column]], link_split_char, 2)
  # and the cell type
  finemapped_colocs_colocing_short[['cell_type']] <- finemapped_colocs_colocing[[cell_type_column]]
  return(finemapped_colocs_colocing_short)
}


#' Extract Beta Values for Cell Type
#'
#' This function extracts beta values for a specific cell type from a given QTL file and confinement table.
#'
#' @param cell_type_file_loc A string specifying the file location of the cell type data.
#' @param confinement_table A data frame containing the confinement table.
#' @param cell_type_name A string specifying the name of the cell type.
#' @param confinement_variant_column A string specifying the column name for variants in the confinement table. Default is 'variant'.
#' @param confinement_feature_column A string specifying the column name for features in the confinement table. Default is 'gene'.
#' @param confinement_celltype_column A string specifying the column name for cell types in the confinement table. Default is 'cell_type'.
#' @param qtl_variant_column A string specifying the column name for variants in the QTL file. Default is 'snp_id'.
#' @param qtl_feature_column A string specifying the column name for features in the QTL file. Default is 'feature_id'.
#' @param qtl_beta_column A string specifying the column name for beta values in the QTL file. Default is 'beta'.
#' @param qtl_effect_column A string specifying the column name for the assessed allele in the QTL file. Default is 'assessed_allele'.
#'
#' @return A data frame containing the cell type, variants, features, beta values, and assessed alleles.
#'
#' @examples
#' cell_type_file_loc <- "path/to/cell_type_file.txt"
#' confinement_table <- data.frame(
#'   cell_type = c('cell1', 'cell2'),
#'   var_feature = c('variant1_gene1', 'variant2_gene2')
#' )
#' get_betas_file(cell_type_file_loc, confinement_table, 'cell1')
get_betas_file <- function(cell_type_file_loc, confinement_table, cell_type_name, confinement_variant_column='variant', confinement_feature_column='gene', confinement_celltype_column='cell_type', qtl_variant_column='snp_id', qtl_feature_column='feature_id', qtl_beta_column='beta', qtl_effect_column='assessed_allele') {
  # read the file
  cell_type_table <- fread(cell_type_file_loc, header = T, sep = '\t')
  # check which var-features we want to get
  var_features_celltype <- confinement_table[confinement_table[[confinement_celltype_column]] == cell_type_name,'var_feature']
  # get this for the qtl file as well
  var_features_qtls <- paste(cell_type_table[[qtl_variant_column]], cell_type_table[[qtl_feature_column]], sep = '_')
  # filter the cell type table for these
  cell_type_table <- cell_type_table[var_features_qtls %in% var_features_celltype, ]
  # extract the data we care about
  cell_type_table_slim <- data.frame('cell_type' = rep(cell_type_name, times = nrow(cell_type_table)), 'variant' = cell_type_table[[qtl_variant_column]], 'feature' = cell_type_table[[qtl_feature_column]], 'effect' = cell_type_table[[qtl_beta_column]], 'allele' = cell_type_table[[qtl_effect_column]])
  return(cell_type_table_slim)
}


#' Extract Beta Values for Multiple Cell Types
#'
#' This function extracts beta values for multiple cell types from QTL files and a confinement table.
#'
#' @param qtl_output_loc A string specifying the directory location of the QTL output files.
#' @param confinement_table A data frame containing the confinement table.
#' @param cell_types A vector of strings specifying the cell types to include. Default is NULL, which includes all cell types.
#' @param confinement_variant_column A string specifying the column name for variants in the confinement table. Default is 'variant'.
#' @param confinement_feature_column A string specifying the column name for features in the confinement table. Default is 'gene'.
#' @param confinement_celltype_column A string specifying the column name for cell types in the confinement table. Default is 'cell_type'.
#' @param qtl_file_name A string specifying the name of the QTL file. Default is 'qtl_results_all.txt.gz'.
#' @param qtl_variant_column A string specifying the column name for variants in the QTL file. Default is 'snp_id'.
#' @param qtl_feature_column A string specifying the column name for features in the QTL file. Default is 'feature_id'.
#' @param qtl_beta_column A string specifying the column name for beta values in the QTL file. Default is 'beta'.
#' @param qtl_effect_column A string specifying the column name for the assessed allele in the QTL file. Default is 'assessed_allele'.
#'
#' @return A data frame containing the cell types, variants, features, beta values, and assessed alleles.
#'
#' @examples
#' qtl_output_loc <- "path/to/qtl_output"
#' confinement_table <- data.frame(
#'   cell_type = c('cell1', 'cell2'),
#'   variant = c('variant1', 'variant2'),
#'   gene = c('gene1', 'gene2')
#' )
#' get_betas_per_celltype(qtl_output_loc, confinement_table)
get_betas_per_celltype <- function(qtl_output_loc, confinement_table, cell_types=NULL, confinement_variant_column='variant', confinement_feature_column='gene', confinement_celltype_column='cell_type', qtl_file_name='qtl_results_all.txt.gz', qtl_variant_column='snp_id', qtl_feature_column='feature_id', qtl_beta_column='beta', qtl_effect_column='assessed_allele') {
  # get the var-gene combinations in the confinement
  confinement_table[['var_feature']] <- paste(confinement_table[[confinement_variant_column]], confinement_table[[confinement_feature_column]], sep = '_')
  # list all cell types
  cell_type_folders <- list.dirs(qtl_output_loc, recursive = F, full.names = F)
  # overlap with cell types supplied if done so
  if (!is.null(cell_types)) {
    cell_type_folders <- intersect(cell_type_folders, cell_types)
  }
  # store the cell type result in a list
  cell_types_list <- list()
  # check cell type
  for (cell_type in cell_type_folders) {
    # paste together the full path
    cell_type_file_loc <- paste0(qtl_output_loc, '/', cell_type, '/', qtl_file_name)
    message(paste('reading', cell_type_file_loc, '\n'))
    # get the betas
    cell_type_table_slim <- get_betas_file(cell_type_file_loc, 
                                           confinement_table, 
                                           cell_type_name=cell_type, 
                                           confinement_variant_column=confinement_variant_column, 
                                           confinement_feature_column=confinement_feature_column, 
                                           confinement_celltype_column=confinement_celltype_column, 
                                           qtl_variant_column=qtl_variant_column, 
                                           qtl_feature_column=qtl_feature_column, 
                                           qtl_beta_column=qtl_beta_column, 
                                           qtl_effect_column=qtl_effect_column)
    # put in the list
    cell_types_list[[cell_type]] <- cell_type_table_slim
  }
  cell_types_all <- do.call('rbind', cell_types_list)
  return(cell_types_all)
}


#' Merge Beta Values from Two Files
#'
#' This function merges beta values from two different files based on matching variants, features, and cell types.
#'
#' @param matching_table A data frame containing the matching table.
#' @param beta_table1 A data frame containing the first set of beta values.
#' @param beta_table2 A data frame containing the second set of beta values.
#' @param matching_variant_column A string specifying the column name for variants in the matching table. Default is 'variant'.
#' @param matching_feature1_column A string specifying the column name for the first feature in the matching table. Default is 'region'.
#' @param matching_feature2_column A string specifying the column name for the second feature in the matching table. Default is 'gene'.
#' @param matching_celltype_column A string specifying the column name for cell types in the matching table. Default is 'cell_type'.
#' @param variant_column1 A string specifying the column name for variants in the first beta table. Default is 'variant'.
#' @param variant_column2 A string specifying the column name for variants in the second beta table. Default is 'variant'.
#' @param feature_column1 A string specifying the column name for features in the first beta table. Default is 'feature'.
#' @param feature_column2 A string specifying the column name for features in the second beta table. Default is 'feature'.
#' @param cell_type_column1 A string specifying the column name for cell types in the first beta table. Default is 'cell_type'.
#' @param cell_type_column2 A string specifying the column name for cell types in the second beta table. Default is 'cell_type'.
#' @param beta_column1 A string specifying the column name for beta values in the first beta table. Default is 'effect'.
#' @param beta_column2 A string specifying the column name for beta values in the second beta table. Default is 'effect'.
#' @param allele_column1 A string specifying the column name for alleles in the first beta table. Default is 'allele'.
#' @param allele_column2 A string specifying the column name for alleles in the second beta table. Default is 'allele'.
#'
#' @return A data frame containing the merged beta values, along with the matching variants, features, and cell types.
#'
#' @examples
#' matching_table <- data.frame(
#'   cell_type = c('cell1', 'cell2'),
#'   variant = c('variant1', 'variant2'),
#'   region = c('region1', 'region2'),
#'   gene = c('gene1', 'gene2')
#' )
#' beta_table1 <- data.frame(
#'   cell_type = c('cell1', 'cell2'),
#'   variant = c('variant1', 'variant2'),
#'   feature = c('region1', 'region2'),
#'   effect = c(0.5, 0.6),
#'   allele = c('A', 'T')
#' )
#' beta_table2 <- data.frame(
#'   cell_type = c('cell1', 'cell2'),
#'   variant = c('variant1', 'variant2'),
#'   feature = c('gene1', 'gene2'),
#'   effect = c(0.4, 0.7),
#'   allele = c('A', 'C')
#' )
#' merge_beta_files(matching_table, beta_table1, beta_table2)
merge_beta_files <- function(
    matching_table,
    beta_table1, 
    beta_table2, 
    matching_variant_column='variant', 
    matching_feature1_column='region', 
    matching_feature2_column='gene', 
    matching_celltype_column='cell_type', 
    variant_column1='variant', 
    variant_column2='variant', 
    feature_column1='feature', 
    feature_column2='feature', 
    cell_type_column1='cell_type', 
    cell_type_column2='cell_type', 
    beta_column1='effect', 
    beta_column2='effect', 
    allele_column1='allele', 
    allele_column2='allele') {
  # get the variant-gene-celltype
  matching_set1 <- paste(matching_table[[matching_celltype_column]], matching_table[[matching_variant_column]], matching_table[[matching_feature1_column]])
  # the same for the first beta table
  beta_table1_set1 <- paste(beta_table1[[cell_type_column1]], beta_table1[[variant_column1]], beta_table1[[feature_column1]])
  # get in the right order from the set
  beta_table1_ordered <- beta_table1[match(matching_set1, beta_table1_set1), ]
  
  # get the variant-gene-celltype
  matching_set2 <- paste(matching_table[[matching_celltype_column]], matching_table[[matching_variant_column]], matching_table[[matching_feature2_column]])
  # the same for the decond beta table
  beta_table2_set2 <- paste(beta_table2[[cell_type_column2]], beta_table2[[variant_column2]], beta_table2[[feature_column2]])
  # get in the right order from the set
  beta_table2_ordered <- beta_table2[match(matching_set2, beta_table2_set2), ]
  
  # get the cell type
  cell_types <- matching_table[[matching_celltype_column]]
  # get the betas for the first
  betas1 <- beta_table1_ordered[[beta_column1]]
  # and the second one
  betas2 <- beta_table2_ordered[[beta_column2]]
  
  # drop empty entries
  indices_empty_entries <- is.na(beta_table1_ordered[[beta_column1]]) | is.na(beta_table2_ordered[[beta_column2]])
  if (length(indices_empty_entries) > 0) {
    warning(paste('removing', as.character(length(indices_empty_entries)), 'entries due to the beta being absent\n'))
    matching_table <- matching_table[!indices_empty_entries, ]
    beta_table1_ordered <- beta_table1_ordered[!indices_empty_entries, ]
    beta_table2_ordered <- beta_table2_ordered[!indices_empty_entries, ]
    cell_types <- cell_types[!indices_empty_entries]
    betas1 <- betas1[!indices_empty_entries]
    betas2 <- betas2[!indices_empty_entries]
  }
  
  # and convert beta2 if the alleles don't match
  betas2[beta_table1_ordered[[allele_column1]] != beta_table2_ordered[[allele_column2]]] <- betas2[beta_table1_ordered[[allele_column1]] != beta_table2_ordered[[allele_column2]]] * -1
  
  # make one table
  beta_comparison_table <- data.frame('effect1' = betas1, 'effect2' = betas2)
  beta_comparison_table <- cbind(matching_table, beta_comparison_table)
  return(beta_comparison_table)
}


####################
# Main Code        #
####################

# location of the finemapped results
finemapped_colocs_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/multiome_eqtl/finemapping/coloc_e_iQTL_vs_caQTL/complete_results_all_cts.txt'
# read the tabel
finemapped_colocs <- fread(finemapped_colocs_loc, header = T, sep = '\t')
# location of the eQTL output
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/UT/'
eqtl_file_name <- 'qtl_results_all_nominally_significant.txt.gz'
# location of the caQTL output
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/combined_output_50kb/FDR_all_effects/'
# get the colocs
significant_finemapped_colocs <- get_coloc_variants(finemapped_colocs)
# get the eqtl betas
eqtl_betas <- get_betas_per_celltype(eqtl_output_loc, significant_finemapped_colocs, qtl_file_name = eqtl_file_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'))
eqtl_betas_unique <- unique(data.frame(eqtl_betas))
# check each cell type
finemapped_colocs_varfeature <- significant_finemapped_colocs
finemapped_colocs_varfeature[['var_feature']] <- paste(finemapped_colocs_varfeature[['variant']], finemapped_colocs_varfeature[['region']], sep = '_')
caqtl_betas_list <- list()
for (cell_type in c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')) {
  # read file
  caqtl_betas_list[[cell_type]] <- get_betas_file(paste0(caqtl_output_loc, '/UT_', cell_type, '_conditions_speific_0001.tsv'), 
                                                  confinement_table = finemapped_colocs_varfeature, 
                                                  confinement_variant_column='variant', 
                                                  confinement_feature_column='gene', 
                                                  confinement_celltype_column='cell_type', 
                                                  cell_type_name=cell_type, 
                                                  qtl_variant_column='snp_id', 
                                                  qtl_feature_column='feature_id', 
                                                  qtl_beta_column='beta', 
                                                  qtl_effect_column='assessed_allele'
  )
}
# merge
caqtl_betas <- do.call('rbind', caqtl_betas_list)
caqtl_betas_unique <- unique(data.frame(caqtl_betas))


# # add cell type to region to variant
# finemapped_colocs_varfeature[['var_feat_ct']] <- paste(finemapped_colocs_varfeature[['variant']], finemapped_colocs_varfeature[['region']], finemapped_colocs_varfeature[['cell_type']])
# # get the full stats
# all_ca_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/combined_output_50kb/FDR_all_effects/ALL_UT__conditions_speific_0001.tsv'
# all_ca <- fread(all_ca_loc, header = T, sep = '\t')
# # add the same in the caQTL table
# all_ca[['var_feat_ct']] <- paste(all_ca[['snp_id']], all_ca[['feature_id']], all_ca[['ct']])
# # now get those specifically
# all_ca_ordered <- all_ca[match(finemapped_colocs_varfeature[['var_feat_ct']], all_ca[['var_feat_ct']]), ]
# # and get the columns we care about
# caqtl_betas <- data.frame('cell_type' = all_ca_ordered[['ct']], 'variant' = all_ca_ordered[['snp_id']], 'feature' = all_ca_ordered[['feature_id']], 'effect' = all_ca_ordered[['beta']], 'allele' = all_ca_ordered[['assessed_allele']])

# remove empty entries
significant_finemapped_colocs <- significant_finemapped_colocs[significant_finemapped_colocs[['variant']] != '-', ]
# get eqtl and caqtl betas
betas_both <- merge_beta_files(significant_finemapped_colocs, caqtl_betas, eqtl_betas)

sum(sign(betas_both$effect1)  == sign(betas_both$effect2)) / nrow(betas_both)
# 0.9270331
cor(betas_both$effect1, betas_both$effect2)
# 0.8112404

# do each cell type
for (cell_type in unique(betas_both[['cell_type']])) {
  # get results for that table
  betas_both_ct <- betas_both[betas_both[['cell_type']] == cell_type, ]
  # and do the same step
  print(cell_type)
  print(paste('con:', sum(sign(betas_both_ct$effect1)  == sign(betas_both_ct$effect2)) / nrow(betas_both_ct)))
  print(paste('cor:', cor(betas_both_ct$effect1, betas_both_ct$effect2)))
}
write.table(betas_both, gzfile('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/mediation/coloc/mo_colocing_betas.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)

# read the file again
betas_both <- read.table('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/mediation/coloc/mo_colocing_betas.tsv.gz', header = T, sep = '\t')

# get the min and max values on the axis
max_beta1 <- max(abs(betas_both[['effect1']]))
min_beta1 <- min(abs(betas_both[['effect1']]))
max_beta2 <- max(abs(betas_both[['effect2']]))
min_beta2 <- min(abs(betas_both[['effect2']]))

# start making the plot
plot(x = betas_both[['effect1']], y = betas_both[['effect2']], 
     ylim = c(max_beta2 * -1, max_beta2),
     xlim = c(max_beta1 * -1, max_beta1),
     xlab = 'beta caQTL',
     ylab = 'beta eQTL',
     pch = 16,
     cex = 0.5,
     main = 'betas of caQTLs versus eQTLs'
)
# the top to bottom rectangle of non-significant effects
rect(xleft = -1 * min_beta1, xright = min_beta1, ybottom = -1 * max_beta2, ytop = max_beta2, border = NA, col = rgb(red = 1, green = 1, blue = 1, alpha = 0.5))
# the left to right rectangle of non-significant effects
rect(xleft = -1 * max_beta1, xright = max_beta1, ybottom = -1 * min_beta2, ytop = min_beta2, border = NA, col = rgb(red = 1, green = 1, blue = 1, alpha = 0.5))
# bottom left concordant
rect(xleft = -1 * max_beta1, xright = -1 *min_beta1, ybottom = -1 * max_beta2, ytop = -1 * min_beta2, border = NA, col = rgb(red = 0, green = .45, blue = .7, alpha = 0.5))
# bottom right disconcordant
rect(xleft = min_beta1, xright = max_beta1, ybottom = -1 * max_beta2, ytop = -1 * min_beta2, border = NA, col = rgb(red = .8, green = .4, blue = 0, alpha = 0.5))
# top left disconcordant
rect(xleft = -1 * max_beta1, xright = -1 *min_beta1, ybottom = min_beta2, ytop = max_beta2, border = NA, col = rgb(red = .8, green = .4, blue = 0, alpha = 0.5))
# top right disconcordant
rect(xleft = min_beta1, xright = max_beta1, ybottom = min_beta2, ytop = max_beta2, border = NA, col = rgb(red = 0, green = .45, blue = .7, alpha = 0.5))
# box to put concordance label in
rect(xleft = max_beta2 * 0.60, max_beta2 * 0.95, ybottom = max_beta2 * -0.8, ytop = max_beta2 * -0.70, col = 'white')
# add concordance label
text(x = max_beta1 * 0.75 , y = max_beta2 * -0.75, labels = c(paste('concordance', 0.92, sep = ':\n')))
# add the annotation of what if concordant and disconcordant
#text(x = max_sig_z_venema * 0.75 , y = max_sig_z_eqtlgen * 0.75, labels = c('concordant'), col = rgb(red = 0, green = .45, blue = .7))
#text(x = max_sig_z_venema * - 0.75 , y = max_sig_z_eqtlgen * 0.75, labels = c('disconcordant'), col = rgb(red = .8, green = .4, blue = 0))
# add dashed lines
lines(c(max_beta1 * -1, max_beta1), c(min_beta2 * -1, min_beta2 * -1), type = "l", lty = 2)
lines(c(max_beta1 * -1, max_beta1), c(min_beta2, min_beta2), type = "l", lty = 2)
lines(c(min_beta1 * -1, min_beta1 * -1), c(max_beta2 * -1, max_beta2), type = "l", lty = 2)
lines(c(min_beta1, min_beta1), c(max_beta2 * -1, max_beta2), type = "l", lty = 2)

