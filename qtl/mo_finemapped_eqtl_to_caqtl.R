#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_finemapped_eqtl_to_caqtl.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(stringr)

####################
# Settings        #
####################


get_coloc_variants <- function(finemapped_colocs, cell_type_column='celltype', feature_column='Gene', variant_column='hit1', PP_H4_abf_column='PP.H4.abf', PP_H4_abf_column_cutoff=0.8, link_split_char='_') {
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


get_betas_per_celltype <- function(qtl_output_loc, cell_types=NULL, confinement_table, confinement_variant_column='variant', confinement_feature_column='gene') {
    # 

}

####################
# Main Code        #
####################

# location of the finemapped results
finemapped_colocs_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/multiome_eqtl/finemapping/coloc_e_iQTL_vs_caQTL/complete_results_all_cts.txt'
# read the tabel
finemapped_colocs <- fread(finemapped_colocs_loc, header = T, sep = '\t')
# location of the eQTL output
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# location of the caQTL output
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'
# get the colocs
significant_finemapped_colocs <- get_coloc_variants(finemapped_colocs)
