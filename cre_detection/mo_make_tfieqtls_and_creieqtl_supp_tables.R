#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_make_tfieqtls_and_creieqtl_supp_tables.R
# Function: make the supplementary tables for the TFa-i-eQTLs and the CRE-i-eQTLs
############################################################################################################################


####################
# libraries        #
####################

# data table
library(data.table)
# for making the Excel file
library(xlsx)


####################
# Functions        #
####################


####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- T


###################
# interested data #
###################

# set location of both TF tables
tfa_ieqtls_annotated_loc <- paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/', 'results_fdr_with_replication.tsv.gz', sep = '/')
tfe_ieqtls_loc <- paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion_tfexpressions/merged/', 'results_fdr.tsv.gz', sep = '/')
# read these
tfa_ieqtls <- fread(tfa_ieqtls_annotated_loc, header = T, sep = '\t')
tfe_ieqtls <- fread(tfe_ieqtls_loc, header = T, sep = '\t')
# rename columns
colnames(tfa_ieqtls) <- gsub('tf', 'tfa', colnames(tfa_ieqtls))
colnames(tfe_ieqtls) <- gsub('region', 'tfe', colnames(tfe_ieqtls))
# add explicit tf column
tfa_ieqtls[['tf']] <- gsub("_(extended|direct).*", "", tfa_ieqtls[['tfa']])
# which will just be the same for the tfe
tfe_ieqtls[['tf']] <- tfe_ieqtls[['tfe']]
# one column will be the same, replace it
colnames(tfe_ieqtls) <- gsub('anova', 'anova_tfe', colnames(tfe_ieqtls))
# add the tfe info to the TFa table
tfa_ieqtls <- merge(
  tfa_ieqtls, 
  tfe_ieqtls[, c('variant', 'tf', 'gene', 'tfe_beta', 'tfe_se', 'tfe_p', 'tfe:genotype_beta', 'tfe:genotype_se', 'tfe:genotype_p', 'tfe_bh', 'tfe:genotype_bh', 'anova_tfe_p', 'anova_tfe_bh')], 
  by = c('variant', 'tf', 'gene'), 
  all.x = T, 
  all.y = F
)
tfa_ieqtls_columns_to_keep <- setdiff(colnames(tfa_ieqtls), c('tfa_tval', 'genotype_tval', 'tfa:genotype_tval', 'chunk', 'family', 'anova', 'tfa_bf', 'genotype_bf', 'anova_bf', 'tfa:genotype_bf', 'gene_tss', 'gene_strand', 'peak_chromosome', 'peak_start', 'peak_end', 'is_ieqtl', 'is_tfaqtl', 'significant'))
# remove some columns to reduce crowdedness
tfa_ieqtls <- tfa_ieqtls[, ..tfa_ieqtls_columns_to_keep]
# annotate the columns
tfa_columns <- list(
  "variant" = "genetic variant tested", 
  "tf" = "transcription factor eregulon describes", 
  "gene" = "gene tested", 
  "tfa" = "eregulon/TF activity tested", 
  "tfa_beta" = "baseline beta of TF activity", 
  "tfa_se" = "baseline standard error of TF activity", 
  "tfa_df" = "degrees of freedom for TF activity", 
  "tfa_p" = "baseline nominal p value of TF activity", 
  "genotype_beta" = "baseline beta of genotype", 
  "genotype_se" = "baseline standard error of genotype", 
  "genotype_df" = "degrees of freedom for genotype", 
  "genotype_p" = "baseline nominal p value of genotype", 
  "tfa:genotype_beta" = "interaction beta of TF activity*genotype interaction term", 
  "tfa:genotype_se" = "interaction standard error of TF activity*genotype interaction term", 
  "tfa:genotype_df" = "degrees of freedom for TF activity*genotype interaction term", 
  "tfa:genotype_p" = "nominal p value for TF activity*genotype interaction term", 
  "anova_p" = "nominal p value of F-test comparing baseline vs interaction model for TF activity", 
  "ncell" = "number of cells used for this test", 
  "nparticipant" = "number of unique genotype samples used for this test", 
  "tfa_bh" = "B&H corrected p value for baseline TF activity", 
  "genotype_bh" = "B&H corrected p value for baseline genotype", 
  "tfa:genotype_bh" = "B&H corrected p value for activity*genotype interaction term", 
  "anova_bh" = "B&H corrected p value of F-test comparing baseline vs interaction model for TF activity", 
  "gene_chromosome" = "chromosome the gene is located", 
  "gene_start" = "chromosomal start position of gene", 
  "gene_end" = "chromosomal end position of gene", 
  "gene_tt" = "type of transcript the gene is", 
  "snp_chromosome" = "chromosome the variant is located", 
  "snp_position" = "chromosomal", 
  "is_eqtl" = "the variant and gene are an eQTL in these cell types", 
  "is_caqtl_variant" = "the variant is a caQTL in these cell types", 
  "tfa:genotype_p_rep" = "the nominal p value of the replication run removing the i-eGene from the TF activity", 
  "tfa:genotype_bh_rep" = "the B&H corrected p value of the replication run removing the i-eGene from the TF activity", 
  "anova_p_rep" = "nominal p value of F-test comparing baseline vs interaction model, of the replication run removing the i-eGene from the TF activity", 
  "anova_bh_rep" = "B&H corrected p value of F-test comparing baseline vs interaction model, of the replication run removing the i-eGene from the TF activity", 
  "tfa:genotype_beta_rep" = "interaction beta of TF activity*genotype interaction term, of the replication run removing the i-eGene from the TF activity", 
  "tfa:genotype_se_rep" = "interaction standard error of TF activity*genotype interaction term, of the replication run removing the i-eGene from the TF activity", 
  "tfa_bh_rep" = "B&H corrected p value for baseline TF activity, of the replication run removing the i-eGene from the TF activity", 
  "tfa_gene_direction" = "whether the TF activity is a repressor (-/+) or an activator (+/+)", 
  "replication_flip" = "whether the TF activity was flipped in the replication, due to being a repressor", 
  "is_cre_ieqtl" = "whether the variant and gene also are part of a CRE-i-eQTL", 
  "tfe_beta" = "baseline beta of TF expression", 
  "tfe_se" = "baseline standard error of TF expression", 
  "tfe_p" = "baseline nominal p value of TF expression", 
  "tfe:genotype_beta" = "interaction beta of TF expression*genotype interaction term", 
  "tfe:genotype_se" = "interaction standard error of TF expression*genotype interaction term", 
  "tfe:genotype_p" = "nominal p value of TF expression*genotype interaction term", 
  "tfe_bh" = "B&H corrected p value of baseline TF expression", 
  "tfe:genotype_bh" = "B&H corrected p value of TF expression*genotype interaction term", 
  "anova_tfe_p" = "nominal p value of F-test comparing baseline vs interaction model for TF expression", 
  "anova_tfe_bh" = "B&H corrected p value of F-test comparing baseline vs interaction model for TF expression"
)
# into df
tfa_columns_df <- data.frame(
  'column' = names(tfa_columns), 
  'description' = as.vector(unlist(tfa_columns[names(tfa_columns)]))
)
# create a new workbook
wb_tfa = createWorkbook()
# create the sheet
tfa_sheet_data = createSheet(wb_tfa, 'TFa_i_eQTL_sumstats')
# add dataframe of actual data to sheet
addDataFrame(tfa_ieqtls, sheet = tfa_sheet_data, startColumn=1, row.names=FALSE)
# create sheet of variable explanations
tfa_sheet_varexps = createSheet(wb_tfa, 'TFa_i_eQTL_variable_explanations')
# add dataframe of variable explanations to sheet
addDataFrame(tfa_columns_df, sheet = tfa_sheet_varexps, startColumn=1, row.names=FALSE)
# set the output location
tfa_excel_output_loc <- paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/', 'results_fdr_with_replication_supptbl.xlsx', sep = '/')
# write the result
saveWorkbook(wb_tfa, tfa_excel_output_loc)


# location of the region-i-ieqtls
region_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/region_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz'
# read the inputs
region_ieqtls <- fread(region_ieqtl_loc, header = T, sep = '\t')
# replace 'region' with 'cre'
colnames(region_ieqtls) <- gsub('region', 'cre', colnames(region_ieqtls))
# remove columns to make smaller
region_ieqtls_columns_to_keep <- setdiff(colnames(region_ieqtls), c('cre_tval', 'genotype_tval', 'cre:genotype_tval', 'family', 'anova', 'chunk', 'cre_bf', 'genotype_bf', 'cre:genotype_bf', 'anova_bf'))
# keep those columns
region_ieqtls <- region_ieqtls[, ..region_ieqtls_columns_to_keep]
# annotate the columns
rieqtl_columns <- list(
  "variant" = "genetic variant tested", 
  "cre" ="Cis Regulatory Element tested", 
  "gene" = "gene tested", 
  "cre_beta" = "baseline beta of CRE accessibility", 
  "cre_se" = "baseline standard error of CRE accessibility", 
  "cre_df" = "degrees of freedom for CRE accessibility", 
  "cre_p" = "baseline nominal p value of CRE accessibility", 
  "genotype_beta" = "baseline beta of genotype", 
  "genotype_se" = "baseline standard error of genotype", 
  "genotype_df" = "degrees of freedom for genotype", 
  "genotype_p" = "baseline nominal p value of genotype", 
  "cre:genotype_beta" = "interaction beta of CRE accessibility*genotype interaction term", 
  "cre:genotype_se " ="standard error of CRE accessibility*genotype interaction term", 
  "cre:genotype_df" = "degrees of freedom for CRE accessibility*genotype interaction term", 
  "cre:genotype_p" = "nominal p value for CRE accessibility*genotype interaction term", 
  "anova_p" = "nominal p value of F-test comparing baseline vs interaction model", 
  "ncell" = "number of cells used for this test", 
  "nparticipant" = "number of unique genotype samples used for this test", 
  "cre_bh" = "baseline B&H corrected p value of CRE accessibility", 
  "genotype_bh" = "baseline B&H corrected p value of genotype", 
  "cre:genotype_bh" = "B&H corrected p value for CRE accessibility*genotype interaction term", 
  "anova_bh" = "B&H corrected p value of F-test comparing baseline vs interaction model"
)
# into df
rieqtl_columns_df <- data.frame(
  'column' = names(rieqtl_columns), 
  'description' = as.vector(unlist(rieqtl_columns[names(rieqtl_columns)]))
)
# create a new workbook
wb_rieqtl = createWorkbook()
# create the sheet
rieqtl_sheet_data = createSheet(wb_rieqtl, 'CRE_i_eQTL_sumstats')
# add dataframe of actual data to sheet
addDataFrame(region_ieqtls, sheet = rieqtl_sheet_data, startColumn=1, row.names=FALSE)
# create sheet of variable explanations
rieqtl_sheet_varexps = createSheet(wb_rieqtl, 'CRE_i_eQTL_variable_explanations')
# add dataframe of variable explanations to sheet
addDataFrame(rieqtl_columns_df, sheet = rieqtl_sheet_varexps, startColumn=1, row.names=FALSE)
# set the output location
rieqtl_excel_output_loc <- paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/region_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/', 'results_fdr_supptbl.xlsx', sep = '/')
# write the result
saveWorkbook(wb_rieqtl, rieqtl_excel_output_loc)


# location of the SCENIC outout
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both_filtered.tsv.gz'
# read that
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# remove what is not plausible
scenic_output <- scenic_output[scenic_output[['Gene_signature_direction']] %in% c('-/+', '+/+'), ]
# get which are non-extended if it was both extended and non-extended
scenic_eregs <- unique(scenic_output[, c('TF', 'Gene_signature_direction', 'Gene_signature_name', 'source')])
scenic_eregs <-scenic_eregs[order(scenic_eregs[['source']]), ]
scenic_eregs_to_keep <- scenic_eregs[!duplicated(paste(scenic_eregs[['TF']], scenic_eregs[['Gene_signature_direction']])), ]
# then use that to filer
scenic_output <- scenic_output[scenic_output[['Gene_signature_name']] %in% scenic_eregs_to_keep[['Gene_signature_name']], ]
# remove columns to make it smaller
scenic_output_columns_to_keep <- setdiff(colnames(scenic_output), c('is_extended', 'topics_membership', 'topics_dar', 'caqtl_celltype', 'eqtl_celltype', 'r2g', 'eqtl_varregion_celltype'))
# keep those columns
scenic_output <- scenic_output[, ..scenic_output_columns_to_keep]
# order
scenic_output <- scenic_output[order(scenic_output[['triplet_rank']]), ]
# create a new workbook
wb_scenic = createWorkbook()
# create the sheet
scenic_sheet_data = createSheet(wb_scenic, 'SCENICplus_triplets')
# add dataframe of actual data to sheet
addDataFrame(scenic_output, sheet = scenic_sheet_data, startColumn=1, row.names=FALSE)
# set the output location
scenic_excel_output_loc <- paste('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/ST6_scenic_triplets.xlsx', sep = '/')
# write the result
saveWorkbook(wb_scenic, scenic_excel_output_loc)

