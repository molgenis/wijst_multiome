#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_overlap_stim_cre_and_tf_ieqtls.R
# Function: merge the QTL results across cell types into an Excels
############################################################################################################################

####################
# libraries        #
####################

library(xlsx)
library(data.table)
library(stringr)


###################
# Functions        #
####################

get_qtls_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all_qval_allchroms_fdr005_significant.txt.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column='pval_nominal_threshold_global', nominal_significance_column='p_value', verbose=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  qtls_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(qtl_output_loc, '/', cell_type, '/', output_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    # filter the results on significance
    if (!is.null(significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by significance column', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[significance_column]]) &
          cell_type_output[[significance_column]] < significance_cutoff, 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by significance column', nrow(cell_type_output)))
      }
    }
    if (!is.null(nominal_cutoff_column) & !is.null(nominal_significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering by nominal cutoff', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[nominal_cutoff_column]]) & 
          !is.na(cell_type_output[[nominal_significance_column]]) &
          cell_type_output[[nominal_significance_column]] <= cell_type_output[[nominal_cutoff_column]], 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering by nominal cutoff', nrow(cell_type_output)))
      }
    }
    # add the celltype as a column
    cell_type_output[['cell_type']] <- cell_type
    # add to the list
    qtls_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(qtls_per_celltype)
}


#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @param gene_column which column to use as the gene identifier
#' @returns a list with the output tables per cell type
#' 
get_output_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all_qval_allchroms_fdr005_significant_cs.tsv.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, verbose=T, add_nominal_cutoff=T, nominal_p_value_column='p_value') {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  egenes_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(qtl_output_loc, '/', cell_type, '/', output_file, sep = '')
    # log if requested
    if (verbose) {
      print(paste('reading', full_cell_type_path))
    }
    # read the file
    cell_type_output <- read.table(full_cell_type_path, sep = '\t', header = T)
    # filter the results on significance
    if (!is.null(significance_column)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
      }
      # filter
      cell_type_output <- cell_type_output[
        !is.na(cell_type_output[[significance_column]]) &
          cell_type_output[[significance_column]] < significance_cutoff, 
      ]
      if (verbose) {
        print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
      }
      # add nominal cutoff by first getting the top effects
      cell_type_output_top <- cell_type_output[order(cell_type_output[[nominal_p_value_column]]), ]
      cell_type_output_top <- cell_type_output_top[!duplicated(cell_type_output_top[[gene_column]]), ]
      # get the highest still significant p value
      max_sig_p <- max(cell_type_output_top[[nominal_p_value_column]])
      # and add that to the output
      cell_type_output[['nominal_p_value_cutoff']] <- max_sig_p
    }
    # add to the list
    egenes_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(egenes_per_celltype)
}



add_top_effect_annotation <- function(qtls_per_celltype, feature_column='feature_id', variant_column='snp_id', significance_column='p_value', decreasing=F) {
  # go through each of the cell types
  for (cell_type in names(qtls_per_celltype)) {
    # extract that table
    qtls_celltype <- qtls_per_celltype[[cell_type]]
    # order by the significancce
    qtls_celltype_ordered <- qtls_celltype[order(qtls_celltype[[significance_column]], decreasing = decreasing), ]
    # keep only the top effect
    qtls_celltype_top <- qtls_celltype_ordered[!duplicated(qtls_celltype_ordered[[feature_column]]), ]
    # and annotate in the original table if the variant-feature combination was the top one
    qtls_celltype[['is_top_variant']] <- paste(qtls_celltype[[variant_column]], qtls_celltype[[feature_column]]) %in% paste(qtls_celltype_top[[variant_column]], qtls_celltype_top[[feature_column]])
    # put that back in the list
    qtls_per_celltype[[cell_type]] <- qtls_celltype
  }
  return(qtls_per_celltype)
}


get_closest_flanks <- function(position_table, left_flank_column1, right_flank_column1, left_flank_column2, right_flank_column2) {
  # get the distance between left flanks
  dist_left_flank1_to_left_flank2 <- position_table[[left_flank_column1]] - position_table[[left_flank_column2]]
  # distance between the right flanks
  dist_right_flank1_to_right_flank2 <- position_table[[right_flank_column1]] - position_table[[right_flank_column2]]
  # distance between left flank 1 and right flank 2
  dist_left_flank1_to_right_flank2 <- position_table[[left_flank_column1]] - position_table[[right_flank_column2]]
  # distance between right flank 1 and left flank 2
  dist_right_flank1_to_left_flank2 <- position_table[[right_flank_column1]] - position_table[[left_flank_column2]]
  # put in a table for convenience sake
  distances_tbl <- data.table(
    'lf1_to_lf2' = dist_left_flank1_to_left_flank2, 
    'rf1_to_rf2' = dist_right_flank1_to_right_flank2, 
    'lf1_to_rf2' = dist_left_flank1_to_right_flank2, 
    'rf1_to_lf2' = dist_right_flank1_to_left_flank2
  )
  # add the minimum absolute distance
  distances_tbl[['min_dist']] <- apply(distances_tbl, 1, function(x) {
    return(min(abs(x)))
  })
  # but set this to zero if any of the flanks end in the bodies
  #              -----
  #                 ++++
  distances_tbl[(distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['lf1_to_rf2']] > 0) |
                  #                   ----
                #                 ++++
                (distances_tbl[['rf1_to_lf2']] > 0 & distances_tbl[['lf1_to_rf2']] < 0) |
                  #                   ----
                #                 +++++++++
                (distances_tbl[['lf1_to_lf2']] > 0 & distances_tbl[['rf1_to_rf2']] < 0) |
                  #                 ---------
                #                   ++++
                (distances_tbl[['lf1_to_lf2']] < 0 & distances_tbl[['rf1_to_rf2']] > 0)
                , 'min_dist'] <- 0
  return(distances_tbl)
}


add_significant_celltypes_as_strings_vectorised <- function(overlap_table, variant1_col='hit1', variant2_col='hit2', trait1_col='trait1', trait2_col='trait2', cell_type_column='cell_type') {
  # convert to dataframe format
  overlap_table_df <- data.frame(overlap_table)
  # add overlap 
  overlap_table_df[['full_overlap']] <- paste(overlap_table_df[[variant1_col]], overlap_table_df[[variant2_col]],
                                              overlap_table_df[[trait1_col]], overlap_table_df[[trait2_col]])
  
  # split into groups
  split_cts <- split(overlap_table_df[[cell_type_column]], overlap_table_df[['full_overlap']])
  
  # summarise to string per group
  ct_strings <- sapply(split_cts, function(x) {
    ct_string <- paste(sort(unique(x)), collapse = ",")
    return(ct_string)
  })
  
  # map back
  overlap_table_df[['other_ct']] <- ct_strings[overlap_table_df[['full_overlap']]]
  # back to data table format
  return(data.table(overlap_table_df))
}



####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- F


#######################################
# Read TF-i-eQTL original/replication #
#######################################

# location of the tf-i-ieqtls
tf_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz'
# location of the region-i-ieqtls
region_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/region_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr.tsv.gz'
# read the inputs
tf_ieqtls <- fread(tf_ieqtl_loc, header = T, sep = '\t')
region_ieqtls <- fread(region_ieqtl_loc, header = T, sep = '\t')

# get the replication set
tf_ieqtl_rep_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion_significant/merged/results_fdr.tsv.gz'
tf_ieqtls_rep <- fread(tf_ieqtl_rep_loc, header = T, sep = '\t')

# remove the TF-i-eGene that was removed from the name (it was added in the second pass)
tf_ieqtls_rep[['region']] <- gsub('\\)_.*', ')', tf_ieqtls_rep[['region']])

# rename region to TF where applicable
colnames(tf_ieqtls) <- gsub('region', 'tf', colnames(tf_ieqtls))
colnames(tf_ieqtls_rep) <- gsub('region', 'tf', colnames(tf_ieqtls_rep))


############################
# Read SCENIC+ information #
############################

# the location of SCENIC output
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# remove what SCENIC thinks is less likely
scenic_output <- scenic_output[scenic_output[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# order by the extended
scenic_output <- scenic_output[order(scenic_output[['is_extended']]), ]
# get which are non-extended if it was both extended and non-extended
scenic_eregs <- unique(scenic_output[, c('TF', 'Gene_signature_direction', 'Gene_signature_name', 'source')])
scenic_eregs <-scenic_eregs[order(scenic_eregs[['source']]), ]
scenic_eregs_to_keep <- scenic_eregs[!duplicated(paste(scenic_eregs[['TF']], scenic_eregs[['Gene_signature_direction']])), ]
# then use that to filer
scenic_output <- scenic_output[scenic_output[['Gene_signature_name']] %in% scenic_eregs_to_keep[['Gene_signature_name']], ]
# rename the regions
scenic_output[['region_cpeaks']] <- gsub(':', '-', scenic_output[['Region']])


############################
# add annotations for gene #
############################

# read the location of the UCSC annotations
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/genome_annotation.tsv'
gene_anno <- fread(gene_anno_loc, header = T, sep = '\t')
# rename columns to be the same as in limix
colnames(gene_anno) <-c('chrom', 'start', 'end', 'strand', 'gs','Transcription_Start_Site','Transcript_type')
# get gene info
tf_ieqtls_gene_info <- gene_anno[match(tf_ieqtls[['gene']], gene_anno[['gs']]), ]
# make columns gene specific
colnames(tf_ieqtls_gene_info) <- c('gene_chromosome', 'gene_start', 'gene_end', 'gene_strand', 'gs', 'gene_tss', 'gene_tt')
# remove gs column
tf_ieqtls_gene_info[['gs']] <- NULL
# add to the overlap table
tf_ieqtls <- cbind(tf_ieqtls, tf_ieqtls_gene_info)

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = '\t')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])
# get cpeak info
tf_ieqtls_peak_info <- cpeaks_anno[match(tf_ieqtls[['region']], cpeaks_anno[['signac_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')]
# make columns gene specific
colnames(tf_ieqtls_peak_info) <- c('peak_chromosome', 'peak_start', 'peak_end')
# add to the overlap table
tf_ieqtls <- cbind(tf_ieqtls, tf_ieqtls_peak_info)

# read the variant annotation
variant_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
variant_anno <- fread(variant_anno_loc, header = T, sep = '\t')
# get info
tf_ieqtls_variant_info <- variant_anno[match(tf_ieqtls[['variant']], variant_anno[['snp_id']]), c('snp_chromosome', 'snp_position')]
# add to the overlap table
tf_ieqtls <- cbind(tf_ieqtls, tf_ieqtls_variant_info)

# same for variant to gene
var_gene_distances <- get_closest_flanks(tf_ieqtls, 'snp_position', 'snp_position', 'gene_start', 'gene_end')
tf_ieqtls[['var_gene_dist']] <- tf_ieqtls[['min_dist']]


########################
# Read QTL information #
########################

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
eqtl_outputs <- get_output_per_celltype_limix(eqtl_output_loc)
# location of the QTL outputs
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/caqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
caqtl_outputs <- get_output_per_celltype_limix(caqtl_output_loc)

# add cell types in eQTLs
for (cell_type in names(eqtl_outputs)) {
  eqtl_outputs[[cell_type]][['cell_type']] <- cell_type
}
# add cell types in caQTLs
for (cell_type in names(caqtl_outputs)) {
  caqtl_outputs[[cell_type]][['cell_type']] <- cell_type
}
# merge results
eqtl_outputs_all <- rbindlist(eqtl_outputs, fill = T)
caqtl_outputs_all <- rbindlist(caqtl_outputs, fill = T)
# in the case we didn't name the column correctly, fix that
caqtl_outputs_all[is.na(caqtl_outputs_all[['pval_nominal_threshold_global']]), ][['pval_nominal_threshold_global']] <- caqtl_outputs_all[is.na(caqtl_outputs_all[['pval_nominal_threshold_global']]), ][['nominal_p_value_cutoff']]
# filter on significance
eqtl_outputs_all_sig <- eqtl_outputs_all[eqtl_outputs_all[['feature_q_value']] < 0.05 &
                             eqtl_outputs_all[['p_value']] <= eqtl_outputs_all[['pval_nominal_threshold_global']], ]
caqtl_outputs_all_sig <- caqtl_outputs_all[caqtl_outputs_all[['feature_q_value']] < 0.05 &
                                            caqtl_outputs_all[['p_value']] <= caqtl_outputs_all[['pval_nominal_threshold_global']], ]


# location of the i-eqtl output
ieqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# get all the QTL output
ieqtl_output <- get_qtls_per_celltype_limix(ieqtl_output_loc, output_file = 'inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', gene_column='feature', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column=NULL, nominal_significance_column=NULL)
# merge all the results
ieqtl_output_all <- do.call('rbind', ieqtl_output)
# filter on signifiacnce
ieqtl_output_all_sig <- ieqtl_output_all[ieqtl_output_all[['feature_q_value']] < 0.05 &
                                           ieqtl_output_all[['feature_bf_eigen']] < 0.05, ]
# location of the i-eqtl output
icaqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# get all the QTL output
icaqtl_output <- get_qtls_per_celltype_limix(icaqtl_output_loc, output_file = 'inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', gene_column='feature', significance_column='feature_q_value', significance_cutoff=0.05, nominal_cutoff_column=NULL, nominal_significance_column=NULL)
# merge all the results
icaqtl_output_all <- do.call('rbind', icaqtl_output)
# filter on signifiacnce
icaqtl_output_all_sig <- icaqtl_output_all[icaqtl_output_all[['feature_q_value']] < 0.05 &
                                            icaqtl_output_all[['feature_bf_eigen']] < 0.05, ]

# add other cell types for each QTL
eqtl_outputs_all_sig <- add_significant_celltypes_as_strings_vectorised(eqtl_outputs_all_sig, variant1_col = 'snp_id', variant2_col = 'snp_id', trait1_col = 'feature_id', trait2_col = 'feature_id', cell_type_column = 'cell_type')
caqtl_outputs_all_sig <- add_significant_celltypes_as_strings_vectorised(caqtl_outputs_all_sig, variant1_col = 'snp_id', variant2_col = 'snp_id', trait1_col = 'feature_id', trait2_col = 'feature_id', cell_type_column = 'cell_type')
ieqtl_output_all_sig <- add_significant_celltypes_as_strings_vectorised(ieqtl_output_all_sig, variant1_col = 'snp_id', variant2_col = 'snp_id', trait1_col = 'feature_id', trait2_col = 'feature_id', cell_type_column = 'cell_type')
icaqtl_output_all_sig <- add_significant_celltypes_as_strings_vectorised(icaqtl_output_all_sig, variant1_col = 'snp_id', variant2_col = 'snp_id', trait1_col = 'feature_id', trait2_col = 'feature_id', cell_type_column = 'cell_type')

# add eQTL celltype infos. The other_ct column already includes the eQTL cell type itself, so we just need the first match
tf_ieqtls[['is_eqtl']] <- eqtl_outputs_all_sig[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['gene']]), paste(eqtl_outputs_all_sig[['snp_id']], eqtl_outputs_all_sig[['feature_id']])), ][['other_ct']]
tf_ieqtls[['is_ieqtl']] <- ieqtl_output_all_sig[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['gene']]), paste(ieqtl_output_all_sig[['snp_id']], ieqtl_output_all_sig[['feature']])), ][['other_ct']]
# add caQTL variant info, we don't have the region, so we just use the variant information
tf_ieqtls[['is_caqtl_variant']] <- caqtl_outputs_all_sig[match(paste(tf_ieqtls[['variant']]), paste(caqtl_outputs_all_sig[['snp_id']])), ][['other_ct']]
tf_ieqtls[['is_icaqtl_variant']] <- icaqtl_output_all[match(paste(tf_ieqtls[['variant']]), paste(icaqtl_output_all[['snp_id']])), ][['other_ct']]


########################
# add replication info #
########################

# and add the second round info
tf_ieqtls[['tf:genotype_p_rep']] <- tf_ieqtls_rep[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']], tf_ieqtls[['gene']]), paste(tf_ieqtls_rep[['variant']], tf_ieqtls_rep[['tf']], tf_ieqtls_rep[['gene']])), ][['tf:genotype_p']]
tf_ieqtls[['tf:genotype_bh_rep']] <- tf_ieqtls_rep[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']], tf_ieqtls[['gene']]), paste(tf_ieqtls_rep[['variant']], tf_ieqtls_rep[['tf']], tf_ieqtls_rep[['gene']])), ][['tf:genotype_bh']]
tf_ieqtls[['anova_p_rep']] <- tf_ieqtls_rep[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']], tf_ieqtls[['gene']]), paste(tf_ieqtls_rep[['variant']], tf_ieqtls_rep[['tf']], tf_ieqtls_rep[['gene']])), ][['anova_p']]
tf_ieqtls[['anova_bh_rep']] <- tf_ieqtls_rep[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']], tf_ieqtls[['gene']]), paste(tf_ieqtls_rep[['variant']], tf_ieqtls_rep[['tf']], tf_ieqtls_rep[['gene']])), ][['anova_bh']]
tf_ieqtls[['tf:genotype_beta_rep']] <- tf_ieqtls_rep[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']], tf_ieqtls[['gene']]), paste(tf_ieqtls_rep[['variant']], tf_ieqtls_rep[['tf']], tf_ieqtls_rep[['gene']])), ][['tf:genotype_beta']]
tf_ieqtls[['tf:genotype_se_rep']] <- tf_ieqtls_rep[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']], tf_ieqtls[['gene']]), paste(tf_ieqtls_rep[['variant']], tf_ieqtls_rep[['tf']], tf_ieqtls_rep[['gene']])), ][['tf:genotype_se']]
tf_ieqtls[['tf_bh_rep']] <- tf_ieqtls_rep[match(paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']], tf_ieqtls[['gene']]), paste(tf_ieqtls_rep[['variant']], tf_ieqtls_rep[['tf']], tf_ieqtls_rep[['gene']])), ][['tf_bh']]

# add info on the direction of the TF
tf_ieqtls[['tf_gene_direction']] <- str_extract(tf_ieqtls[['tf']], '\\+\\/\\+|\\-\\/\\-|\\+\\/\\-|\\-\\/\\+')
# and whether it was flipped for the replication
tf_ieqtls[['replication_flip']] <- tf_ieqtls[['tf_gene_direction']] == '-/+'


#######################
# add CRE-i-eQTL info #
#######################

# get significant CRE-i-eQTLs
region_ieqtls_sig <- region_ieqtls[
  !is.na(region_ieqtls[['region:genotype_bh']]) &
    region_ieqtls[['region:genotype_bh']] < 0.05 & 
    region_ieqtls[['anova_bh']] < 0.05 & 
    # region_ieqtls[['region_bh']] < 0.05 &
    # region_ieqtls[['genotype_bh']] < 0.05,
    region_ieqtls[['region_bh']] < 0.05,
]
# add whether the variant-gene pair is also a significant CRE-i-eQTL
tf_ieqtls[['is_cre_ieqtl']] <- paste(tf_ieqtls[['variant']], tf_ieqtls[['gene']]) %in% paste(region_ieqtls_sig[['variant']], region_ieqtls_sig[['gene']])


##########################
# Read TFQTL information #
##########################

# location of TFQTL output
# tfqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/tfqtl/limix_sc/output/tf/all/merged/results_fdr.tsv.gz'
tfqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/tfqtl/limix_sc/output/tf_noiegene/all/merged/results_fdr.tsv.gz'
# read table
tfqtl_output <- fread(tfqtl_output_loc, header = T, sep = '\t')
# subset to what is significant
tfqtl_output_sig <- tfqtl_output[tfqtl_output[['genotype_bf']] < 0.05, ]
# remove the gene that was removed
tfqtl_output_sig[['tf']]<- gsub('\\)_.*', ')', tfqtl_output_sig[['tf']])
# add whether the variant-gene pair is also a significant CRE-i-eQTL
tf_ieqtls[['is_tfqtl']] <- paste(tf_ieqtls[['variant']], tf_ieqtls[['tf']]) %in% paste(tfqtl_output_sig[['variant']], tfqtl_output_sig[['tf']])


#########################
# get significance info #
#########################

# add significance
tf_ieqtls[['significant']] <- 
  (tf_ieqtls[['tf:genotype_bh']] < 0.05 & 
  tf_ieqtls[['anova_bh']] < 0.05 & 
  tf_ieqtls[['tf_bh']] < 0.05 & 
  !is.na(tf_ieqtls[['tf:genotype_bh_rep']] < 0.05) & 
  tf_ieqtls[['tf:genotype_bh_rep']] < 0.05 & 
  !is.na(tf_ieqtls[['tf:genotype_bh_rep']]) &
  tf_ieqtls[['anova_bh_rep']] < 0.05 & 
  !is.na(tf_ieqtls[['anova_bh_rep']]) &
  tf_ieqtls[['tf_bh_rep']] < 0.05)


################
# save results #
################

# set location of both tables
tf_ieqtls_annotated_loc <- paste(dirname(tf_ieqtl_loc), 'results_fdr_with_replication.tsv.gz', sep = '/')
tf_ieqtls_annotated_significant_loc <- paste(dirname(tf_ieqtl_loc), 'results_fdr_with_replication_significant.tsv.gz', sep = '/')
# save these
write.table(tf_ieqtls, gzfile(tf_ieqtls_annotated_loc), col.names = T, sep = '\t', quote = F, row.names = F)
mdfiver::create_sha256_for_file(tf_ieqtls_annotated_loc)
write.table(tf_ieqtls[tf_ieqtls[['significant']], ], gzfile(tf_ieqtls_annotated_significant_loc), col.names = T, sep = '\t', quote = F, row.names = F)
mdfiver::create_sha256_for_file(tf_ieqtls_annotated_significant_loc)


#####################
# check enrichments #
#####################

# relative enrichment of TF-i-eQTLs in stimulation-i-eQTLs
tf_ieqtls[['is_any_ieqtl']] <- !is.na(tf_ieqtls[['is_ieqtl']])
table(tf_ieqtls[, c('significant', 'is_any_ieqtl')])
fisher.test(table(tf_ieqtls[, c('significant', 'is_any_ieqtl')]))
# 
# Fisher's Exact Test for Count Data
# 
# data:  table(tf_ieqtls[, c("significant", "is_any_ieqtl")])
# p-value = 6.141e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.014701 3.111693
# sample estimates:
# odds ratio 
#    2.50783 

# location of eQTLgen cis-trans pairs
cistrans_eqtlgen_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/replication/bulkeqtl_positives_annotated.csv.gz'
# read that file
cistrans_eqtlgen <- fread(cistrans_eqtlgen_loc, header = T, sep = ',')
# read the location of the genes
gene_anno_arc_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno_arc <- fread(gene_anno_arc_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno_arc) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')
# add gene symbols to the genes in scenic
tf_ieqtls[['gene_ens']] <- gene_anno_arc[match(tf_ieqtls[['gene']], gene_anno_arc[['gs']]), ][['ens']]
# remove the annotations behind the TF to just get the TF
tf_ieqtls[['tf_short']] <- gsub('_.*', '', tf_ieqtls[['tf']])
# add that gene symbol as well
tf_ieqtls[['tf_ens']] <- gene_anno_arc[match(tf_ieqtls[['tf_short']], gene_anno_arc[['gs']]), ][['ens']]
# add info on whether this is a cis-trans TF
tf_ieqtls[['tf2g_cistrans_eqtlgen']] <- paste(tf_ieqtls[['tf_ens']], tf_ieqtls[['gene_ens']]) %in% paste(cistrans_eqtlgen[['CisGene']], cistrans_eqtlgen[['TransGene']])
tf_ieqtls[['g2tf_cistrans_eqtlgen']] <- paste(tf_ieqtls[['tf_ens']], tf_ieqtls[['gene_ens']]) %in% paste(cistrans_eqtlgen[['TransGene']], cistrans_eqtlgen[['CisGene']])

