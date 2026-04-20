#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_overlap_scenic_qtls_with_gwas.R
# Function: overlap SCENIC+ data with eQTLs/caQTLs
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(stringr)
library(UpSetR)
library(progress)


####################
# Functions        #
####################


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


get_top_effect_per_cs <- function(cs_output_per_ct, feature_column='feature_id', variant_column='snp_id', sort_column='p_value', cs_column='CS', decreasing_sort=F) {
  # store the variant-feature confinement per ct
  var_feature_per_ct <- list()
  # check each cell type
  for (cell_type in names(cs_output_per_ct)) {
    # get output for this cell type
    cs_output_ct <- cs_output_per_ct[[cell_type]]
    # order
    cs_output_ct <- cs_output_ct[order(cs_output_ct[[sort_column]], decreasing = decreasing_sort), ]
    # check if there are entries where the CS is not known
    cs_output_ct_cs_na <- is.na(cs_output_ct[[cs_column]])
    # if there are more than zero
    if (length(cs_output_ct_cs_na) > 0) {
      # set the CS as 'U'
      cs_output_ct[cs_output_ct_cs_na, ][[cs_column]] <- 'U'
      # get the top variant per feature and cs
      cs_output_ct <- cs_output_ct[!duplicated(paste(cs_output_ct[[feature_column]], cs_output_ct[[cs_column]])), ]
      # get the entries where the CS is 'U'
      cs_output_ct_u <- cs_output_ct[cs_output_ct[[cs_column]] == 'U', ]
      # count the credible sets with unknown credible sets
      cs_output_ct_u_ct_counts <- data.frame(table(cs_output_ct_u[[feature_column]]))
      # get the ones that have a credible set in addition to the U one
      cs_output_ct_u_ct_also_not_u <- cs_output_ct_u_ct_counts[cs_output_ct_u_ct_counts[['Freq']] > 1, ]
      # check if there are any
      if (nrow(cs_output_ct_u_ct_also_not_u) > 1) {
        # get the features
        cs_output_ct_u_ct_also_not_u_features <- cs_output_ct_u_ct_also_not_u[[feature_id]]
        # then keep only what is not U and at the same time also not-U
        cs_output_ct <- cs_output_ct[
          !(cs_output_ct[[cs_column]] == 'U' & cs_output_ct[[feature_column]] %in% cs_output_ct_u_ct_also_not_u_features), 
        ]
      }
    }
    # get the variant and feature
    top_var_feature_per_cs <- data.frame('cell_type' = rep(cell_type, times = nrow(cs_output_ct)), 'variant' = cs_output_ct[[variant_column]], 'feature' = cs_output_ct[[feature_column]])
    # put in the list
    var_feature_per_ct[[cell_type]] <- top_var_feature_per_cs
  }
  # merge all
  var_feature_all <- do.call('rbind', var_feature_per_ct)
  return(var_feature_all)
}


add_associated_region <- function(eqtl_outputs, scenic_output, eqtl_gene_column='feature_id', scenic_gene_column='Gene', scenic_region_column='Region') {
  # keep a list per cell type
  cs_output_per_ct_region <- list()
  # check each cell type
  for (cell_type in names(eqtl_outputs)) {
    # extract for cell type
    cs_output_ct <- eqtl_outputs[[cell_type]]
    # subset scenic
    scenic_relevant <- scenic_output[, c(..scenic_gene_column, ..scenic_region_column)]
    # with standardised columns
    colnames(scenic_relevant) <- c('feature', 'region')
    # add associated region to gene
    cs_output_ct <- merge(
      cs_output_ct, 
      scenic_relevant, 
      by.x = eqtl_gene_column, 
      by.y = 'feature'
    )
    # put back in the list
    cs_output_per_ct_region[[cell_type]] <- cs_output_ct
  }
  return(cs_output_per_ct_region)
}


add_all_gwas_traits <- function(variants, ld_proxies, variant_to_trait, ld_proxy_index_column='index_variant', ld_proxy_ld_column='ld_variant', ld_r2_column='R2', variant_to_trait_var_column='gwas_variant', variant_to_trait_trait_column='MAPPED_TRAIT') {
  # we need the number of variants, so we can check each of them
  n_variants <- length(variants)
  # we'll make a dataframe of all variants and the traits we mapped to them
  variant_trait_df <- data.frame('variant' = variants, trait = NA)
  # initialize progress bar
  pb <- progress_bar$new(total = n_variants)
  pb$tick(0)
  # check each variant
  for (i in 1:n_variants) {
    # extract variant
    variant <- variants[i]
    # we'll keep a vector of traits we find
    traits_variant_direct <- NULL
    traits_variant_indirect <- NULL
    # check if it is in the variant to trait mapping
    if (variant %in% variant_to_trait[[variant_to_trait_var_column]]) {
      traits_variant_direct <- variant_to_trait[variant_to_trait[[variant_to_trait_var_column]] == variant, ][[variant_to_trait_trait_column]]
    }
    # now check if we have LD proxies for this variant
    if (variant %in% ld_proxies[[ld_proxy_index_column]]) {
      # get the ld proxies for this variant
      ld_proxies_variant <- ld_proxies[ld_proxies[[ld_proxy_index_column]] == variant, ][[ld_proxy_ld_column]]
      # check each of those for traits
      for (j in 1:length(ld_proxies_variant)) {
        ld_proxy <- ld_proxies_variant[j]
        # check if this proxy is in the variant to trait mapping
        if (ld_proxy %in% variant_to_trait[[variant_to_trait_var_column]]) {
          traits_variant_indirect_proxy <- variant_to_trait[variant_to_trait[[variant_to_trait_var_column]] == ld_proxy, ][[variant_to_trait_trait_column]]
          # add to vector if not empty
          if (!is.null(traits_variant_indirect) & is.null(traits_variant_indirect)) {
            traits_variant_indirect <- traits_variant_indirect_proxy
          } else if (!is.null(traits_variant_indirect_proxy) & !is.null(traits_variant_indirect)) {
            traits_variant_indirect <- c(traits_variant_indirect, traits_variant_indirect_proxy)
          }
        }
      }
    }
    # all traits
    traits_variant <- NULL
    # put together
    if (!is.null(traits_variant_direct) & !is.null(traits_variant_indirect)) {
      traits_variant <- c(traits_variant_direct, traits_variant_indirect)
    } else if (!is.null(traits_variant_direct) & is.null(traits_variant_indirect)) {
      traits_variant <- traits_variant_direct
    } else if (is.null(traits_variant_direct) & !is.null(traits_variant_indirect)) {
      traits_variant <- traits_variant_indirect
    } else {
      # nothing
    }
    # add to table if there is information
    if (!is.null(traits_variant)) {
      # make the traits unique
      traits_variant <- unique(traits_variant)
      # sort
      traits_variant <- traits_variant[order(traits_variant)]
      # make into string
      traits_variant_string <- paste(traits_variant, collapse = ';')
      # put into dataframe
      variant_trait_df[i, 'trait'] <- traits_variant_string
    }
    
    # update progress bar
    pb$tick()
  }
  return(variant_trait_df)
}


###################
# Settings        #
###################


# set location of both tables
tf_ieqtls_annotated_significant_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr_with_replication_significant.tsv.gz'
# load those
tf_ieqtls <- fread(tf_ieqtls_annotated_significant_loc, header = T, sep = '\t')

# gwas data location
gwas_catalogue_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full-chromposrefalt.tsv.gz'
gwas_catalogue <- fread(gwas_catalogue_loc, header = T, sep = '\t')
# keep only what is GWS
gwas_catalogue[['P-VALUE']] <- gsub(',', '.', gwas_catalogue[['P-VALUE']])
gwas_catalogue[['P-VALUE']] <- gsub('E', 'e', gwas_catalogue[['P-VALUE']])
gwas_catalogue[['P-VALUE']] <- as.numeric(gwas_catalogue[['P-VALUE']])
gwas_catalogue <- gwas_catalogue[gwas_catalogue[['P-VALUE']] < 5e10-8, ]
# get the ld data
ld_data_gwas_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/moldpairs/1000G_HC/eurpop/immune-gwas-catalog-download-associations-alt-full-moldpairs_ld_pairs_075.tsv.gz'
ld_data_gwas <- fread(ld_data_gwas_loc, header = T, sep = '\t')
# keep only where one of the ld variants is the GWAS variant
ld_data_gwas <- ld_data_gwas[ld_data_gwas[['ld_variant']] %in% gwas_catalogue[['chromposaltref']], ]
# and where the index variant is one of our QTL variants
# ld_data_gwas <- ld_data_gwas[ld_data_gwas[['index_variant']] %in% tf_ieqtls[['variant']], ]
# add info on whether the variant was in the gwas
tf_ieqtls[['variant_in_imm_gwas']] <- tf_ieqtls[['variant']] %in% gwas_catalogue[['chromposaltref']]
# or had an LD variant
tf_ieqtls[['variant_ld_imm_gwas']] <- tf_ieqtls[['variant']] %in% ld_data_gwas[['ld_variant']]
# and if either
tf_ieqtls[['variant_imm_gwas']] <- tf_ieqtls[['variant_in_imm_gwas']] | tf_ieqtls[['variant_ld_imm_gwas']]
# keep ones that have some gwas signal
tf_ieqtls_wgwas <- tf_ieqtls[tf_ieqtls[['variant_imm_gwas']] == T, ]
# get the outer join of variants we are talking about
tf_ieqtls_gwas_variants_associated <- unique(c(
  tf_ieqtls[tf_ieqtls[['variant_imm_gwas']], ][['variant']],
  ld_data_gwas[ld_data_gwas[['index_variant']] %in% tf_ieqtls[tf_ieqtls[['variant_ld_imm_gwas']] == T, ][['variant']], ][['ld_variant']]
))
# get the traits for those
tf_ieqtls_gwas_variants_associated_traits <- merge(x = data.frame(gwas_variant = tf_ieqtls_gwas_variants_associated), y = gwas_catalogue[, c('chromposaltref', 'MAPPED_TRAIT')], by.x = 'gwas_variant', by.y = 'chromposaltref', all.x = T)
# now use that to get the gwas info for each variant
tf_ieqtls_gwas_variants_to_traits <- add_all_gwas_traits(tf_ieqtls_gwas_variants_associated, ld_data_gwas, tf_ieqtls_gwas_variants_associated_traits)
# and add to the table
tf_ieqtls[['gwas']] <- tf_ieqtls_gwas_variants_to_traits[match(tf_ieqtls[['variant']], tf_ieqtls_gwas_variants_to_traits[['variant']]), ][['trait']]
# write the result
tf_ieqtls_wgwas_loc <- '~/multiome/tables/tf_ieqtls_gwas_overlap.tsv.gz'
write.table(tf_ieqtls, gzfile(tf_ieqtls_wgwas_loc), col.names = T, row.names = F, sep = '\t')
mdfiver::create_sha256_for_file(tf_ieqtls_wgwas_loc)


# location of the SCENIC outout
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both_filtered.tsv.gz'
# read that
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
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


# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# location of the QTL outputs
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/caqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
eqtl_outputs <- get_output_per_celltype_limix(eqtl_output_loc)
# get the top effects per credible set
eqtl_top_var_feature_per_cs <- get_top_effect_per_cs(eqtl_outputs)
# read the eQTL output
caqtl_outputs <- get_output_per_celltype_limix(caqtl_output_loc)
# get the top effects per credible set
caqtl_top_var_feature_per_cs <- get_top_effect_per_cs(caqtl_outputs)


# add scenic region info
eqtl_outputs_scenic_regions <- add_associated_region(eqtl_outputs, scenic_output, scenic_region_column = 'region_cpeaks')
# add cell types
for (ct in names(eqtl_outputs_scenic_regions)) {
  eqtl_outputs_scenic_regions[[ct]] <- cbind(data.frame('cell_type' = rep(ct, times = nrow(eqtl_outputs_scenic_regions[[ct]]))), eqtl_outputs_scenic_regions[[ct]])
}
# then merge
eqtl_outputs_scenic_regions_all <- rbindlist(eqtl_outputs_scenic_regions)
# location of the variant annotations
variant_annotation_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
variant_annotation <- fread(variant_annotation_loc, header = T, sep = '\t')
# add the cpeaks region
eqtl_outputs_scenic_regions_all[['variant_region']] <- variant_annotation[match(eqtl_outputs_scenic_regions_all[['snp_id']], variant_annotation[['snp_id']]), ][['overlapping_feature']]
# then keep all where the associated region is the same as the region where the variant is located
eqtl_outputs_scenic_regions_all_varinregion <- eqtl_outputs_scenic_regions_all[
  !is.na(eqtl_outputs_scenic_regions_all[['variant_region']]) & 
    !is.na(eqtl_outputs_scenic_regions_all[['region']]) & 
    eqtl_outputs_scenic_regions_all[['variant_region']] == eqtl_outputs_scenic_regions_all[['region']], 
  
]
# keep what is significant at variant level as well
eqtl_outputs_scenic_regions_all_varinregion <- eqtl_outputs_scenic_regions_all_varinregion[eqtl_outputs_scenic_regions_all_varinregion[['p_value']] <= eqtl_outputs_scenic_regions_all_varinregion[['pval_nominal_threshold_global']], ]
# get the unique eqtl variants
eqtl_outputs_scenic_regions_all_varinregion_variants <- unique(eqtl_outputs_scenic_regions_all_varinregion[['snp_id']])
# now use that to get the gwas info for each variant
eqtl_outputs_scenic_regions_all_varinregion_variants_to_traits <- add_all_gwas_traits(eqtl_outputs_scenic_regions_all_varinregion_variants, ld_data_gwas, gwas_catalogue[, c('chromposaltref', 'MAPPED_TRAIT')], variant_to_trait_var_column = 'chromposaltref', variant_to_trait_trait_column = 'MAPPED_TRAIT')
# and add to the table
eqtl_outputs_scenic_regions_all_varinregion[['gwas']] <- eqtl_outputs_scenic_regions_all_varinregion_variants_to_traits[match(eqtl_outputs_scenic_regions_all_varinregion[['snp_id']], eqtl_outputs_scenic_regions_all_varinregion_variants_to_traits[['variant']]), ][['trait']]
# add r2g
eqtl_outputs_scenic_regions_all_varinregion[['r2g']] <- paste(eqtl_outputs_scenic_regions_all_varinregion[['region']], eqtl_outputs_scenic_regions_all_varinregion[['feature_id']], sep = '_')
scenic_output[['r2g']] <- paste(scenic_output[['region_cpeaks']], scenic_output[['Gene']], sep = '_')
# add the TF info
eqtl_outputs_scenic_regions_all_varinregion <- merge(eqtl_outputs_scenic_regions_all_varinregion, scenic_output[, c('r2g', 'TF')])
# remove the r2g
eqtl_outputs_scenic_regions_all_varinregion[['r2g']] <- NULL
# store the overlapping data
eqtl_output_scenic_overlap_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/mo_eqtl_overlap_scenic_then_gwas.tsv.gz'
write.table(eqtl_outputs_scenic_regions_all_varinregion, gzfile(eqtl_output_scenic_overlap_loc), col.names = T, row.names = F, sep = '\t')
mdfiver::create_sha256_for_file(eqtl_output_scenic_overlap_loc)

# add cell types
for (ct in names(eqtl_outputs)) {
  eqtl_outputs[[ct]] <- cbind(data.frame('cell_type' = rep(ct, times = nrow(eqtl_outputs[[ct]]))), eqtl_outputs[[ct]])
}
# add cell types
for (ct in names(caqtl_outputs)) {
  caqtl_outputs[[ct]] <- cbind(data.frame('cell_type' = rep(ct, times = nrow(caqtl_outputs[[ct]]))), caqtl_outputs[[ct]])
}
# merge all of them
eqtl_outputs_all <- rbindlist(eqtl_outputs, fill = T)
caqtl_outputs_all <- rbindlist(caqtl_outputs, fill = T)
# remove variants not at threshold
eqtl_outputs_all <- eqtl_outputs_all[eqtl_outputs_all[['p_value']] <= eqtl_outputs_all[['pval_nominal_threshold_global']], ]
# fix significance column
caqtl_outputs_all[is.na(caqtl_outputs_all[['pval_nominal_threshold_global']]), ][['pval_nominal_threshold_global']] <- caqtl_outputs_all[is.na(caqtl_outputs_all[['pval_nominal_threshold_global']]), ][['nominal_p_value_cutoff']]
caqtl_outputs_all <- caqtl_outputs_all[caqtl_outputs_all[['p_value']] <= caqtl_outputs_all[['pval_nominal_threshold_global']], ]
# extract variant id
eqtl_outputs_all_variant_id <- eqtl_outputs_all[['snp_id']]
caqtl_outputs_all_variant_id <- caqtl_outputs_all[['snp_id']]
# and cell type
eqtl_outputs_all_cell_type <- eqtl_outputs_all[['cell_type']]
caqtl_outputs_all_cell_type <- caqtl_outputs_all[['cell_type']]
# rename columns
colnames(eqtl_outputs_all) <- paste('e', colnames(eqtl_outputs_all), sep = '_')
colnames(caqtl_outputs_all) <- paste('ca', colnames(caqtl_outputs_all), sep = '_')
# add variant id back
eqtl_outputs_all <- cbind(data.frame('variant_id' = eqtl_outputs_all_variant_id, 'cell_type' = eqtl_outputs_all_cell_type), eqtl_outputs_all)
caqtl_outputs_all <- cbind(data.frame('variant_id' = caqtl_outputs_all_variant_id, 'cell_type' = caqtl_outputs_all_cell_type), caqtl_outputs_all)
# merge the outputs based on overlapping variant (for that cell type)
qtl_outputs_all_overlapping <- merge(x = caqtl_outputs_all, y = eqtl_outputs_all, by = c('cell_type', 'variant_id'))
# get the unique eqtl variants
dual_qtl_outputs_variants <- unique(qtl_outputs_all_overlapping[['variant_id']])
# now use that to get the gwas info for each variant
dual_qtl_outputs_variants_to_traits <- add_all_gwas_traits(dual_qtl_outputs_variants, ld_data_gwas, gwas_catalogue[, c('chromposaltref', 'MAPPED_TRAIT')], variant_to_trait_var_column = 'chromposaltref', variant_to_trait_trait_column = 'MAPPED_TRAIT')
# and add to the table
qtl_outputs_all_overlapping[['gwas']] <- dual_qtl_outputs_variants_to_traits[match(qtl_outputs_all_overlapping[['variant_id']], dual_qtl_outputs_variants_to_traits[['variant']]), ][['trait']]
# add r2g
qtl_outputs_all_overlapping[['r2g']] <- paste(qtl_outputs_all_overlapping[['ca_feature_id']], qtl_outputs_all_overlapping[['e_feature_id']], sep = '_')
# add TF (where possible)
qtl_outputs_all_overlapping <- merge(qtl_outputs_all_overlapping, scenic_output[, c('r2g', 'TF')], all.x = T)
# store the overlapping data
dual_qtl_outputs_variants_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/mo_dual_qtl_overlap_then_gwas.tsv.gz'
write.table(qtl_outputs_all_overlapping, gzfile(dual_qtl_outputs_variants_loc), col.names = T, row.names = F, sep = '\t')
mdfiver::create_sha256_for_file(dual_qtl_outputs_variants_loc)
