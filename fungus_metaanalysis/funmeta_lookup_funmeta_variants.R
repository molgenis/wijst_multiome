#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: funmeta_lookup_funmeta_variants.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(qvalue)

####################
# Functions        #
####################


add_qvalue <- function(qtl_table, p_value_column='empirical_feature_p_value', feature_id_column='feature_id') {
  # sort by significance
  qtl_table <- qtl_table[order(qtl_table[[p_value_column]]), ]
  # now only get the first entry
  qtl_table_unique <- qtl_table[!duplicated(qtl_table[[feature_id_column]]), ]
  # extract the features
  features <- qtl_table_unique[[feature_id_column]]
  # and p values
  p_values <- qtl_table_unique[[p_value_column]]
  # get q values for those
  q_values <- qvalue(p_values)$qvalues
  # make into new table
  q_value_table <- data.table('features' = features, 'qval' = q_values)
  colnames(q_value_table) <- c(feature_id_column, 'qval')
  return(q_value_table)
}


get_qtls_from_files <- function(qtls_loc, variants, qtl_prepend='', qtl_append='.tsv.gz', variant_column='snp_id', add_qvalue=T, p_value_column='empirical_feature_p_value', feature_id_column='feature_id', alpha_param_clean=T, alpha_param_column='alpha_param') {
  # we'll store a results per cell type
  qtls_celltype <- list()
  # list all the files
  qtl_files <- list.files(qtls_loc, full.names = F, recursive = F, include.dirs = F)
  # paste expression together
  file_pattern <- paste0(qtl_prepend, '.*', qtl_append)
  # and filter those files
  qtl_files <- qtl_files[grep(pattern = file_pattern, x = qtl_files)]
  # read each file
  for (qtl_file in qtl_files) {
    # get the full path
    qtl_file_full <- paste0(qtls_loc, '/', qtl_file)
    # read the file
    qtls <- fread(qtl_file_full, header = T, sep = '\t')
    # filter to only include the variants we care about
    if (add_qvalue) {
      # remove entries with weird alpha params
      if(alpha_param_clean) {
        qtls <- qtls[!(qtls[[alpha_param_column]] < .2 | qtls[[alpha_param_column]] > 5), ]
      }
      # get qvalues
      qvalue_table <- add_qvalue(qtl_table = qtls, p_value_column = p_value_column, feature_id_column = feature_id_column)
      # add the tables together
      qtls <- merge(qtls, qvalue_table, by = feature_id_column, all = T)
    }
    # filter on only the variants we care about
    qtls <- qtls[qtls[[variant_column]] %in% variants, ]
    # get the cell type from the name
    cell_type <- gsub(qtl_prepend, '', qtl_file)
    cell_type <- gsub(qtl_append, '', cell_type)
    # add cell type as column
    qtls <- cbind(data.table('cell_type' = rep(cell_type, times = nrow(qtls))), qtls)
    # put in the list
    qtls_celltype[[cell_type]] <- qtls
  }
  # merge all together
  qtls_all <- do.call('rbind', qtls_celltype)
  return(qtls_all)
}

####################
# Main Code        #
####################

# we will look up the variants, and their proxy variants
variant1_ld_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/proxy_variants/proxy_rs3807184.txt'
variant2_ld_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/proxy_variants/proxy_rs413524.txt'
variant3_ld_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/funmeta/proxy_variants/proxy_rs113472423.txt'

# here are the sc-eQTLgen results
sceqtlgen_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta/Out_202406/'

# with a cutoff of 0.8
ld_r2_cutoff <- 0.8

# get the tables
variant1_ld <- fread(variant1_ld_loc, header = T, sep = '\t')
variant2_ld <- fread(variant2_ld_loc, header = T, sep = '\t')
variant3_ld <- fread(variant3_ld_loc, header = T, sep = '\t')

# let's filter first
variant1_ld <- variant1_ld[variant1_ld[['R2']] >= ld_r2_cutoff, ]
variant2_ld <- variant2_ld[variant2_ld[['R2']] >= ld_r2_cutoff, ]
variant3_ld <- variant3_ld[variant3_ld[['R2']] >= ld_r2_cutoff, ]

# add the original variant
variant1_ld <- cbind(data.table('original_snp' = rep('rs3807184', times = nrow(variant1_ld))), variant1_ld)
variant2_ld <- cbind(data.table('original_snp' = rep('rs413524', times = nrow(variant2_ld))), variant2_ld)
variant3_ld <- cbind(data.table('original_snp' = rep('rs113472423', times = nrow(variant3_ld))), variant3_ld)

# we'll keep a set of columns
ld_columns_to_keep <- c('original_snp', 'RS_Number', 'Coord', 'Alleles', 'Distance', 'R2', 'Correlated_Alleles')
# and subset to those
variant1_ld <- variant1_ld[, ..ld_columns_to_keep]
variant2_ld <- variant2_ld[, ..ld_columns_to_keep]
variant3_ld <- variant3_ld[, ..ld_columns_to_keep]
# and we'll give them better names
ld_columns_to_rename <- c('original_snp', 'proxy_snp', 'coordinate', 'alleles', 'distance', 'R2', 'cor_alleles')
# which we'll use
colnames(variant1_ld) <- ld_columns_to_rename
colnames(variant2_ld) <- ld_columns_to_rename
colnames(variant3_ld) <- ld_columns_to_rename

# get qtls for all variants
sceqtlgen_eqtls <- get_qtls_from_files(sceqtlgen_loc, unique(c(variant1_ld[['proxy_snp']], variant2_ld[['proxy_snp']], variant3_ld[['proxy_snp']])))

# there is no need to keep all columns, so let's reduce a bit
#qtl_columns_to_deep <- 
