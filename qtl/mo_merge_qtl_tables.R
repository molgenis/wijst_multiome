#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_merge_qtl_tables.R
# Function: merge the QTL results across cell types into an Excels
############################################################################################################################

####################
# libraries        #
####################

library(xlsx)


####################
# Functions        #
####################


get_label_dict <- function(){
  label_dict <- list()
  # condition combinations
  label_dict[['UT24hCA']] <- 'UT-24hCA'
  # conditions
  label_dict[['UT']] <- 'C'
  label_dict[['C']] <- 'C'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't6-8w'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  # do the datasets
  label_dict[["mo"]] <- "multiome"
  label_dict[["1m"]] <- "NC 2022"
  label_dict[["1M"]] <- "NC 2022"
  return(label_dict)
}


rename_labels <- function(vector_to_rename) {
  # get the labels that are present
  label_renaming <- get_label_dict()
  # convert vector
  vector_to_rename <- as.character(vector_to_rename)
  # now check which labels we are missing
  missing_renames <- setdiff(unique(vector_to_rename), names(label_renaming))
  # add those renames as not being renames
  for (missing_rename in missing_renames) {
    label_renaming[[missing_rename]] <- missing_rename
  }
  # now replace each value with the rename
  renamed_vector <- as.vector(unlist(label_renaming[vector_to_rename]))
  # and return that
  return(renamed_vector)
}


write_output_to_excel <- function(results_per_celltype, excel_output_loc, cell_type_column='cell_type', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, remove_columns=c(), filter_list=NULL) {
  # get unique cell types
  cell_types_tbl <- unique(results_per_celltype[[cell_type_column]])
  # filter on which we want
  if (!is.null(cell_types)) {
    cell_types_tbl <- intersect(cell_types_tbl, cell_types)
  }
  # create a new workbook
  wb = createWorkbook()
  # check each cell type
  for (celltype in cell_types_tbl) {
    # create the name for the sheet
    sheet_name <- celltype
    # replace with a better label if requested
    if (use_label_dict) {
      sheet_name <- rename_labels(celltype)
    }
    # create the sheet
    sheet = createSheet(wb, sheet_name)
    # fetch the dataframe
    celltype_result <- results_per_celltype[!is.na(results_per_celltype[[cell_type_column]]) & results_per_celltype[[cell_type_column]] == celltype, ]
    # again, replace labels if requested
    if (use_label_dict) {
      celltype_result[['cell_type']] <- as.vector(unlist(get_label_dict()[celltype_result[['cell_type']]]))
    }
    if (!is.null(remove_columns)) {
      celltype_result <- celltype_result[, setdiff(colnames(celltype_result), remove_columns)]
    }
    # filter if requested
    if (!is.null(filter_list)) {
      for (filter_column in names(filter_list)) {
        celltype_result <- celltype_result[
          celltype_result[[filter_column]] < filter_list[[filter_column]], 
        ]
      }
    }
    # add dataframe to sheet
    addDataFrame(celltype_result, sheet = sheet, startColumn=1, row.names=FALSE)
  }
  # write the result
  saveWorkbook(wb, excel_output_loc)
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


add_credible_sets_output <- function(output_table_per_celltype, feature_column='feature_id', cs_column='CS', item_column='cell_type') {
  # filter on pval_nominal_threshold_global
  eqtl_outputs_filtered <- list()
  for (ct in names(output_table_per_celltype)) {
    # extract this eQTL output
    ct_eqtl_output <- output_table_per_celltype[[ct]]
    # # filter on the pval threshold
    # ct_eqtl_output <- ct_eqtl_output[ct_eqtl_output[['p_value']] < ct_eqtl_output[['pval_nominal_threshold_global']], ]
    # get the unique combinations of features and credible sets
    ct_eqtl_feature_cs <- unique(ct_eqtl_output[, c(feature_column, cs_column)])
    # get features with NA values
    ct_eqtl_feature_cs_na <- ct_eqtl_feature_cs[is.na(ct_eqtl_feature_cs[[cs_column]]), feature_column]
    # and the ones without
    ct_eqtl_feature_cs_nona <- ct_eqtl_feature_cs[!is.na(ct_eqtl_feature_cs[[cs_column]]), feature_column]
    # get the ones which are only with na
    ct_eqtl_feature_cs_onlyna <- setdiff(ct_eqtl_feature_cs_na, ct_eqtl_feature_cs_nona)
    if (length(ct_eqtl_feature_cs_onlyna)) {
      # for those, set the CS to simply be L1
      ct_eqtl_output[ct_eqtl_output[[feature_column]] %in% ct_eqtl_feature_cs_onlyna, ][[cs_column]] <- 'L1'
    }
    # get the ones which are partly na
    ct_eqtl_feature_cs_somena <- intersect(ct_eqtl_feature_cs_nona, ct_eqtl_feature_cs_na)
    if (length(ct_eqtl_feature_cs_somena) > 0) {
      ct_eqtl_output[(ct_eqtl_output[[feature_column]] %in% ct_eqtl_feature_cs_somena) & is.na(ct_eqtl_output[[cs_column]]), ][[cs_column]] <- 'L0'
    }
    # add ct
    ct_eqtl_output[[item_column]] <- ct
    # put back in list
    eqtl_outputs_filtered[[ct]] <- ct_eqtl_output
  }
  return(eqtl_outputs_filtered)
}

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

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/finemapping/eqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
eqtl_outputs <- get_output_per_celltype_limix(eqtl_output_loc)
# location of the QTL outputs
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/finemapping/caqtl/sc-eqtlgen/combined_with_qtl/combined/L1/'
# read the eQTL output
caqtl_outputs <- get_output_per_celltype_limix(caqtl_output_loc)

# add cs info
eqtl_outputs_cs <- add_credible_sets_output(eqtl_outputs)
caqtl_outputs_cs <- add_credible_sets_output(caqtl_outputs)

# merge all of them
eqtl_outputs_cs_all <- rbindlist(eqtl_outputs_cs, fill = T)
caqtl_outputs_cs_all <- rbindlist(caqtl_outputs_cs, fill = T)
# remove variants not at threshold
eqtl_outputs_cs_all <- eqtl_outputs_cs_all[eqtl_outputs_cs_all[['p_value']] <= eqtl_outputs_cs_all[['pval_nominal_threshold_global']], ]
caqtl_outputs_cs_all <- caqtl_outputs_cs_all[caqtl_outputs_cs_all[['p_value']] <= caqtl_outputs_cs_all[['pval_nominal_threshold_global']], ]

# remove confusing column
eqtl_outputs_cs_all[['pval_nominal_threshold_local']] <- NULL
caqtl_outputs_cs_all[['pval_nominal_threshold_local']] <- NULL
eqtl_outputs_cs_all[['nominal_p_value_cutoff']] <- NULL
caqtl_outputs_cs_all[['nominal_p_value_cutoff']] <- NULL

# set the Excel output location
eqtl_output_excel_loc <- paste(eqtl_output_loc, 'eqtls_significant.xlsx', sep = '/')
caqtl_output_excel_loc <- paste(caqtl_output_loc, 'caqtls_significant.xlsx', sep = '/')
# and write the result
write_output_to_excel(eqtl_outputs_cs_all, eqtl_output_excel_loc)
write_output_to_excel(caqtl_outputs_cs_all, caqtl_output_excel_loc)
