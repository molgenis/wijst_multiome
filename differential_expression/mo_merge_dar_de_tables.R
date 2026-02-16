#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_merge_dar_de_tables.R
# Function: merge the DE results across cell types into an Excel, DARs as well
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

get_output_per_comparison <- function(output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), file_append='_condition_final.tsv.gz') {
  # store the results per cell type
  results_per_celltype <- list()
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the output location
    output_loc_full <- paste(output_loc, cell_type, file_append, sep = '')
    # read the table
    output <- read.table(output_loc_full, header = T, sep = '\t', row.names = 1)
    # add to the comparison list
    results_per_celltype[[cell_type]] <- output
  }
  return(results_per_celltype)
}

write_output_to_tsv <- function(de_output_loc, tsv_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, remove_columns=c(), file_append='_condition_final.tsv.gz', filter_list=NULL) {
  # get the output
  results_per_celltype <- get_output_per_comparison(de_output_loc, cell_types = cell_types, file_append = file_append)
  # check each of these
  for (ct in names(results_per_celltype)) {
    # extract
    results_ct <- results_per_celltype[[ct]]
    # add the cell type as an explicit column
    results_ct <- cbind(data.frame('cell_type' = rep(ct, times = nrow(results_ct))), results_ct)
    # put back into list
    results_per_celltype[[ct]] <-  results_ct
  }
  # now merge all
  results_all_celltype <- do.call('rbind', results_per_celltype)
  # set output loc as the tsv
  output_loc <- tsv_output_loc
  # gz file ends with .gz
  if (grepl('.gz$', tsv_output_loc)) {
    # gzip if ends with .gz
    output_loc <- gzfile(tsv_output_loc)
  }
  write.table(results_all_celltype, output_loc, sep = '\t', row.names = F, col.names = T, quote = F)
  # make a checksum
  mdfiver::create_sha256_for_file(tsv_output_loc)
}


write_output_to_excel <- function(de_output_loc, excel_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, remove_columns=c(), file_append='_condition_final.tsv.gz', filter_list=NULL) {
  # get the output
  results_per_celltype <- get_output_per_comparison(de_output_loc, cell_types = cell_types, file_append = file_append)
  # create a new workbook
  wb = createWorkbook()
  # check each cell type
  for (celltype in names(results_per_celltype)) {
    # create the name for the sheet
    sheet_name <- celltype
    # replace with a better label if requested
    if (use_label_dict) {
      sheet_name <- rename_labels(celltype)
    }
    # create the sheet
    sheet = createSheet(wb, sheet_name)
    # fetch the dataframe
    celltype_result <- results_per_celltype[[celltype]]
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


####################
# Main code        #
####################

# de loc
de_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/'
# excel loc
de_excel_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/mo_de_stim.xlsx'
# run the summary
write_output_to_excel(
  de_output_loc, 
  de_excel_loc, 
  cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), 
  use_label_dict=T, 
  remove_columns=c(), 
  file_append='_condition_final.tsv.gz'
)
# also do the significant only
de_excel_sig_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/mo_de_stim_significant.xlsx'
# run the summary
write_output_to_excel(
  de_output_loc, 
  de_excel_sig_loc, 
  cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), 
  use_label_dict=T, 
  remove_columns=c(), 
  file_append='_condition_final.tsv.gz', 
  filter_list = list('p.bonferroni' = 0.05)
)
# also a tsv
de_tsv_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/mo_de_stim.tsv.gz'
write_output_to_tsv(
  de_output_loc, 
  de_tsv_loc, 
  cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), 
  use_label_dict=T, 
  remove_columns=c(), 
  file_append='_condition_final.tsv.gz'
)

# DAR loc
dar_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/stimulation/pct01/'
# also do the significant only
dar_excel_sig_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/stimulation/pct01//mo_dar_stim_significant.xlsx'
# run the summary
write_output_to_excel(
  dar_output_loc, 
  dar_excel_sig_loc, 
  cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), 
  use_label_dict=T, 
  remove_columns=c('permuted', 'seed', 'perm.FDR'), 
  file_append='_condition_final.wpermfdr.tsv.gz', 
  filter_list = list('p.bonferroni' = 0.05)
)
# tsv as well
dar_tsv_sig_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/stimulation/pct01//mo_dar_stim_significant.tsv.gz'
# run the summary
write_output_to_tsv(
  dar_output_loc, 
  dar_tsv_sig_loc, 
  cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), 
  use_label_dict=T, 
  remove_columns=c('permuted', 'seed', 'perm.FDR'), 
  file_append='_condition_final.wpermfdr.tsv.gz', 
  filter_list = list('p.bonferroni' = 0.05)
)
# full sumstats as well
dar_tsv_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/stimulation/pct01//mo_dar_stim.tsv.gz'
# run the summary
write_output_to_tsv(
  dar_output_loc, 
  dar_tsv_loc, 
  cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), 
  use_label_dict=T, 
  remove_columns=c('permuted', 'seed', 'perm.FDR'), 
  file_append='_condition_final.wpermfdr.tsv.gz'
)
