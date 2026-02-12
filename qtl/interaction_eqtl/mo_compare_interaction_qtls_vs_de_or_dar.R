#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_compare_interaction_qtls_vs_de_or_dar.R
# Function: plot the number of DARs vs i-caQTLs and DE vs i-eQTLs
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(ggplot2)


####################
# Functions        #
####################

#' get the number eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @param gene_column which column to use as the gene identifier
#' @returns a dataframe with the cell type in the 'cell_type' column, and the number of eGenes in the 'nr' column
#' 
get_egene_numbers_per_celltype_limix <- function(qtl_output_loc, output_file='qtl_results_all.txt.gz', gene_column='feature_id', significance_column='feature_q_value', significance_cutoff=0.05, verbose=T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(qtl_output_loc, full.names = F, recursive = F)
  # we will store the results in a list for now
  numbers_per_celltype <- list()
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
      # check each column
      for (significance_column_i in 1 : length(significance_column)) {
        cell_type_output <- cell_type_output[
          !is.na(cell_type_output[[significance_column[[significance_column_i]]]]) &
            cell_type_output[[significance_column[[significance_column_i]]]] < significance_cutoff[[significance_column_i]], 
        ]
      }
      if (verbose) {
        print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
      }
    }
    # get the unique genes in this file
    unique_genes <- unique(cell_type_output[[gene_column]])
    # get how many these are
    nr_genes <- length(unique_genes)
    # add to the list
    numbers_per_celltype[[cell_type]] <- nr_genes
  }
  # turn into a dataframe
  numbers_per_celltype_df <- data.frame(cell_type = names(numbers_per_celltype), nr = as.vector(unlist(numbers_per_celltype)))
  return(numbers_per_celltype_df)
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




de_genes_number_to_table <- function(limma_output_mo_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column = 'p.bonferroni', pval_cutoff=0.05, lfc_column = 'logFC', lfc_cutoff = 0.1, only_positive = F, only_negative = F, file_append='_condition_final.tsv.gz'){
  # get DE outputs
  results_per_celltype <- get_output_per_comparison(limma_output_mo_loc, cell_types=cell_types, file_append=file_append)
  # get numbers per cell type
  de_nr_per_ct_l <- list()
  # check each of these outputs
  for (ct in names(results_per_celltype)) {
    # grab for cell type
    results_celltype <- results_per_celltype[[ct]]
    # filter by significance
    if (!is.null(pval_column)) {
      results_celltype <- results_celltype[
        results_celltype[[pval_column]] < pval_cutoff, 
      ]
    }
    # filter by lfc
    if (!is.null(lfc_column) & !is.null(lfc_cutoff)) {
      results_celltype <- results_celltype[
        abs(results_celltype[[lfc_column]]) >= lfc_cutoff, 
      ]
    }
    # and direction
    if (only_positive) {
      results_celltype <- results_celltype[
        results_celltype[[lfc_column]] > 0, 
      ]
    }
    if (only_negative) {
      results_celltype <- results_celltype[
        results_celltype[[lfc_column]] < 0, 
      ]
    }
    # get the number of genes
    de_nr_per_ct_l[[ct]] <- data.frame('cell_type' = c(ct), 'nr' = c(nrow(results_celltype)))
  }
  # merge all
  de_nr_per_ct <- do.call('rbind', de_nr_per_ct_l)
  return(de_nr_per_ct)
}


wide_to_long_table <- function(wide_table) {
  # turn into plottable dataframe
  category_combinations_n <- length(unique(colnames(wide_table))) * length(unique(rownames(wide_table)))
  # make new table
  long_table <- data.frame(x = rep(NA, times = category_combinations_n), y = rep(NA, times = category_combinations_n), number = rep(NA, times = category_combinations_n))
  # keep an index
  i_new <- 1
  # check each row
  for (i_row in 1:nrow(wide_table)) {
    # check each column
    for (i_column in 1:ncol(wide_table)) {
      # get the x, the row
      x <- rownames(wide_table)[i_row]
      # get the y, the column
      y <- colnames(wide_table)[i_column]
      # get the value
      number <- wide_table[i_row, i_column]
      # add in the new table
      long_table[i_new, 'x'] <- x
      long_table[i_new, 'y'] <- y
      long_table[i_new, 'number'] <- number
      # update new index
      i_new <- i_new + 1
    }
  }
  return(long_table)
}


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


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UTBaseline"]] <- "khaki2"
  color_coding[["UTt24h"]] <- "khaki4"
  color_coding[["UTt8w"]] <- "paleturquoise1"
  color_coding[["Baselinet24h"]] <- "paleturquoise3"
  color_coding[["Baselinet8w"]] <- "rosybrown1"
  color_coding[["t24ht8w"]] <- "rosybrown3"
  color_coding[["UT\nBaseline"]] <- "khaki2"
  color_coding[["UT\nt24h"]] <- "khaki4"
  color_coding[["UT\nt8w"]] <- "paleturquoise1"
  color_coding[["Baseline\nt24h"]] <- "paleturquoise3"
  color_coding[["Baseline\nt8w"]] <- "rosybrown1"
  color_coding[["t24h\nt8w"]] <- "rosybrown3"
  color_coding[["UT-Baseline"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["Baseline-t24h"]] <- "paleturquoise3"
  color_coding[["Baseline-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  color_coding[["UT-t0"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["HC-t0"]] <- "khaki2"
  color_coding[["t0-HC"]] <- "khaki2"
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["t24h-HC"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t8w-HC"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t24h-t0"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t0-t8w"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t8w-t0"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t24h-t8w"]] <- "#C00040" #"rosybrown3"
  color_coding[["t8w-t24h"]] <- "#C00040" #"rosybrown3"
  # set condition colors
  color_coding[["HC"]] <- "grey"
  color_coding[["C"]] <- "grey"
  color_coding[["Controls"]] <- "grey"
  color_coding[["t0"]] <- "pink"
  color_coding[["t24h"]] <- "red"
  color_coding[["t8w"]] <- "purple"
  color_coding[["t6-8w"]] <- "purple"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # up and down
  color_coding[["Bulk up"]] <- "black"
  color_coding[["CD4T up"]] <- "#153057"
  color_coding[["CD8T up"]] <- "#009DDB"
  color_coding[["monocyte up"]] <- "#EDBA1B"
  color_coding[["NK up"]] <- "#E64B50"
  color_coding[["B up"]] <- "#71BC4B"
  color_coding[["DC up"]] <- "#965EC8"
  color_coding[["CD4+ T up"]] <- "#153057"
  color_coding[["CD8+ T up"]] <- "#009DDB"
  # percentage of whitening 
  pct_whitening=40
  color_coding[["Bulk down"]] <- colorRampPalette(c(color_coding[["Bulk up"]], "white"))(100)[pct_whitening]
  color_coding[["CD4T down"]] <- colorRampPalette(c(color_coding[["CD4T up"]], "white"))(100)[pct_whitening]
  color_coding[["CD8T down"]] <- colorRampPalette(c(color_coding[["CD8T up"]], "white"))(100)[pct_whitening]
  color_coding[["monocyte down"]] <- colorRampPalette(c(color_coding[["monocyte up"]], "white"))(100)[pct_whitening]
  color_coding[["NK down"]] <- colorRampPalette(c(color_coding[["NK up"]], "white"))(100)[pct_whitening]
  color_coding[["B down"]] <- colorRampPalette(c(color_coding[["B up"]], "white"))(100)[pct_whitening]
  color_coding[["DC down"]] <- colorRampPalette(c(color_coding[["DC up"]], "white"))(100)[pct_whitening]
  color_coding[["CD4+ T down"]] <- colorRampPalette(c(color_coding[["CD4+ T up"]], "white"))(100)[pct_whitening]
  color_coding[["CD8+ T down"]] <- colorRampPalette(c(color_coding[["CD8+ T up"]], "white"))(100)[pct_whitening]
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}


####################
# Main code        #
####################

# location of interaction-eQTLs
ieqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/'
# location of interaction-caQTLs
icaqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_caqtl/sc-eqtlgen/output/combined_significant/L1/'
# get the number of interactions with caQTLs
icaqtl_peak_numbers <- get_egene_numbers_per_celltype_limix(icaqtl_output_loc, output_file = 'inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', significance_column = c('feature_q_value', 'feature_bf_eigen'), gene_column = 'feature', significance_cutoff = c(0.05, 0.05))
# get the number of interactions with eQTLs
ieqtl_gene_numbers <- get_egene_numbers_per_celltype_limix(ieqtl_output_loc, output_file = 'inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', significance_column = c('feature_q_value', 'feature_bf_eigen'), gene_column = 'feature', significance_cutoff = c(0.05, 0.05))

# get the locations of the DE output
limma_output_mo_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_expression/limma_dream/output/stimulation/'
# get a table of the number of DE genes
de_numbers_table_mo_limma <- de_genes_number_to_table(limma_output_mo_loc, pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_positive = F)
# get the locations of the DE output
limma_dar_output_mo_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/stimulation/pct01/'
# get a table of the number of DE genes
dar_numbers_table_mo_limma <- de_genes_number_to_table(limma_dar_output_mo_loc, pval_column = 'p.bonferroni', lfc_column = 'logFC', lfc_cutoff = 0.1, only_positive = F, file_append = '_condition_final.wpermfdr.wgene.tsv.gz')

# merge the tables
de_vs_ieqtl_numbers <- merge(de_numbers_table_mo_limma, ieqtl_gene_numbers, by = 'cell_type')
# set column names
colnames(de_vs_ieqtl_numbers) <- c('cell_type', 'n_de', 'n_ieqtls')
# merge the tables
dar_vs_icaqtl_numbers <- merge(dar_numbers_table_mo_limma, icaqtl_peak_numbers, by = 'cell_type')
# set column names
colnames(dar_vs_icaqtl_numbers) <- c('cell_type', 'n_dar', 'n_icaqtls')

# rename the cell types
de_vs_ieqtl_numbers[['cell_type']] <- rename_labels(de_vs_ieqtl_numbers[['cell_type']])
dar_vs_icaqtl_numbers[['cell_type']] <- rename_labels(dar_vs_icaqtl_numbers[['cell_type']])

# plot these
p_de_vs_ieqtl <- ggplot(
  data = de_vs_ieqtl_numbers, 
  mapping = aes(
    x = n_de, 
    y = n_ieqtls, 
    colour = cell_type
  )
) +
  geom_point(size = 5) +
  scale_colour_manual(values = get_color_coding_dict()) + 
  xlab('Number of DE genes') + 
  ylab('Number of interaction-eGenes') + 
  labs(colour = 'Cell type') + 
  xlim(c(0, 1.1 * max(de_vs_ieqtl_numbers[['n_de']]))) +
  ylim(c(0, 1.1 * max(de_vs_ieqtl_numbers[['n_ieqtls']]))) +
  theme(legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_de_vs_ieqtl
ggsave('~/plots/mo_multiome_de_vs_iegene_numbers.pdf', width = 8, height = 6, plot = p_de_vs_ieqtl)

# nex the DARs
p_dar_vs_icaqtl <- ggplot(
  data = dar_vs_icaqtl_numbers, 
  mapping = aes(
    x = n_dar, 
    y = n_icaqtls, 
    colour = cell_type
  )
) +
  geom_point(size = 5) +
  scale_colour_manual(values = get_color_coding_dict()) + 
  xlab('Number of DARs') + 
  ylab('Number of interaction-caPeaks') + 
  labs(colour = 'Cell type') + 
  xlim(c(0, 1.1 * max(dar_vs_icaqtl_numbers[['n_dar']]))) +
  ylim(c(0, 1.1 * max(dar_vs_icaqtl_numbers[['n_icaqtls']]))) +
  theme(legend.title = element_text(size=14), 
        legend.text = element_text(size=12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12),
        strip.text.x = element_text(size=12)) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
p_dar_vs_icaqtl
ggsave('~/plots/mo_multiome_dar_vs_icapeak_numbers.pdf', width = 8, height = 6, plot = p_de_vs_ieqtl)
