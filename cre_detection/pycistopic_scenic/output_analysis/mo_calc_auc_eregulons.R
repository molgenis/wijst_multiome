#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_calc_auc_eregulons.R
# Function: recalculate the AUC excluding the eGenes
#
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Matrix)
library(data.table)
library(AUCell)
library(ggplot2)
library(stringr)


####################
# Functions        #
####################

#' get a label dict that replaces the posix safe names into printable versions
#' 
#' @returns a label dict that replaces the posix safe names into printable versions
#' label_dict_names <- get_label_dict()
get_label_dict <- function() {
  label_dict <- list()
  label_dict[['CD4_T_cells']] <- 'CD4+ T cells'
  label_dict[['CD8_T_cells']] <- 'CD8+ T cells'
  label_dict[['CD4T']] <- 'CD4+ T'
  label_dict[['CD8T']] <- 'CD8+ T'
  label_dict[['CD4_T']] <- 'CD4+ T'
  label_dict[['CD8_T']] <- 'CD8+ T'
  label_dict[['Dendritic_cells']] <- 'Dendritic cells'
  label_dict[['Endothelial_cells']] <- 'Endothelial cells'
  label_dict[['Fibroblasts']] <- 'Fibroblasts'
  label_dict[['Glia_cells']] <- 'Glia cells'
  label_dict[['Mast_cells']] <- 'MAST cells'
  label_dict[['Mature_absorptive_enterocytes']] <- 'Mature absorptive enterocytes'
  label_dict[['Mature_secretory_enterocytes']] <- 'Mature secretory enterocytes'
  label_dict[['Memory_B']] <- 'Memory B cells'
  label_dict[['Monocytes']] <- 'Monocytes'
  label_dict[['monocyte']] <- 'Monocyte'
  label_dict[['Mono']] <- 'Monocyte'
  label_dict[['Plasma_cells']] <- 'Plasma cells'
  label_dict[['Stem_cells']] <- 'Stem cells'
  label_dict[['Stromal_cells']] <- 'Stromal cells'
  label_dict[['T_others']] <- 'other T cells'
  label_dict[['Transit_amplifying_cells']] <- 'Transit amplifying cells'
  label_dict[['AI']] <- 'Actively Inflamed'
  label_dict[['NI']] <- 'Non-Inflamed'
  return(label_dict)
}


rename_labels <- function(vector_to_rename) {
  # get the labels that are present
  label_renaming <- get_label_dict()
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


remap_with_label_dict <- function(vector_of_names) {
  # get the label dict
  relabels <- get_label_dict()
  # get the labels available for renaming
  labels_available <- names(relabels)
  # get the ones we cant remap
  unmappable <- setdiff(unique(vector_of_names), labels_available)
  # report on those
  if (length(unmappable) > 0) {
    print(paste('cannot remap the following names, they will be returned unchanged:', paste(unmappable, collapse = ',')))
    # and put those in our remapping list as their originals
    relabels[unmappable] <- unmappable
  }
  # now actually do the remapping
  remapped <- as.vector(unlist(relabels[vector_of_names]))
  return(remapped)
}


get_color_coding_dict <- function() {
  # medhigh
  color_coding_dict <- list()
  color_coding_dict[["B"]] <- "#71BC4B"
  #color_coding_dict[['CD4_T_cells']] <- '#7FC97F'
  color_coding_dict[['CD4_T_cells']] <- '#153057'
  color_coding_dict[['CD4T']] <- '#153057'
  #color_coding_dict[['CD8_T_cells']] <- '#BEAED4'
  color_coding_dict[['CD8_T_cells']] <- '#009DDB'
  color_coding_dict[['CD8T']] <- '#009DDB'
  #color_coding_dict[['Dendritic_cells']] <- '#FDC086'
  color_coding_dict[['Dendritic_cells']] <- '#965EC8'
  color_coding_dict[['DC']] <- '#965EC8'
  color_coding_dict[['Endothelial_cells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia_cells']] <- '#F0027F'
  color_coding_dict[['Mast_cells']] <- '#BF5B17'
  color_coding_dict[['Mature_absorptive_enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature_secretory_enterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory_B']] <- '#D95F02'
  color_coding_dict[['Microfold_cell']] <- '#BEAED4'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Monocyte']] <- '#EDBA1B'
  color_coding_dict[['Naive_B_cells']] <- '#FDC086'
  color_coding_dict[['NK']] <- '#E64B50'
  #color_coding_dict[['Plasma_cells']] <- '#E7298A'
  color_coding_dict[['Plasma_cells']] <- '#DB8E00'
  color_coding_dict[['Stem_cells']] <- '#66A61E'
  color_coding_dict[['Stromal_cells']] <- '#8DD3C7'
  #color_coding_dict[['T_others']] <- '#A6761D'
  color_coding_dict[['T_others']] <- '#FF63B6'
  color_coding_dict[['T_other']] <- '#FF63B6'
  color_coding_dict[['Transit_amplifying_cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  #color_coding_dict[['CD4+ T cells']] <- '#7FC97F'
  color_coding_dict[['CD4+ T cells']] <- '#153057'
  color_coding_dict[['CD4+ T']] <- '#153057'
  #color_coding_dict[['CD8+ T cells']] <- '#BEAED4'
  color_coding_dict[['CD8+ T cells']] <- '#009DDB'
  color_coding_dict[['CD8+ T']] <- '#009DDB'
  #color_coding_dict[['Dendritic cells']] <- '#FDC086'
  color_coding_dict[['Dendritic cells']] <- '#965EC8'
  color_coding_dict[['Endothelial cells']] <- '#FFFFB3'
  color_coding_dict[['Endothelial\ncells']] <- '#FFFFB3'
  color_coding_dict[['Fibroblasts']] <- '#386CB0'
  color_coding_dict[['Glia cells']] <- '#F0027F'
  color_coding_dict[['MAST cells']] <- '#BF5B17'
  color_coding_dict[['Mature absorptive enterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature\nabsorptive\nenterocytes']] <- '#A6CEE3'
  color_coding_dict[['Mature secretory enterocytes']] <- '#1B9E77'
  color_coding_dict[['Mature secretory\nenterocytes']] <- '#1B9E77'
  color_coding_dict[['Memory B cells']] <- '#D95F02'
  #color_coding_dict[['Monocytes']] <- '#7570B3'
  color_coding_dict[['Microfold cells']] <- '#BEAED4'
  color_coding_dict[['Monocytes']] <- '#EDBA1B'
  color_coding_dict[['Naive B cells']] <- '#FDC086'
  #color_coding_dict[['Plasma cells']] <- '#E7298A'
  color_coding_dict[['Plasma cells']] <- '#DB8E00'
  color_coding_dict[['Plasmablast']] <- '#DB8E00'
  color_coding_dict[['plasmablast']] <- '#DB8E00'
  color_coding_dict[['Stem cells']] <- '#66A61E'
  color_coding_dict[['Stromal cells']] <- '#8DD3C7'
  #color_coding_dict[['other T cells']] <- '#A6761D'
  color_coding_dict[['other T cells']] <- '#FF63B6'
  color_coding_dict[['Transit amplifying cells']] <- '#FF7F00'
  color_coding_dict[['Transit\namplifying cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  # up and down regulation will be added to, we need a whitening percentage
  pct_whitening <- 40
  # then we will check each cell type
  for (cell_type in names(color_coding_dict)) {
    # the up color is the same as the regular one
    color_coding_dict[[paste(cell_type, 'up')]] <- color_coding_dict[[cell_type]]
    # but the down one will have a more faded colour
    color_coding_dict[[paste(cell_type, 'down')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    # we'll do something similiar when we have multiple conditions
    color_coding_dict[[paste(cell_type, 'combined')]] <- color_coding_dict[[cell_type]]
    color_coding_dict[[paste(cell_type, 'UT')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "white"))(100)[pct_whitening]
    color_coding_dict[[paste(cell_type, '24hCA')]] <- colorRampPalette(c(color_coding_dict[[cell_type]], "black"))(100)[pct_whitening]
  }
  # general
  color_coding_dict[['AI']] <- 'darkblue'
  color_coding_dict[['NI']] <- 'darkred'
  color_coding_dict[['Actively Inflamed']] <- 'darkblue'
  color_coding_dict[['Non-Inflamed']] <- 'darkred'
  return(color_coding_dict)
}




####################
# Main code        #
####################

# location of the AUC matrix
auc_mtx_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc.mtx.gz'
# location of the barcodes
auc_barcodes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_barcodes.txt.gz'
# and the eregulon names
ereg_names_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_eregnames.txt.gz'

# location of the SCENIC outout
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both_filtered.tsv.gz'
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
# read the matrix
auc_mtx <- Matrix::readMM(auc_mtx_loc)
# read the barcodes and eregnames
barcodes <- fread(auc_barcodes_loc, header = F)[[1]]
eregnames <- fread(ereg_names_loc, header = F)[[1]]
# keep only the eregulons that we kept in our output
auc_mtx <- auc_mtx[eregnames %in% scenic_output[['Gene_signature_name']], ]
eregnames <- eregnames[eregnames %in% scenic_output[['Gene_signature_name']]]
# construct matrix in Seurat format
colnames(auc_mtx) <- barcodes
rownames(auc_mtx) <- eregnames

# location of the Seurat object
seurat_object_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240223_seuratv5_normalized.rds'
# laod the data
seurat_object <- readRDS(seurat_object_loc)

# make a list of all the genes per ereg
ereg_to_genes <- list()
# for (ereg in unique(scenic_output[['Gene_signature_name']])) {
#   # get the genes
#   genes_ereg <- unique(scenic_output[
#     scenic_output[['Gene_signature_name']] == ereg, 
#   ][['Gene']])
#   # put in the list
#   ereg_to_genes[[ereg]] <- genes_ereg
# }
# location of the eregulon from scenics
scenic_eregs_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_signatures.tsv.gz'
# read that
scenic_eregs <- fread(scenic_eregs_loc, header = T, sep = '\t')
# check each of the gene based ones
for (ereg in unique(scenic_eregs[scenic_eregs[['modality']] == 'Gene_based', ][['signature_name']])) {
    # get the genes
    genes_ereg <- unique(scenic_eregs[
      scenic_eregs[['signature_name']] == ereg,
    ][['gene_or_region']])
    # put in the list
    ereg_to_genes[[ereg]] <- unique(genes_ereg)
}

# extract expression matrix
expr_mat <- seurat_object@assays$RNA@layers$counts
# set the barcodes as column names
colnames(expr_mat) <- colnames(seurat_object)
# extract genes
gene_names <- data.frame(seurat_object@assays$RNA@features)
gene_names_counts <- rownames(gene_names[gene_names[['counts']], , drop = F])
# set those as row names
rownames(expr_mat) <- gene_names_counts

# Calculate enrichment scores
cells_AUC <- AUCell_run(expr_mat, ereg_to_genes)
# set the path to save this one
cells_AUC_loc <- '~/multiome/rds/mo_auc_mtx_full.rds'
# save the file
saveRDS(cells_AUC, cells_AUC_loc)
# make a checksum
mdfiver::create_sha256_for_file(cells_AUC_loc)

# subset the expression matrix
expr_mat_included <- expr_mat[, colnames(expr_mat) %in% barcodes]
# Calculate enrichment scores
cells_AUC_included <- AUCell_run(expr_mat_included, ereg_to_genes)
# set the path to save this one
cells_AUC_included_loc <- '~/multiome/rds/mo_auc_mtx_included_fromscenic.rds'
# save the file
saveRDS(cells_AUC_included, cells_AUC_included_loc)
# make a checksum
mdfiver::create_sha256_for_file(cells_AUC_included_loc)

# reload data
cells_AUC_loc <- '~/multiome/rds/mo_auc_mtx_full.rds'
cells_AUC_included_loc <- '~/multiome/rds/mo_auc_mtx_included_fromscenic.rds'
cells_AUC <- readRDS(cells_AUC_loc)
cells_AUC_included <- readRDS(cells_AUC_included_loc)

# first do AUC we had done before vs the new full one
barcodes_before_vs_full <- intersect(barcodes, colnames(cells_AUC))
# then subset
cells_AUC_vs_before <- getAUC(cells_AUC[, barcodes_before_vs_full])
cells_before_vs_AUC <- auc_mtx[, barcodes_before_vs_full]
# then do AUC we had before vs the new limited one
barcodes_before_vs_included <- intersect(barcodes, colnames(cells_AUC_included))
# then subset
cells_AUC_included_vs_before <- getAUC(cells_AUC_included[, barcodes_before_vs_included])
cells_before_vs_AUC_included <- auc_mtx[, barcodes_before_vs_included]

# keep the eregs also in SCENIC+
ereg_to_genes <- ereg_to_genes[names(ereg_to_genes) %in% scenic_eregs_to_keep[['Gene_signature_name']]]

# create a list to store each df
auc_cors <- list()
# get correlation for each gene set
for (gene_set in names(ereg_to_genes)) {
  print(gene_set)
  # do the full matrix vs what we had before
  auc_vs_before_cor <- cor(
    x = as.vector(unlist(cells_AUC_vs_before[gene_set, ])), 
    y = as.vector(unlist(cells_before_vs_AUC[gene_set, ]))
  )
  # do the limited matrix vs what we had before
  auc_included_vs_before_cor <- cor(
    x = as.vector(unlist(cells_AUC_included_vs_before[gene_set, ])), 
    y = as.vector(unlist(cells_before_vs_AUC_included[gene_set, ]))
  )
  # and how many genes where in here
  auc_n_genes <- length(ereg_to_genes[[gene_set]])
  # make into a df
  cor_df_gs <- data.frame('regulon' = c(gene_set), 'cor_full' = c(auc_vs_before_cor), 'cor_included' = c(auc_included_vs_before_cor), 'n_genes' = c(auc_n_genes))
  # put in the list
  auc_cors[[gene_set]] <- cor_df_gs
}
# merge all
auc_cors_all <- do.call('rbind', auc_cors)
# add the gene list number
auc_cors_all[['n_gene_original']] <- apply(auc_cors_all, 1, function(x) {
  # extract value
  str_res <- str_extract(x[['regulon']], '_\\(\\d+g\\)$')[[1]]
  # remove the characters that are not numeric
  str_res <- gsub('\\(|\\)|g|_', '', str_res)
  # # then make numeric
  num_res <- as.numeric(str_res)
  return(num_res)
})
# add pct
auc_cors_all[['n_genes_pct']] <- auc_cors_all[['n_genes']] / auc_cors_all[['n_gene_original']]
# show the correlation versus the percentage
ggplot(data = auc_cors_all, mapping = aes(x = n_genes_pct, y = cor_included)) +
  geom_point() +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('% of genes of original gene set') + 
  ylab('Correlation of SCENIC+ vs AUCell values')
# show the correlation versus the percentage again for 
ggplot(data = auc_cors_all, mapping = aes(x = n_genes_pct, y = cor_full)) +
  geom_point() +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('% of genes of original gene set') + 
  ylab('Correlation of SCENIC+ vs AUCell values')
ggplot(data = auc_cors_all, mapping = aes(x = regulon, y = cor_full)) +
  geom_point() +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('eregulon name') + 
  ylab('Correlation of SCENIC+ vs AUCell values') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# read the TF-interaction-QTL output
tf_i_eqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls/merged/results_fdr.tsv.gz'
tf_i_eqtl <- fread(tf_i_eqtl_loc,  header = T, sep = '\t')
# get what is significant
sc_tf_ieqtl_sig <- tf_i_eqtl[
  tf_i_eqtl[['region:genotype_bh']] < 0.05 & 
    tf_i_eqtl[['anova_bh']] < 0.05 & 
    tf_i_eqtl[['region_bh']] < 0.05 & 
    tf_i_eqtl[['genotype_bh']] < 0.05, 
]
# get the top TF for each gene
sc_tf_ieqtl_sig_top_tf_per_gene <- sc_tf_ieqtl_sig[order(sc_tf_ieqtl_sig[['region:genotype_p']]), ]
sc_tf_ieqtl_sig_top_tf_per_gene <- sc_tf_ieqtl_sig_top_tf_per_gene[!duplicated(sc_tf_ieqtl_sig_top_tf_per_gene[['gene']]), ]
# remove the top gene from the eregulon lists
ereg_to_genes_notop <- list()
for (ereg in names(ereg_to_genes)) {
  # check if there is an interaction
  if (ereg %in% sc_tf_ieqtl_sig_top_tf_per_gene[['region']]) {
    # get the genes
    genes_ereg <- ereg_to_genes[[ereg]]
    # get the top gene
    top_gene <- unlist(as.vector(sc_tf_ieqtl_sig_top_tf_per_gene[sc_tf_ieqtl_sig_top_tf_per_gene[['region']] == ereg, ][['gene']][[1]]))
    # and remove one
    genes_ereg_notop <- setdiff(genes_ereg, top_gene)
    # put in the list
    ereg_to_genes_notop[[ereg]] <- genes_ereg_notop
  }
}
# get the rankings
expr_mat_included_ranks <- AUCell_buildRankings(expr_mat_included, plotStats=TRUE)
# Calculate enrichment scores
cells_AUC_notop <- AUCell_run(expr_mat_included, ereg_to_genes_notop)
# set the path to save this one
cells_AUC_notop_loc <- '~/multiome/rds/mo_auc_mtx_notop_fromscenic.rds'
# save the file
saveRDS(cells_AUC_notop, cells_AUC_notop_loc)
# make a checksum
mdfiver::create_sha256_for_file(cells_AUC_notop_loc)


# first do AUC we had done before vs the new full one
barcodes_full_vs_notop <- intersect(colnames(cells_AUC), colnames(cells_AUC_notop))
# then subset
cells_AUC_notop_vs_full <- getAUC(cells_AUC_notop[, barcodes_full_vs_notop])
cells_AUC_full_vs_notop <- getAUC(cells_AUC[, barcodes_full_vs_notop])

# first do AUC we had done before vs the new full one
barcodes_included_vs_notop <- intersect(colnames(cells_AUC_included), colnames(cells_AUC_notop))
# then subset
cells_AUC_notop_vs_included <- getAUC(cells_AUC_notop[, barcodes_included_vs_notop])
cells_AUC_included_vs_notop <- getAUC(cells_AUC_included[, barcodes_included_vs_notop])
# get correlation for each gene set
auc_cors_notop <- list()
for (gene_set in names(ereg_to_genes_notop)) {
  print(gene_set)
  # do the full matrix vs what we had before
  auc_notop_vs_full_cor <- cor(
    x = as.vector(unlist(cells_AUC_notop_vs_full[gene_set, ])), 
    y = as.vector(unlist(cells_AUC_full_vs_notop[gene_set, ]))
  )
  auc_notop_vs_included_cor <- cor(
    x = as.vector(unlist(cells_AUC_notop_vs_included[gene_set, ])), 
    y = as.vector(unlist(cells_AUC_included_vs_notop[gene_set, ]))
  )
  # and how many genes where in here
  auc_n_genes <- length(ereg_to_genes[[gene_set]])
  # make into a df
  cor_df_gs <- data.frame('regulon' = c(gene_set), 'cor_full_notop' = c(auc_notop_vs_before_cor), 'cor_included_notop' = c(auc_notop_vs_included_cor), 'n_genes' = c(auc_n_genes))
  # put in the list
  auc_cors_notop[[gene_set]] <- cor_df_gs
}
# merge all
auc_cors_notop_all <- do.call('rbind', auc_cors_notop)
ggplot(data = auc_cors_notop_all, mapping = aes(x = regulon, y = cor_included_notop )) +
  geom_point() +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('eregulon name') + 
  ylab('Correlation of AUCell values before and after top TF-i-eQTL gene') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# we'll do some things related to metadata
metadata_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/mo_celllevel_metadata_nounannotated.tsv.gz'
metadata <- fread(metadata_loc, header = T, sep = '\t')
# add nicer cell type name
metadata[['cell_type_nice']] <- rename_labels(as.character(metadata[['celltype_imputed_lowerres']]))
# add a new column of ct to stim
metadata[['cell_type_stim']] <- paste(as.character(metadata[['cell_type_nice']]), as.character(metadata[['condition_final']]))
# subset to some columns, so we can use those to create an ordering
ct_to_stim <- unique(metadata[, c('condition_final', 'cell_type_nice', 'cell_type_stim')])
# order by condition in the way that we like
ct_to_stim[['condition_final']] <- factor(ct_to_stim[['condition_final']], levels = c('UT', '24hCA'))
# now order the cell types and conditions, having made sure that UT is before 24hCA
ct_to_stim <- ct_to_stim[order(ct_to_stim[['cell_type_nice']], ct_to_stim[['condition_final']]), ]
# and make the cell type stim order like this in the full metadata
metadata[['cell_type_stim']] <- factor(metadata[['cell_type_stim']], levels = ct_to_stim[['cell_type_stim']])
# subset the metadata to the AUC data we have
metadata_auc_barcodes <- intersect(barcodes, metadata[['barcode_lane']])
metadata_to_auc <- metadata[match(metadata_auc_barcodes, metadata[['barcode_lane']]), ]
# and the AUC data as well
auc_to_metadata <- auc_mtx[, metadata_auc_barcodes]
# get a list of eregulons to plot
eregs_to_plot <- c('EBF1_direct_+/+_(142g)', 'PAX5_direct_+/+_(115g)')
# store the plots
ereg_plots <- list()
# do each ereg
for (ereg in eregs_to_plot) {
  # add this ereg the metadata
  metadata_to_auc[['tf_activity']] <- as.vector(unlist(auc_to_metadata[ereg, ]))
  # and plot
  p_ereg <- ggplot(data = metadata_to_auc, mapping = aes(x = cell_type_stim, y = tf_activity, fill = cell_type_stim)) + 
    geom_boxplot(outlier.shape = NA) + 
    # and add jitter
    geom_jitter(size = 0.5, alpha = 0.1) +
    # and labels
    xlab('Cell type and condition') + 
    ylab(paste(ereg, 'activity')) + 
    theme(legend.position = 'none') +
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual(values = get_color_coding_dict())
  # put in list
  ereg_plots[[ereg]] <- p_ereg
}

# remove the top gene from the eregulon lists
ereg_to_genes_iegenes <- list()
for (ereg in names(ereg_to_genes)) {
  # check if there is an interaction
  if (ereg %in% sc_tf_ieqtl_sig[['region']]) {
    # get the genes
    genes_ereg <- ereg_to_genes[[ereg]]
    # get the genes
    tf_i_egenes <- unlist(as.vector(sc_tf_ieqtl_sig_top_tf_per_gene[sc_tf_ieqtl_sig_top_tf_per_gene[['region']] == ereg, ][['gene']]))
    # check each gene
    for (tf_i_egene in tf_i_egenes) {
      # and thate one
      genes_ereg_no_tf_i_egene <- setdiff(genes_ereg, tf_i_egene)
      # put in the list
      ereg_to_genes_iegenes[[paste(ereg, tf_i_egene, sep = '_')]] <- genes_ereg_no_tf_i_egene
    }
  }
}
# recalculate enrichment scores
cells_AUC_notop <- AUCell_run(expr_mat_included, ereg_to_genes_iegenes)
# extract the auc matrix
cells_AUC_notop_mtx <- getAUC(cells_AUC_notop)
# as matrix
cells_AUC_notop_mtx <- as.matrix(cells_AUC_notop_mtx)
# extract dimension names
cells_AUC_notop_eregs <- rownames(cells_AUC_notop_mtx)
cells_AUC_notop_barcodes <- colnames(cells_AUC_notop_mtx)
# make into data.table
cells_AUC_notop_dt <- data.table(cells_AUC_notop_mtx)
# but add the eregulon as the first column
cells_AUC_notop_dt <- cbind(data.table('eregulon' = cells_AUC_notop_eregs), cells_AUC_notop_dt)
# set export location
cells_AUC_notop_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/mo_tf_interaction_eregulon_gene_removed_auc_all_nonsparse_transposed.tsv.gz'
# write the output
write.table(cells_AUC_notop_dt, gzfile(cells_AUC_notop_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# with a checksum
mdfiver::create_sha256_for_file(cells_AUC_notop_loc)

# make a new confinement file
sc_tf_ieqtl_sig_confinement <- sc_tf_ieqtl_sig[, c('variant', 'region', 'gene')]
# and replace the region with the eregulon name and the gene, as we generated for the new AUC matrix
sc_tf_ieqtl_sig_confinement[['region']] <- paste(sc_tf_ieqtl_sig[['region']], sc_tf_ieqtl_sig[['gene']], sep = '_')
# update the column names to actually reflect what we do
colnames(sc_tf_ieqtl_sig_confinement) <- c('variant', 'eregulon', 'gene')
# set output loc
sc_tf_ieqtl_sig_confinement_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/mo_var_tf_gene_confinement_inclcaqtls_significant.tsv.gz'
# write output
write.table(sc_tf_ieqtl_sig_confinement, gzfile(sc_tf_ieqtl_sig_confinement_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# with a checksum
mdfiver::create_sha256_for_file(sc_tf_ieqtl_sig_confinement_loc)
