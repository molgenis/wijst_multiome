#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_differential_tf_activity.R
# Function: get the differential inferred TF activity between lineages and cell types, and plot it in a waterfall plot
#
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(cowplot)

####################
# Functions        #
####################


postprocess_tf_markers <- function(tf_markers_table) {
  # add the eregulon as the first column
  tf_markers_table <- cbind('eregulon' = rownames(tf_markers_table), tf_markers_table)
  # but fix the names
  tf_markers_table[['eregulon']] <- gsub('--', ' -', tf_markers_table[['eregulon']])
  tf_markers_table[['eregulon']] <- gsub('-\\+', ' +', tf_markers_table[['eregulon']])
  # set the p value to the smallest value if it is zero
  if (sum(tf_markers_table[['p_val']] == 0) > 0) {
    tf_markers_table[tf_markers_table[['p_val']] == 0, ][['p_val']] <- .Machine$double.xmin
  }
  # redo MTC
  tf_markers_table[['p_val_adj']] <- p.adjust(tf_markers_table[['p_val']], method = 'BH')
  # take -log10 of that p value
  tf_markers_table[['p_val_adj_neg_log10']] <- -log10(tf_markers_table[['p_val_adj']])
  # add a directed p value
  tf_markers_table[['p_val_adj_neg_log10_directed']] <- sign(tf_markers_table[['avg_log2FC']]) * tf_markers_table[['p_val_adj_neg_log10']]
  # order by alphabet, then lfc, then p value
  tf_markers_table <- tf_markers_table[order(tf_markers_table[['eregulon']]), ]
  tf_markers_table <- tf_markers_table[order(tf_markers_table[['p_val']], decreasing = F), ]
  tf_markers_table <- tf_markers_table[order(tf_markers_table[['avg_log2FC']], decreasing = T), ]
  # now set the eregulon as a factor, so when we plot it, we have that order
  tf_markers_table[['eregulon']] <- factor(tf_markers_table[['eregulon']], levels = tf_markers_table[['eregulon']])
  return(tf_markers_table)
}

plot_waterfall <- function(postprocessed_tf_markers_table, ereg_column='eregulon', logfc_column='avg_log2FC', pval_column='p_val_adj_neg_log10_directed', low_color='darkblue', mid_color='white', high_color='darkred', legendless=F, paper_style=T, angle_x_labels=T) {
  # set some column names so we can use them consistently
  postprocessed_tf_markers_table[['eregulon']] <- postprocessed_tf_markers_table[[ereg_column]]
  postprocessed_tf_markers_table[['log2FC']] <- postprocessed_tf_markers_table[[logfc_column]]
  postprocessed_tf_markers_table[['significance']] <- postprocessed_tf_markers_table[[pval_column]]
  # create a color scale
  color_scale <- scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, midpoint = 0)
  p <- ggplot(data = postprocessed_tf_markers_table, mapping = aes(x = eregulon, y = log2FC, fill = significance)) + geom_bar(stat = 'identity') + 
    color_scale + 
    labs(fill = pval_column)
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  if (paper_style) {
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  if (angle_x_labels) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  return(p)
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
# and keep only what is non-extended if it was both extended and non-extended
scenic_output <- scenic_output[!duplicated(paste(scenic_output[['TF']], scenic_output[['Gene_signature_direction']])), ]


# the location of the metadata
metadata_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/mo_celllevel_metadata_nounannotated.tsv.gz'
# read the metadata
metadata <- fread(metadata_loc, header = T, sep = '\t')

# read the matrix
auc_mtx <- Matrix::readMM(auc_mtx_loc)
# read the barcodes and eregnames
barcodes <- fread(auc_barcodes_loc, header = F)[[1]]
eregnames <- fread(ereg_names_loc, header = F)[[1]]

# keep only the eregulons that we kept in our output
auc_mtx <- auc_mtx[eregnames %in% scenic_output[['Gene_signature_name']], ]
eregnames <- eregnames[eregnames %in% scenic_output[['Gene_signature_name']]]
# rename the eregs to remove some info that is not relevant
eregnames <- gsub('_direct|_extended', '', eregnames)
eregnames <- gsub('_\\(\\d+g)$', '', eregnames)

# check which barcodes are in the metadata and this object
barcodes_both <- intersect(barcodes, metadata[['barcode_lane']])
# keep the barcodes that are in the 'both'
auc_mtx <- auc_mtx[, barcodes %in% barcodes_both]
barcodes <- barcodes[barcodes %in% barcodes_both]
# and make sure the metadata is in the same order now
metadata <- metadata[match(barcodes, metadata[['barcode_lane']]), ]

# construct matrix in Seurat format
colnames(auc_mtx) <- barcodes
rownames(auc_mtx) <- eregnames
auc_assay <- CreateAssay5Object(data = auc_mtx)
# then a Seurat object
mo_tf <- CreateSeuratObject(counts = auc_assay, assay = 'TF', meta.data = metadata, project = 'wijst_multiome')
# add to the metadata myeloid and lymphoid info
mo_tf@meta.data[['lineage']] <- NA
# get classification of lymphoid or myeloid
myeloid_cts <- c('DC', 'monocyte')
lymphoid_cts <- c('B', 'CD4T', 'CD8T', 'NK', 'plasmablast', 'T_other')
# then do actual annotation
mo_tf@meta.data[!is.na(mo_tf@meta.data[['celltype_imputed_lowerres']]) & mo_tf@meta.data[['celltype_imputed_lowerres']] %in% lymphoid_cts, ][['lineage']] <- 'lymphoid'
mo_tf@meta.data[!is.na(mo_tf@meta.data[['celltype_imputed_lowerres']]) & mo_tf@meta.data[['celltype_imputed_lowerres']] %in% myeloid_cts, ][['lineage']] <- 'myeloid'

# do DE-like analysis
tf_lineage_markers <- FindMarkers(object = mo_tf, ident.1 = 'lymphoid', ident.2 = 'myeloid', group.by = 'lineage', assay = 'TF', test.use = 'wilcox', logfc.threshold = 0, min.pct = 0)
# do postprocessing before we plot
tf_lineage_markers <- postprocess_tf_markers(tf_lineage_markers)
# show this for the markers
plot_grid(
  plot_waterfall(tf_lineage_markers, legendless = T) + 
    xlab('eRegulon') + 
    ylab('log2FC'),
  get_legend(plot_waterfall(tf_lineage_markers, pval_column = 'p_val_adj_neg_log10', high_color = 'black') + labs(fill = '-log10(p)')), 
  nrow = 1, 
  ncol = 2, 
  rel_widths = c(7,1)
)
# store results per lineage
ct_vs_nonlineage <- list()
# do this repeatedly for each of the cell typed
for (cell_type in unique(mo_tf@meta.data[['celltype_imputed_lowerres']])) {
  # get group it is in
  lineage_ct <- NA
  if (cell_type %in% lymphoid_cts) {
    lineage_ct <- 'lymphoid'
  }
  else if (cell_type %in% myeloid_cts) {
    lineage_ct <- 'myeloid'
  }
  # subset the data to that cell type and the opposite lineage
  mo_tf_ct_vs_lineage <- mo_tf[, 
      !is.na(mo_tf@meta.data[['celltype_imputed_lowerres']]) &
      !is.na(mo_tf@meta.data[['lineage']]) & 
      (mo_tf@meta.data[['celltype_imputed_lowerres']] == cell_type | mo_tf@meta.data[['lineage']] != lineage_ct)
  ]
  # do the comparison
  tf_ct_vs_opposite_lineage_markers <- FindMarkers(object = mo_tf_ct_vs_lineage, ident.1 = cell_type, ident.2 = NULL, group.by = 'celltype_imputed_lowerres', assay = 'TF', test.use = 'wilcox', logfc.threshold = 0, min.pct = 0)
  # postprocess
  tf_ct_vs_opposite_lineage_markers <- postprocess_tf_markers(tf_ct_vs_opposite_lineage_markers)
  # set the ident names
  tf_ct_vs_opposite_lineage_markers <- cbind(data.frame('ident.1' = rep(cell_type, times = nrow(tf_ct_vs_opposite_lineage_markers)), 'ident.2' = rep(paste('not', lineage_ct), times = nrow(tf_ct_vs_opposite_lineage_markers))), tf_ct_vs_opposite_lineage_markers)
  # put in the list
  ct_vs_nonlineage[[paste(cell_type, 'vs_non', lineage_ct, sep = '_')]] <- tf_ct_vs_opposite_lineage_markers
}
# show this for the markers
plot_grid(
  plot_waterfall(ct_vs_nonlineage[['B_vs_non_lymphoid']], legendless = T) + 
    xlab('eRegulon') + 
    ylab('log2FC'),
  get_legend(plot_waterfall(ct_vs_nonlineage[['B_vs_non_lymphoid']], pval_column = 'p_val_adj_neg_log10', high_color = 'black') + labs(fill = '-log10(p)')), 
  nrow = 1, 
  ncol = 2, 
  rel_widths = c(7,1)
)

# merge all markers
ct_vs_nonlineage_all <- do.call('rbind', ct_vs_nonlineage)
# keep only most significant ones
ct_vs_nonlineage_all_lfccut <- ct_vs_nonlineage_all[ct_vs_nonlineage_all[['avg_log2FC']] > 2.5, ]
# save these two tables
ct_vs_nonlineage_all_loc <- '~/multiome/tables/mo_differential_ereg_activity.tsv.gz'
ct_vs_nonlineage_all_lfccut_loc <- '~/multiome/tables/mo_differential_ereg_activity_lfc25.tsv.gz'
write.table(ct_vs_nonlineage_all, gzfile(ct_vs_nonlineage_all_loc), row.names = F, col.names = T, sep = '\t')
write.table(ct_vs_nonlineage_all_lfccut, gzfile(ct_vs_nonlineage_all_lfccut_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_sha256_for_file(ct_vs_nonlineage_all_loc)
mdfiver::create_sha256_for_file(ct_vs_nonlineage_all_lfccut_loc)
