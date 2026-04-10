#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_correlate_expression_to_activity.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(Seurat)
library(progress)
library(ggplot2)
library(cowplot)
library(stringr)
library(ggpubr)


####################
# Functions        #
####################

normalize_mj <- function(seurat_object) {
  # get the count matrix where we have the correct cell type
  count_matrix <- GetAssayData(seurat_object, slot = "counts")
  # ignore genes that are never expressed
  count_matrix <-  count_matrix[which(rowSums(count_matrix) != 0), ]
  # create new object to store the counts in
  norm_count_matrix <- count_matrix
  # do mean sample-sum normalization
  sample_sum_info = colSums(norm_count_matrix)
  mean_sample_sum = mean(sample_sum_info)
  sample_scale = sample_sum_info / mean_sample_sum
  # divide each column by sample_scale
  norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
  if ('layers' %in% slotNames(seurat_object[['RNA']])) {
    print('using Seurat v5 style \'layer\'')
    seurat_object[['MJ']] <- CreateAssay5Object(data = norm_count_matrix)
    
  } else {
    print('using Seurat v3/4 style \'slot\'')
    seurat_object[['MJ']] <- CreateAssayObject(data = norm_count_matrix)
  }
  return(seurat_object)
}


correlate_ereg_to_tf <- function(seurat_object, mapping, expression_assay='MJ', activity_assay='tf_activity', expression_slot='data', activity_slot='data') {
  # get the number of eregulons
  n_tf <- nrow(mapping)
  # store result
  cor_tf_ereg <- data.frame(
    matrix(, nrow = n_tf, ncol = 3, dimnames = list(mapping$eRegulon, c('tf', 'eregulon', 'correlation')))
  )
  # make a progress bar
  pb <- progress_bar$new(total = n_tf)
  # initialize the progress bar
  pb$tick(0)
  # check each eregulon
  for (i in 1:n_tf) {
    # update the progress bar
    pb$tick()
    # extract the TF name and the eRegulon name
    tf <- mapping$TF[i]
    ereg <- mapping$eRegulon[i]
    # set these in the matrix
    cor_tf_ereg[['tf']][i] <- tf
    cor_tf_ereg[['eregulon']][i] <- ereg
    # check if there is a TF in the MJ assay that matches the TF in the AUC assay
    if (tf %in% rownames(seurat_object@assays[[expression_assay]]@features)) {
      # get the index of the TF
      tf_index <- which(rownames(seurat_object@assays[[expression_assay]]@features) == tf)
      # get the index of the eregulon
      ereg_index <- which(rownames(seurat_object@assays[[activity_assay]]@features) == ereg)
      # get the expression and activity values
      expr_values <- seurat_object@assays$MJ@layers[[expression_slot]][tf_index, ]
      activity_values <- seurat_object@assays[[activity_assay]]@layers[[activity_slot]][ereg_index, ]
      # calculate the correlation
      cor_value <- cor(expr_values, activity_values, method = 'spearman')
      # put in the matrix
      cor_tf_ereg[['correlation']][i] <- cor_value
    } else {
      cor_tf_ereg[['correlation']][i] <- NA
    }
  }
  return(cor_tf_ereg)
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
  color_coding_dict[['Stem cells']] <- '#66A61E'
  color_coding_dict[['Stromal cells']] <- '#8DD3C7'
  #color_coding_dict[['other T cells']] <- '#A6761D'
  color_coding_dict[['other T cells']] <- '#FF63B6'
  color_coding_dict[['Transit amplifying cells']] <- '#FF7F00'
  color_coding_dict[['Transit\namplifying cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  color_coding_dict[['UT']] <- 'gray'
  color_coding_dict[['24hCA']] <- 'darkgreen'
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


################################
# Main code #
################################

# location of the Seurat object
seurat_object_loc <- paste0('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated_agesexcovid.rds')
# read the object
seurat_object <- readRDS(seurat_object_loc)
# add MJ normalization
seurat_object <- normalize_mj(seurat_object)

# location of the AUC matrix
auc_mtx_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc.mtx.gz'
# location of the barcodes
auc_barcodes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_barcodes.txt.gz'
# and the eregulon names
ereg_names_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eregulon_gene_auc_eregnames.txt.gz'
# read the matrix
auc_mtx <- Matrix::readMM(auc_mtx_loc)
# read the barcodes and eregnames
barcodes <- fread(auc_barcodes_loc, header = F)[[1]]
eregnames <- fread(ereg_names_loc, header = F)[[1]]
# construct matrix in Seurat format
colnames(auc_mtx) <- barcodes
rownames(auc_mtx) <- eregnames
auc_assay <- CreateAssay5Object(data = auc_mtx)

# keep in the seurat object only what we have in the AUC matrix
seurat_object <- seurat_object[, barcodes]
# add the AUC assay to the seurat object
seurat_object[['tf_activity']] <- auc_assay
# create a mapping
mapping <- data.frame(
  'TF' = gsub("-(extended|direct).*", "", rownames(data.frame(data.frame(seurat_object@assays$tf_activity@features)))), 
  'eRegulon' = rownames(data.frame(data.frame(seurat_object@assays$tf_activity@features)))
)
# make a dataframe to store results
cor_tf_ereg <- correlate_ereg_to_tf(seurat_object, mapping)

# subset the object to stim
seurat_object_24hca <- seurat_object[, seurat_object$condition_final == '24hCA']
# add stim specific as well
cor_tf_ereg_stim <- correlate_ereg_to_tf(seurat_object_24hca, mapping)

# subset the object to stim
seurat_object_ut <- seurat_object[, seurat_object$condition_final == 'UT']
# add stim specific as well
cor_tf_ereg_ut <- correlate_ereg_to_tf(seurat_object_ut, mapping)

# add condition
cor_tf_ereg[['condition']] <- 'both'
cor_tf_ereg_stim[['condition']] <- '24hCA'
cor_tf_ereg_ut[['condition']] <- 'UT'
# and merge all
cor_tf_ereg_all <- do.call('rbind', list(cor_tf_ereg, cor_tf_ereg_stim, cor_tf_ereg_ut))

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
# rename the regions
scenic_output[['region_cpeaks']] <- gsub(':', '-', scenic_output[['Region']])

# get TF-i-eQTLs
tf_ieqtl_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cre_detection/limix_sc/output/tf_interaction/sc/all/lane_donor_countrna_inclcaqtls_varinregion/merged/results_fdr_with_replication_significant.tsv.gz'
tf_ieqtls <- fread(tf_ieqtl_loc, header = T, sep = '\t')
tf_ieqtls_sig <- tf_ieqtls[tf_ieqtls[['significant']], ]

# add that info to the correlations
cor_tf_ereg_all[['is_tf_ieqtl']] <- cor_tf_ereg_all[['eregulon']] %in% gsub('_', '-', tf_ieqtls_sig[['tf']])
# add if it is in SCENIC+
cor_tf_ereg_all[['is_scenic_plus']] <- cor_tf_ereg_all[['eregulon']] %in% gsub('_', '-', scenic_output[['Gene_signature_name']])
# add gene signature direction
cor_tf_ereg_all[['gene_signature_direction']] <- str_extract(cor_tf_ereg_all[['eregulon']], '\\+\\/\\+|\\-\\/\\-|\\+\\/\\-|\\-\\/\\+')
# add ordered conditoin
cor_tf_ereg_all[['condition_ordered']] <- factor(cor_tf_ereg_all[['condition']], levels = c('UT', '24hCA', 'both'))
# add corrected correlation
cor_tf_ereg_all[['correlation_corrected']] <- cor_tf_ereg_all[['correlation']]
cor_tf_ereg_all[cor_tf_ereg_all[['gene_signature_direction']] %in% c('-/-', '-/+'), ][['correlation_corrected']] <- -1 * cor_tf_ereg_all[cor_tf_ereg_all[['gene_signature_direction']] %in% c('-/-', '-/+'), ][['correlation']]
# add a nice name by going back to the original name and replacing the underscore
cor_tf_ereg_all[['eregulon_nice']] <- gsub('_', ' ', eregnames[match(cor_tf_ereg_all[['eregulon']], gsub('_', '-', eregnames))])

# plot these number
p_scenic_tf_expression_activity_correlations <- plot_grid(
  # the non TF-i-eQTLs
  ggplot(
    data = cor_tf_ereg_all[cor_tf_ereg_all[['is_scenic_plus']] & cor_tf_ereg_all[['is_tf_ieqtl']] == F & cor_tf_ereg_all[['condition']] %in% c('UT', '24hCA'), ],
    mapping = aes(
      x=condition_ordered, 
      y=eregulon_nice, 
      fill=correlation_corrected)
  ) + geom_tile() + scale_fill_gradient2(low='darkblue', mid = 'white', high='darkred', midpoint = 0)  + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
    xlab('Condition') + ylab('eRegulon') + ggtitle('Not TF-i-eQTLs') + 
    theme(legend.position = 'none'), 
  # the TF-i-eQTLs
  ggplot(
    data = cor_tf_ereg_all[cor_tf_ereg_all[['is_scenic_plus']] & cor_tf_ereg_all[['is_tf_ieqtl']] & cor_tf_ereg_all[['condition']] %in% c('UT', '24hCA'), ],
    mapping = aes(
      x=condition_ordered, 
      y=eregulon_nice, 
      fill=correlation_corrected)
  ) + geom_tile() + scale_fill_gradient2(low='darkblue', mid = 'white', high='darkred', midpoint = 0) + 
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
    xlab('Condition') + ylab('eRegulon') + ggtitle('TF-i-eQTLs') + 
    theme(legend.position = 'none'), 
  # a separate legend
  as_ggplot(
    get_legend(
      ggplot(
        data = cor_tf_ereg_all, mapping = aes(x=condition_ordered, y=eregulon, fill=correlation_corrected)
      ) + geom_tile() + scale_fill_gradient2(low='darkblue', mid = 'white', high='darkred', midpoint = 0) + labs(fill = 'Rho')
    )
  ), 
  # on one line
  nrow = 1, ncol = 3, rel_widths = c(4,4,1)
)
# show the plot
p_scenic_tf_expression_activity_correlations
# and save the result
ggsave(filename = '~/multiome/plots/mo_scenic_tf_expression_activity_correlations.pdf', plot = p_scenic_tf_expression_activity_correlations, width = 10, height = 14)

# also check the correlation difference between the conditions
p_scenic_tf_expression_activity_correlation_density <- ggplot(
  data = cor_tf_ereg_all[cor_tf_ereg_all[['is_scenic_plus']] & cor_tf_ereg_all$condition %in% c('UT', '24hCA'), ], 
  aes(x = correlation_corrected, fill = condition)
) + geom_density(alpha = .5) + scale_fill_manual(values = get_color_coding_dict()) + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  xlab('Rho') + ylab('Density') + labs(fill = 'Condition')
# show the plot
p_scenic_tf_expression_activity_correlation_density
# and save the result
ggsave(filename = '~/multiome/plots/mo_scenic_tf_expression_activity_correlation_density.pdf', plot = p_scenic_tf_expression_activity_correlation_density, width = 6, height = 4)

# check the correlations that we could get in both
cor_tf_ereg_horizontal <- merge(cor_tf_ereg_all[cor_tf_ereg_all[['condition']] == 'UT', c('eregulon', 'correlation_corrected')], cor_tf_ereg_all[cor_tf_ereg_all[['condition']] == '24hCA', c('eregulon', 'correlation_corrected')], by = 'eregulon')
# check for a difference
wilcox.test(x = cor_tf_ereg_horizontal[['correlation_corrected.x']], cor_tf_ereg_horizontal[['correlation_corrected.y']], paired = T, alternative = 'less')
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  cor_tf_ereg_horizontal[["correlation_corrected.x"]] and cor_tf_ereg_horizontal[["correlation_corrected.y"]]
# V = 25814, p-value = 1.282e-14
# alternative hypothesis: true location shift is less than 0

# make this df again, but this time keep everything
cor_tf_ereg_horizontal_all <- merge(cor_tf_ereg_all[cor_tf_ereg_all[['condition']] == 'UT', c('eregulon', 'correlation_corrected')], cor_tf_ereg_all[cor_tf_ereg_all[['condition']] == '24hCA', c('eregulon', 'correlation_corrected')], by = 'eregulon', all.x = T, all.y = T)
# set column names
colnames(cor_tf_ereg_horizontal_all) <- c('eregulon', 'rho_ut', 'rho_24hca') 
# order by name
cor_tf_ereg_horizontal_all <- cor_tf_ereg_horizontal_all[order(cor_tf_ereg_horizontal_all[['eregulon']]), ]
# save this result
cor_tf_ereg_horizontal_all_loc <- '~/multiome/tables/mo_tfe_vs_tfa_correlation_conditions.tsv.gz'
write.table(cor_tf_ereg_horizontal_all, gzfile(cor_tf_ereg_horizontal_all_loc), row.names = F, col.names = T, sep = '\t')
# with a checksum
mdfiver::create_sha256_for_file(cor_tf_ereg_horizontal_all_loc)
# and ones where the rho for either is higher than .5
cor_tf_ereg_horizontal_all_either5 <- cor_tf_ereg_horizontal_all[cor_tf_ereg_horizontal_all[['rho_ut']] > 0.5 | cor_tf_ereg_horizontal_all[['rho_24hca']] > 0.5, ]
# and where the result is in SCENIC
cor_tf_ereg_horizontal_all_either5 <- cor_tf_ereg_horizontal_all_either5[gsub('-', '_', cor_tf_ereg_horizontal_all_either5[['eregulon']]) %in% scenic_output[['Gene_signature_name']], ]
# then write that as well
cor_tf_ereg_horizontal_all_either5_loc <- '~/multiome/tables/mo_tfe_vs_tfa_correlation_conditions_cor5scenic.tsv.gz'
write.table(cor_tf_ereg_horizontal_all_either5, gzfile(cor_tf_ereg_horizontal_all_either5_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_sha256_for_file(cor_tf_ereg_horizontal_all_either5_loc)
