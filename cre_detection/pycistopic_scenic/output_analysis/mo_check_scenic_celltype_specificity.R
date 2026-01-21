#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_scenic_celltype_specificity.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(ggplot2)


####################
# Functions        #
####################


rss_table_to_long <- function(rss_table) {
  # get the categories
  categories <- rss_table[[1]]
  # get the regulons
  regulons <- colnames(rss_table)[-1]
  # create the new table
  long_table <- data.frame('regulon' = regulons)
  # next go through the categories
  for (category in categories) {
    # add that to the table
    long_table[[category]] <- as.vector(unlist(rss_table[rss_table[[1]] == category, c(-1)]))
  }
  return(long_table)
}


#' get the eGenes per cell type from QTL output
#' 
#' @param qtl_output_loc base location of the QTL output per cell type
#' @param output_file which output file to read for the results
#' @returns a list with the egenes per cell type
#' 
get_qtls_per_celltype_limix <- function(qtl_output_loc, output_file='/inflammation_final/iqtl_results_all_eigenmt_qval.tsv.gz', significance_cutoffs=list('feature_q_value' = 0.05, 'total_bf_eigen' = 0.05), verbose=T) {
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
    if (!is.null(significance_cutoffs)) {
      # print progress if requested
      if (verbose) {
        print(paste('variant+phenotype before filtering', nrow(cell_type_output)))
      }
      # check each cutoff
      for(significance_column in names(significance_cutoffs)) {
        # get the cutoff
        significance_cutoff <- significance_cutoffs[[significance_column]]
        # filter
        cell_type_output <- cell_type_output[
          !is.na(cell_type_output[[significance_column]]) &
            cell_type_output[[significance_column]] < significance_cutoff, 
        ]
      }
      if (verbose) {
        print(paste('variant+phenotype after filtering', nrow(cell_type_output)))
      }
    }
    # add to the list
    qtls_per_celltype[[cell_type]] <- cell_type_output
  }
  # turn into a dataframe
  return(qtls_per_celltype)
}


get_de_outputs <- function(mast_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', lfc_cutoff=NULL, to_ens=F, symbols.to.ensg.mapping='genes.tsv', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w'), de_method='mast_seurat', condition_name='timepoint'){
  # set up per cell type
  de_per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # set up per stim combination
    de_per_condition <- list()
    # check each stim
    for(stim in stims){
      for(stim2 in stims){
        try({
          if(stim != stim2){
            print(paste(cell_type, stim, stim2, sep = ' '))
            # paste the filepath together
            filepath <- NULL
            if ("limma" == de_method) {
              #filepath <- paste(mast_output_loc, cell_type, '_', condition_name,  stim, '_', stim2, '.tsv.gz', sep = '')
              filepath <- paste(mast_output_loc, cell_type, '_', condition_name, '.tsv.gz', sep = '')
            }
            else {
              filepath <- paste(mast_output_loc, cell_type, stim,stim2, '.tsv.gz', sep = '')
            }
            if (file.exists(filepath)) {
              # read the mast output
              mast <- read.table(filepath, header=T)
              # filter to only include the significant results
              mast <- mast[mast[[pval_column]] <= 0.05, ]
              # filter for only the positive lfc if required
              if(only_positive){
                mast <- mast[mast[[lfc_column]] < 0, ]
              }
              # filter for only the positive lfc if required
              if(only_negative){
                mast <- mast[mast[[lfc_column]] > 0, ]
              }
              # filter only ones with strong enough effect if required
              if (!is.null(lfc_cutoff)) {
                mast <- mast[abs(mast[[lfc_column]]) >= lfc_cutoff, ]
              }
              # confine in some way if reporting a max number of genes
              if(!is.null(max)){
                # by p if required
                if(max_by_pval){
                  mast <- mast[order(mast[[p_val_column]]), ]
                }
                # by lfc otherwise
                else{
                  mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
                }
                # subset to the number we requested if max was set
                mast <- mast[1:max,]
              }
              # grab the genes from the column names
              genes <- rownames(mast)
              if ("limma" == de_method) {
                genes <- mast[['feature']]
              }
              # convert the symbols to ensemble IDs
              if (to_ens) {
                mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
                mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
                genes <- mapping[match(genes, mapping$V2),"V1"]
              }
              # otherwise change the Seurat replacement back
              else{
                #genes <- gsub("-", "_", genes)
              }
              # change gene names
              if ("limma" == de_method) {
                mast[['feature']] <- genes
              }
              else {
                rownames(mast) <- genes
              }
              de_per_condition[[paste(stim, stim2, sep = '')]] <- mast
            }
          }
        })
      }
    }
    de_per_ct[[cell_type]] <- de_per_condition
  }
  return(de_per_ct)
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
  color_coding_dict[['monocyte']] <- '#EDBA1B'
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
  color_coding_dict[['Stem cells']] <- '#66A61E'
  color_coding_dict[['Stromal cells']] <- '#8DD3C7'
  #color_coding_dict[['other T cells']] <- '#A6761D'
  color_coding_dict[['other T cells']] <- '#FF63B6'
  color_coding_dict[['Transit amplifying cells']] <- '#FF7F00'
  color_coding_dict[['Transit\namplifying cells']] <- '#FF7F00'
  color_coding_dict[['disconcordant']] <- 'gray'
  color_coding_dict[['lymphoid']] <- '#5B7573'
  color_coding_dict[['myeloid']] <- '#C28C72'
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
  color_coding_dict[['UT']] <- 'lightgrey'
  color_coding_dict[['24hCA']] <- 'forestgreen'
  return(color_coding_dict)
}



####################
# Settings         #
####################

# whether we are in debug mode
debug <- F


####################
# Main code        #
####################

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# get the enrichment statistics
rss_condition_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_condition.tsv.gz'
rss_cell_type_major_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_cell_type_lowerres.tsv.gz'
rss_cell_type_major_merged_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_cell_type_lowerres_merged.tsv.gz'
rss_cell_type_minor_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_cell_type_highres.tsv.gz'
rss_cell_type_lineage_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/tf_to_metadata/mo_tf_rss_cell_type_lineage.tsv.gz'

# read the scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# and enrichments
rss_condition <- fread(rss_condition_loc, header = T, sep = '\t')
rss_cell_type_lineage <- fread(rss_cell_type_lineage_loc, header = T, sep = '\t')
rss_cell_type_major_merged <- fread(rss_cell_type_major_merged_loc, header = T, sep = '\t')
# convert rss table to a long table
rss_condition_long <- rss_table_to_long(rss_condition)
rss_cell_type_lineage_long <- rss_table_to_long(rss_cell_type_lineage)
rss_cell_type_major_merged_long <- rss_table_to_long(rss_cell_type_major_merged)


# also make a plottble table
rss_condition_plottable <- rbind(
  data.frame('rss' = rss_condition_long[['UT']], 'condition' = rep('UT', times = nrow(rss_condition_long))), 
  data.frame('rss' = rss_condition_long[['24hCA']], 'condition' = rep('24hCA', times = nrow(rss_condition_long)))
)
rss_cell_type_lineage_plottable <- rbind(
  data.frame('rss' = rss_cell_type_lineage_long[['myeloid']], 'condition' = rep('myeloid', times = nrow(rss_cell_type_lineage_long))), 
  data.frame('rss' = rss_cell_type_lineage_long[['lymphoid']], 'condition' = rep('lymphoid', times = nrow(rss_cell_type_lineage_long)))
)

# and plot that
ggplot(data = rss_cell_type_lineage_plottable, mapping = aes(x = rss, fill = condition)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = get_color_coding_dict()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('Myeloid and Lymphoid enrichment scores') +
  xlab('RSS eRegulon enrichment score') + 
  ylab('Density') +
  labs(fill = 'Lineage')

# calculate the difference between the two lineages
rss_cell_type_lineage_long[['delta']] <- abs(rss_cell_type_lineage_long[['myeloid']] - rss_cell_type_lineage_long[['lymphoid']])
rss_cell_type_lineage_long[['diff_mye_lym']] <- (rss_cell_type_lineage_long[['myeloid']] - rss_cell_type_lineage_long[['lymphoid']])
# calculate the difference between the two conditions
rss_condition_long[['delta']] <- abs(rss_condition_long[['24hCA']] - rss_condition_long[['UT']])
rss_condition_long[['diff_24hCA_UT']] <- (rss_condition_long[['24hCA']] - rss_condition_long[['UT']])

# plot the delta as wel
ggplot(data = rss_cell_type_lineage_long, mapping = aes(x = delta)) +
  geom_density(alpha = 0.5, fill = 'lightgray') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('Myeloid and Lymphoid enrichment score deltas') +
  xlab('RSS eRegulon enrichment score') + 
  ylab('Density') +
  labs(fill = 'Lineage Delta RSS')
# plot the non-absolute difference
ggplot(data = rss_cell_type_lineage_long, mapping = aes(x = diff_mye_lym)) +
  geom_density(alpha = 0.5, fill = 'lightgray') +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('Myeloid and Lymphoid enrichment score deltas') +
  xlab('RSS eRegulon enrichment score') + 
  ylab('Density') +
  labs(fill = 'Lineage difference RSS')

# merge the linage specificity with the stim specificity
rss_diff_cond_to_lin <- merge(rss_condition_long, rss_cell_type_lineage_long, by = 'regulon')
# now plot these against one another
ggplot(data = rss_diff_cond_to_lin, mapping = aes(x = diff_24hCA_UT, y = diff_mye_lym)) +
  geom_point() +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('Myeloid and Lymphoid enrichment diff vs 24hCA and UT enrichment diff') +
  ylab('RSS eRegulon enrichment score diff lineage') + 
  xlab('RSS eRegulon enrichment score diff condition') +
  geom_smooth(method='lm', formula= y~x)

# merge lineage with cell type
rss_ct_vs_lin <- merge(rss_cell_type_major_merged_long, rss_cell_type_lineage_long, by = 'regulon')
# check the major cell types in the lineages
ggplot(data = rss_ct_vs_lin, mapping = aes(x = monocyte, y = CD4T)) +
  geom_point() +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('Myeloid vs lineage RSS enrichment scores') +
  ylab('RSS eRegulon enrichment CD4+ T') + 
  xlab('RSS eRegulon enrichment monocyte') +
  geom_smooth(method='lm', formula= y~x)
# and the lineages themselves, to see if they are similiar
ggplot(data = rss_ct_vs_lin, mapping = aes(x = myeloid, y = lymphoid)) +
  geom_point() +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('Myeloid vs lineage RSS enrichment scores') +
  ylab('RSS eRegulon enrichment lymphoid') + 
  xlab('RSS eRegulon enrichment myeloid') +
  geom_smooth(method='lm', formula= y~x)

# make a plottable version of the major cell type RSS distributions
rss_cell_type_major_merged_plottable_list <- list()
# by checking each of the columns
for (ct in setdiff(colnames(rss_cell_type_major_merged_long), c('regulon'))) {
  # make a df for the regulon and the cell type
  df_ct <- data.table('regulon' = rss_cell_type_major_merged_long[['regulon']], 'cell_type' = rep(ct, times = nrow(rss_cell_type_major_merged_long)), 'rss' = rss_cell_type_major_merged_long[[ct]])
  # put that in the list
  rss_cell_type_major_merged_plottable_list[[ct]] <- df_ct
}
# and put them all together in the end
rss_cell_type_major_merged_plottable <- do.call('rbind', rss_cell_type_major_merged_plottable_list)
# plot those distributions
ggplot(data = rss_cell_type_major_merged_plottable, mapping = aes(x = rss, fill = cell_type)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = get_color_coding_dict()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('Major cell type enrichment scores') +
  xlab('RSS eRegulon enrichment score') + 
  ylab('Density') +
  labs(fill = 'Cell type')
# also plot them as bars
ggplot(data = rss_cell_type_major_merged_plottable, mapping = aes(x = regulon, fill = cell_type, y = rss)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = get_color_coding_dict()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('eRegulon enrichment scores') +
  xlab('eRegulon') + 
  ylab('RSS enrichment score') +
  labs(fill = 'Cell type') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# add the regulons to the plottable lineages
rss_cell_type_lineage_plottable <- cbind(rbind(rss_cell_type_lineage_long[, c('regulon'), drop = F], rss_cell_type_lineage_long[, c('regulon'), drop = F]), rss_cell_type_lineage_plottable)
# and plot bar of lineage scores
ggplot(data = rss_cell_type_lineage_plottable, mapping = aes(x = regulon, fill = condition, y = rss)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = get_color_coding_dict()) +
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  ggtitle('eRegulon enrichment scores') +
  xlab('eRegulon') + 
  ylab('RSS enrichment score') +
  labs(fill = 'Lineage') +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

