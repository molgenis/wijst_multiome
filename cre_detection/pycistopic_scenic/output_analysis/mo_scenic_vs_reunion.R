#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_scenic_vs_reunion.R
# Function: compare SCENIC+ output to REUNION output
############################################################################################################################


####################
# libraries        #
####################

library(data.table)
library(ggplot2)
library(ggvenn)
library(cowplot)
library(stringr)
library(IRanges)


####################
# Functions        #
####################


chromatin_feature_to_positions <- function(feature_names, split_regex='-|:') {
  # make sure to have no duplicates
  feature_names <- unique(feature_names)
  # create table of locations
  feature_locations_table <- data.frame(
    'feature_id' = feature_names, 
    'chromosome' = rep(NA, times = length(feature_names)), 
    'start' = rep(NA, times = length(feature_names)), 
    'end' = rep(NA, times = length(feature_names))
  )
  # now split the features to get the locations
  feature_locations_table[c('chromosome', 'start', 'end')] <- str_split_fixed(feature_names, split_regex, 3)
  # and make the start and end numeric
  feature_locations_table[['start']] <- as.numeric(feature_locations_table[['start']])
  feature_locations_table[['end']] <- as.numeric(feature_locations_table[['end']])
  return(feature_locations_table)
}


overlap_all_regions <- function(region_table1, region_table2, feature_id_column1='feature_id', feature_id_column2='feature_id', chromosome_column1='chromosome', chromosome_column2='chromosome', start_column1='start', start_column2='start', end_column1='end', end_column2='end'){
  # make chromsomes characters
  region_table1[[chromosome_column1]] <- as.character(region_table1[[chromosome_column1]])
  region_table2[[chromosome_column2]] <- as.character(region_table2[[chromosome_column2]])
  # make data framke
  region_table1 <- data.frame(region_table1)
  region_table2 <- data.frame(region_table2)
  # get the unique chromosomes
  chromosomes <- intersect(unique(region_table1[[chromosome_column1]]), unique(region_table2[[chromosome_column2]]))
  # we'll read the LCL output per chromosome, so we'll store that in a list first
  overlap_per_chromosome <- list()
  # we'll check each chromosome
  for (chromosome in chromosomes) {
    # subset to this chromosome
    table1_regions_chromosome <- region_table1[!is.na(region_table1[[chromosome_column1]]) & region_table1[[chromosome_column1]] == chromosome, ]
    table2_regions_chromosome <- region_table2[!is.na(region_table2[[chromosome_column2]]) & region_table2[[chromosome_column2]] == chromosome, , ]
    # turn into iranges objects
    table1_chromosome_iranges <- IRanges(start = table1_regions_chromosome[[start_column1]], end = table1_regions_chromosome[[end_column1]])
    table2_chromosome_iranges <- IRanges(start = table2_regions_chromosome[[start_column2]], end = table2_regions_chromosome[[end_column2]])
    # find overlaps
    feature_chromosome_overlaps <- findOverlaps(table1_chromosome_iranges, table2_chromosome_iranges)
    # extract overlapping ranges
    overlapping_ranges <- pintersect(table2_chromosome_iranges[subjectHits(feature_chromosome_overlaps)], table1_chromosome_iranges[queryHits(feature_chromosome_overlaps)])
    # create a  table for the overlaps
    overlaps_table <- data.table(
      'feature_id' = table2_regions_chromosome[[feature_id_column2]][subjectHits(feature_chromosome_overlaps)],
      'overlapping_feature' = table1_regions_chromosome[[feature_id_column1]][queryHits(feature_chromosome_overlaps)],
      'overlap_start' = start(overlapping_ranges),
      'overlap_end' = end(overlapping_ranges)
    )
    # add these positions
    table2_regions_chromosome <- merge(table2_regions_chromosome, overlaps_table, by.x = feature_id_column2, by.y = 'feature_id', all.x = T, allow.cartesian=TRUE)
    # and put in the list
    overlap_per_chromosome[[as.character(chromosome)]] <- table2_regions_chromosome
  }
  overlap_all <- do.call('rbind', overlap_per_chromosome)
  return(overlap_all)
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


###################
# Settings        #
###################


####################
# Main Code        #
####################

# location of the REUNION output
reunion_output_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/reunion_yang2024_cres/df_pbmc_region_tf_gene_link.2_2.txt.gz'

# location of our scenic+ output
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# and cpeaks annotation we used
cpeaks_annotation_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/cPeaks/cPeaks_info.tsv'

# read the REUNION table
reunion_output <- fread(reunion_output_loc, header = T, sep = '\t')
# get the scenic+ table
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')

# keep only +/-
scenic_output <- scenic_output[scenic_output[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# order by the extended
scenic_output <- scenic_output[order(scenic_output[['is_extended']]), ]
# get which are non-extended if it was both extended and non-extended
scenic_eregs <- unique(scenic_output[, c('TF', 'Gene_signature_direction', 'Gene_signature_name', 'source')])
scenic_eregs <-scenic_eregs[order(scenic_eregs[['source']]), ]
scenic_eregs_to_keep <- scenic_eregs[!duplicated(paste(scenic_eregs[['TF']], scenic_eregs[['Gene_signature_direction']])), ]
# then use that to filer
scenic_output <- scenic_output[scenic_output[['Gene_signature_name']] %in% scenic_eregs_to_keep[['Gene_signature_name']], ]


# get cpeaks
cpeaks_annotation <- fread(cpeaks_annotation_loc, header = T, sep = ' ')
# add scenic style annotation
cpeaks_annotation[['scenic_hg38']] <- paste0(cpeaks_annotation[['chr_hg38']], ':', cpeaks_annotation[['start_hg38']], '-', cpeaks_annotation[['end_hg38']])

# read gene annotations
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/genome_annotation.tsv'
gene_anno <- fread(gene_anno_loc, header = T, sep = '\t')
# rename columns to be the same as in limix
colnames(gene_anno) <-c('chrom', 'start', 'end', 'strand', 'gs','Transcription_Start_Site','Transcript_type')

# get the positions in reunion
reunion_feature_locations <- data.frame(str_split_fixed(reunion_output[['peak_id']], ':|-', 3))
# add the peak name itself
reunion_feature_locations <- cbind(data.frame(x = reunion_output[['peak_id']]), reunion_feature_locations)
# and set column names
colnames(reunion_feature_locations) <- c('feature_id', 'chromosome', 'start', 'end')
# and remove any duplicates
reunion_feature_locations <- unique(reunion_feature_locations)
# and make numeric for the start and end
reunion_feature_locations[['start']] <- as.numeric(reunion_feature_locations[['start']])
reunion_feature_locations[['end']] <- as.numeric(reunion_feature_locations[['end']])

# subset the search space of scenic, so cpeaks, to the columns we need
scenic_feature_locations <- cpeaks_annotation[, c('chr_hg38', 'start_hg38', 'end_hg38')]
# add an id as well here
scenic_feature_locations <- cbind(data.frame(x = paste0(scenic_feature_locations[['chr_hg38']], ':', scenic_feature_locations[['start_hg38']], '-', scenic_feature_locations[['end_hg38']])), scenic_feature_locations)
# and rename columns again
colnames(scenic_feature_locations) <- c('feature_id', 'chromosome', 'start', 'end')

# get matched regions
reunion_to_mo_regions <- overlap_all_regions(scenic_feature_locations, reunion_feature_locations)
# calculate the overlap size
reunion_to_mo_regions[['overlap_size']] <- abs(reunion_to_mo_regions[['overlap_end']] - reunion_to_mo_regions[['overlap_start']])
# set overlap size to zero for anything without an overlap
reunion_to_mo_regions[is.na(reunion_to_mo_regions[['overlap_size']]), ][['overlap_size']] <- 0

# now add the overlap information
reunion_output <- merge(reunion_output, data.table(reunion_to_mo_regions[, c('feature_id', 'overlapping_feature', 'overlap_size')]), all.x = T, by.x = 'peak_id', by.y = 'feature_id', allow.cartesian=TRUE)

# add location for SCENIC+
scenic_output <- cbind(scenic_output, cpeaks_annotation[match(scenic_output[['Region']], cpeaks_annotation[['scenic_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# and the locations of the genes
scenic_output <- cbind(scenic_output, gene_anno[match(scenic_output[['Gene']], gene_anno[['gs']]), c('chrom', 'start', 'end')])
# get the distances again
scenic_distances <- get_closest_flanks(scenic_output, 'start_hg38', 'end_hg38', 'start', 'end')
# add that to the original table
scenic_output[['distance']] <- scenic_distances[['min_dist']]
# keep only +/-
scenic_output <- scenic_output[scenic_output[['Gene_signature_direction']] %in% c('+/+', '-/+'), ]
# order by the extended
scenic_output <- scenic_output[order(scenic_output[['is_extended']]), ]
# get which are non-extended if it was both extended and non-extended
scenic_eregs <- unique(scenic_output[, c('TF', 'Gene_signature_direction', 'Gene_signature_name', 'source')])
scenic_eregs <-scenic_eregs[order(scenic_eregs[['source']]), ]
scenic_eregs_to_keep <- scenic_eregs[!duplicated(paste(scenic_eregs[['TF']], scenic_eregs[['Gene_signature_direction']])), ]
# then use that to filer
scenic_output <- scenic_output[scenic_output[['Gene_signature_name']] %in% scenic_eregs_to_keep[['Gene_signature_name']], ]

# add location for SCENIC+
reunion_output <- cbind(reunion_output, cpeaks_annotation[match(reunion_output[['overlapping_feature']], cpeaks_annotation[['scenic_hg38']]), c('chr_hg38', 'start_hg38', 'end_hg38')])
# and the locations of the genes
reunion_output <- cbind(reunion_output, gene_anno[match(reunion_output[['gene_id']], gene_anno[['gs']]), c('chrom', 'start', 'end')])
# get the distances again
reunion_distances <- get_closest_flanks(reunion_output, 'start_hg38', 'end_hg38', 'start', 'end')
# add that to the original table
reunion_output[['distance']] <- reunion_distances[['min_dist']]
# and region-gene overlaps
reunion_output_unfiltered <- reunion_output
reunion_output <- reunion_output[reunion_output[['distance']] > 0 & reunion_output[['distance']] <= 150000, ]


# now check how it looks if you add these
ggvenn::ggvenn(
  data = list('multiome' = unique(scenic_output[['Region']]), 'REUNION' = unique(reunion_output[['overlapping_feature']]))
) + ggtitle('Overlap of CREs between multiome and REUNION')
ggvenn::ggvenn(
  data = list('multiome' = unique(scenic_output[['Gene']]), 'REUNION' = unique(reunion_output[['gene_id']]))
) + ggtitle('Overlap of CRE genes between multiome and REUNION')
ggvenn::ggvenn(
  data = list('multiome' = unique(paste(scenic_output[['Region']], scenic_output[['Gene']])), 'REUNION' = unique(paste(reunion_output[['overlapping_feature']], reunion_output[['gene_id']])))
) + ggtitle('Overlap of region-gene links between multiome and REUNION')
ggvenn::ggvenn(
  data = list('multiome' = unique(paste(scenic_output[['Region']], scenic_output[['Gene']], scenic_output[['TF']])), 'REUNION' = unique(paste(reunion_output[['overlapping_feature']], reunion_output[['gene_id']], reunion_output[['motif_id']])))
) + ggtitle('Overlap of triplets between multiome and REUNION')

# write the unfiltered reunion output as well
write.table(reunion_output_unfiltered, gzfile('/groups/umcg-franke-scrna/tmp02/external_datasets/reunion_yang2024_cres/df_pbmc_region_tf_gene_link.2_2_cpeaksmatched.tsv.gz'), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_md5_for_file('/groups/umcg-franke-scrna/tmp02/external_datasets/reunion_yang2024_cres/df_pbmc_region_tf_gene_link.2_2_cpeaksmatched.tsv.gz')

# add columns that are region-to-gene
reunion_output[['r2g2tf']] <- paste(reunion_output[['overlapping_feature']], reunion_output[['gene_id']], reunion_output[['motif_id']])
scenic_output[['r2g2tf']] <- paste(scenic_output[['Region']], scenic_output[['Gene']], scenic_output[['TF']])
# merge these
matched_output <- merge(reunion_output[, c('r2g2tf', 'peak_gene_corr', 'gene_tf_corr')], scenic_output[, c('r2g2tf', 'rho_R2G', 'rho_TF2G')], by = 'r2g2tf')
# rename columns
colnames(matched_output) <- c('r2g2tf', 'cor_reunion', 'tfcor_reunion', 'cor_scenic', 'tfcor_scenic')
# and make unique
matched_output <- unique(matched_output)

# calculate the concordance
mo_reunion_rho_concordance <- sum(sign(matched_output[['cor_scenic']]) == sign(matched_output[['cor_reunion']])) / nrow(matched_output)
# get the minimal and maximum correlations
min_sig_rho_10x <- min(abs(matched_output[['cor_reunion']]))
min_sig_rho_mo <- min(abs(matched_output[['cor_scenic']]))
max_sig_rho_10x <- max(abs(matched_output[['cor_reunion']]))
max_sig_rho_mo <- max(abs(matched_output[['cor_scenic']]))
# plot the correlations
p_mo_vs_reunion <- ggplot(data = matched_output, mapping = aes(x = cor_scenic, y = cor_reunion)) + 
  geom_point(size = .1) +
  xlab('R2G Rho in multiome') + 
  ylab('R2G Rho in reunion') + 
  ggtitle('R2G correlations in multiome vs reunion (matched TF only)') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # left to right block of non-significant effects
  geom_rect(aes(xmin = -1 * max_sig_rho_mo, xmax = max_sig_rho_mo, ymin = -1 * min_sig_rho_10x, ymax = min_sig_rho_10x), 
            fill = "white", alpha = 0.005) +
  # bottom to top block of non-significant effects
  geom_rect(aes(xmin = -1 * min_sig_rho_mo, xmax = min_sig_rho_mo, ymin = -1 * max_sig_rho_10x, ymax = max_sig_rho_10x), 
            fill = "white", alpha = 0.005) +
  # bottom left block
  geom_rect(aes(xmin = -1 * max_sig_rho_mo, xmax = -1 *min_sig_rho_mo, ymin = -1 * max_sig_rho_10x, ymax = -1 * min_sig_rho_10x), 
            fill = "#0072B2", alpha = 0.01) +
  # bottom right block
  geom_rect(aes(xmin = min_sig_rho_mo, xmax = max_sig_rho_mo, ymin = -1 * max_sig_rho_10x, ymax = -1 * min_sig_rho_10x), 
            fill = "#D55E00", alpha = 0.01) +
  # top left block
  geom_rect(aes(xmin = -1 * max_sig_rho_mo, xmax = -1 *min_sig_rho_mo, ymin = min_sig_rho_10x, ymax = max_sig_rho_10x), 
            fill = "#D55E00", alpha = 0.01) + 
  # top right block
  geom_rect(aes(xmin = min_sig_rho_mo, xmax = max_sig_rho_mo, ymin = max_sig_rho_10x, ymax = min_sig_rho_10x), 
            fill = "#0072B2", alpha = 0.01) + 
  # add the concordance
  annotate("label", x = max_sig_rho_mo * 0.75 , y = max_sig_rho_10x * -0.75, label = paste('concordance', round(mo_reunion_rho_concordance, digits = 2), sep = ':\n')) +
  # add the names of the concordant and non-concordant blocks
  annotate("text", x = max_sig_rho_mo * -0.70 , y = max_sig_rho_10x * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold') +
  # add the names of the concordant and non-concordant blocks
  annotate("text", x = max_sig_rho_mo * 0.70 , y = max_sig_rho_10x * 0.75, label = 'concordant', colour = '#0072B2', fontface = 'bold')

# show the plot
p_mo_vs_reunion
# save the plot
ggsave(filename = '~/multiome/plots/mo_multiome_vs_reunion_r2g.pdf', plot = p_mo_vs_reunion, width = 5, height = 5)


# calculate the concordance of TF
mo_reunion_rhotf_concordance <- sum(sign(matched_output[['tfcor_scenic']]) == sign(matched_output[['tfcor_reunion']])) / nrow(matched_output)
# get the minimal and maximum correlations
min_sig_rhotf_10x <- min(abs(matched_output[['tfcor_reunion']]))
min_sig_rhotf_mo <- min(abs(matched_output[['tfcor_scenic']]))
max_sig_rhotf_10x <- max(abs(matched_output[['tfcor_reunion']]))
max_sig_rhotf_mo <- max(abs(matched_output[['tfcor_scenic']]))
# plot the correlations
p_mo_vs_reunion_tf <- ggplot(data = matched_output, mapping = aes(x = tfcor_scenic, y = tfcor_reunion)) + 
  geom_point(size = .1) +
  xlab('TF2G Rho in multiome') + 
  ylab('TF2G Rho in reunion') + 
  ggtitle('TF2G correlations in multiome vs reunion (matched TF only)') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # left to right block of non-significant effects
  geom_rect(aes(xmin = -1 * max_sig_rhotf_mo, xmax = max_sig_rhotf_mo, ymin = -1 * min_sig_rhotf_10x, ymax = min_sig_rhotf_10x), 
            fill = "white", alpha = 0.005) +
  # bottom to top block of non-significant effects
  geom_rect(aes(xmin = -1 * min_sig_rhotf_mo, xmax = min_sig_rhotf_mo, ymin = -1 * max_sig_rhotf_10x, ymax = max_sig_rhotf_10x), 
            fill = "white", alpha = 0.005) +
  # bottom left block
  geom_rect(aes(xmin = -1 * max_sig_rhotf_mo, xmax = -1 *min_sig_rhotf_mo, ymin = -1 * max_sig_rhotf_10x, ymax = -1 * min_sig_rhotf_10x), 
            fill = "#0072B2", alpha = 0.01) +
  # bottom right block
  geom_rect(aes(xmin = min_sig_rhotf_mo, xmax = max_sig_rhotf_mo, ymin = -1 * max_sig_rhotf_10x, ymax = -1 * min_sig_rhotf_10x), 
            fill = "#D55E00", alpha = 0.01) +
  # top left block
  geom_rect(aes(xmin = -1 * max_sig_rhotf_mo, xmax = -1 *min_sig_rhotf_mo, ymin = min_sig_rhotf_10x, ymax = max_sig_rhotf_10x), 
            fill = "#D55E00", alpha = 0.01) + 
  # top right block
  geom_rect(aes(xmin = min_sig_rhotf_mo, xmax = max_sig_rhotf_mo, ymin = max_sig_rhotf_10x, ymax = min_sig_rhotf_10x), 
            fill = "#0072B2", alpha = 0.01) + 
  # add the concordance
  annotate("label", x = max_sig_rhotf_mo * 0.75 , y = max_sig_rhotf_10x * -0.75, label = paste('concordance', round(mo_reunion_rhotf_concordance, digits = 2), sep = ':\n')) +
  # add the names of the concordant and non-concordant blocks
  annotate("text", x = max_sig_rhotf_mo * -0.70 , y = max_sig_rhotf_10x * 0.75, label = 'discordant', colour = '#D55E00', fontface = 'bold') +
  # add the names of the concordant and non-concordant blocks
  annotate("text", x = max_sig_rhotf_mo * 0.70 , y = max_sig_rhotf_10x * 0.75, label = 'concordant', colour = '#0072B2', fontface = 'bold')

# show the plot
p_mo_vs_reunion_tf
# save the plot
ggsave(filename = '~/multiome/plots/mo_multiome_vs_reunion_tf2g.pdf', plot = p_mo_vs_reunion_tf, width = 5, height = 5)
