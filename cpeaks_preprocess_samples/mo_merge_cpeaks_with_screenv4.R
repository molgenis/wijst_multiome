#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_merge_cpeaks_with_screenv4.R
# Function: merge cpeaks annotation with the encode SCREEN v4 annotation
#
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(IRanges)
library(txtplot)
library(ggplot2)


####################
# Functions        #
####################


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


add_annotations_from_list <- function(table_to_add_to, annotation_table_list) {
  # check each of the tables
  for (table_name in names(annotation_table_list)) {
    print(table_name)
    # extract that table
    cre_annotation_table <- annotation_table_list[[table_name]]
    # make the overlap
    cre_overlap <- overlap_all_regions(cre_annotation_table, table_to_add_to,  feature_id_column2='name', feature_id_column1='name', chromosome_column2='chr_hg38', chromosome_column1='chromosome', start_column2='start_hg38', start_column1='start', end_column2='end_hg38', end_column1='end')
    # remove where there is no overlap
    cre_overlap <- cre_overlap[!is.na(cre_overlap[['overlapping_feature']]), ]
    # add the overlapping feature info
    cre_overlap[['type']] <- cre_annotation_table[match(cre_overlap[['overlapping_feature']], cre_annotation_table[['name']]), ][['type']]
    # subset to just the overlapping feature and the type
    cre_overlap <- cre_overlap[, c('name', 'type')]
    # make sure the overlap is a data.table
    cre_overlap <- data.table(cre_overlap)
    # make into a wide table
    cre_overlap_wide <- NULL
#    # check each individual type
#    for (dtype in unique(cre_overlap[['type']])) {
#      # subset to that type
#      cre_overlap_type <- cre_overlap[!is.na(cre_overlap[['type']]) & cre_overlap[['type']] == dtype, ]
#      # rename columns
#      colnames(cre_overlap_type) <- c('name', dtype)
#      # merge onto table
#      if (is.null(cre_overlap_wide)) {
#        cre_overlap_wide <- cre_overlap_type
#      }
#      else {
#        cre_overlap_wide <- merge(cre_overlap_wide, cre_overlap_type, all = T, by = 'name')
#      }
#    }
#    # add a merged table
#    cre_overlap_wide[[table_name]] <- apply(cre_overlap_wide, 1, function(x) {
#      # get the names of x
#      x_names <- names(x)
#      # remove 'name'
#      x_filtered <- x[names(x) != 'name']
#      # remove all na values
#      x_filtered <- x[!is.na(x)]
#      # paste these values all together
#      x_string <- paste(x_filtered, collapse = ';')
#    })
    # order by type, so they are always pasted the same way
    cre_overlap <- cre_overlap[order(cre_overlap[['type']]), ]
    # aggregate the different types for the same region
    cre_overlap_wide <- cre_overlap[, .(y = paste(unique(type), collapse = ";")), by = name]
    # make the name of the column specifically this type
    colnames(cre_overlap_wide) <- c('name', table_name)
    # finally add it onto the original table
    table_to_add_to <- merge(table_to_add_to, cre_overlap_wide[, c('name', ..table_name)], by = 'name', all.x = T, allow.cartesian = T)
  }
  return(table_to_add_to)
}


####################
# Main code        #
####################

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_info.tsv'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = ' ')
# read the screen annotation
screen_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/GRCh38-cCREs.bed'
screen_anno <- fread(screen_anno_loc, header = F, sep = '\t')
colnames(screen_anno) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
# add name based on the location
cpeaks_anno[['name']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
screen_anno[['name']] <- paste(screen_anno[['chromosome']], screen_anno[['start']], screen_anno[['end']], sep = '-')

# get where table 1 overlaps with table 2
cpeaks_overlaps_screen <- overlap_all_regions(screen_anno, cpeaks_anno,  feature_id_column2='name', feature_id_column1='name', chromosome_column2='chr_hg38', chromosome_column1='chromosome', start_column2='start_hg38', start_column1='start', end_column2='end_hg38', end_column1='end')
screen_overlaps_cpeaks <- overlap_all_regions(cpeaks_anno, screen_anno,  feature_id_column1='name', feature_id_column2='name', chromosome_column1='chr_hg38', chromosome_column2='chromosome', start_column1='start_hg38', start_column2='start', end_column1='end_hg38', end_column2='end')

# get the specific screen annotations now
screen_els_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/GRCh38-cCREs.ELS.bed.gz'
screen_pls_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/GRCh38-cCREs.PLS.bed.gz'
screen_ctcf_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/GRCh38-CTCF.bed.gz'
screen_els <- fread(screen_els_loc, header = F, sep = '\t')
screen_pls <- fread(screen_pls_loc, header = F, sep = '\t')
screen_ctcf <- fread(screen_ctcf_loc, header = F, sep = '\t')
# update column names
colnames(screen_els) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
colnames(screen_pls) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
colnames(screen_ctcf) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
# add name
screen_els[['name']] <- paste(screen_els[['chromosome']], screen_els[['start']], screen_els[['end']], sep = '-')
screen_pls[['name']] <- paste(screen_pls[['chromosome']], screen_pls[['start']], screen_pls[['end']], sep = '-')
screen_ctcf[['name']] <- paste(screen_ctcf[['chromosome']], screen_ctcf[['start']], screen_ctcf[['end']], sep = '-')

# get all the silencers
screen_silence_cai_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/Cai-Fullwood-2021.Silencer-cCREs.bed.gz'
screen_silence_huan_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/Huan-Ovcharenko-2019.Silencer-cCREs.bed.gz'
screen_silence_jayavelu_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/Jayavelu-Hawkins-2020.Silencer-cCREs.bed.gz'
screen_silence_rest_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/REST-Silencers.bed.gz'
screen_silence_starr_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/STARR-Silencers.Robust.bed.gz'
screen_silence_cai <- fread(screen_silence_cai_loc, header = F, sep = '\t')
screen_silence_huan <- fread(screen_silence_huan_loc, header = F, sep = '\t')
screen_silence_jayavelu <- fread(screen_silence_jayavelu_loc, header = F, sep = '\t')
screen_silence_rest <- fread(screen_silence_rest_loc, header = F, sep = '\t')
screen_silence_starr <- fread(screen_silence_starr_loc, header = F, sep = '\t')
# update screen_silence_cai names
colnames(screen_silence_cai) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
colnames(screen_silence_huan) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
colnames(screen_silence_jayavelu) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
colnames(screen_silence_rest) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
colnames(screen_silence_starr) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
# add name
screen_silence_cai[['name']] <- paste(screen_silence_cai[['chromosome']], screen_silence_cai[['start']], screen_silence_cai[['end']], sep = '-')
screen_silence_huan[['name']] <- paste(screen_silence_huan[['chromosome']], screen_silence_huan[['start']], screen_silence_huan[['end']], sep = '-')
screen_silence_jayavelu[['name']] <- paste(screen_silence_jayavelu[['chromosome']], screen_silence_jayavelu[['start']], screen_silence_jayavelu[['end']], sep = '-')
screen_silence_rest[['name']] <- paste(screen_silence_rest[['chromosome']], screen_silence_rest[['start']], screen_silence_rest[['end']], sep = '-')
screen_silence_starr[['name']] <- paste(screen_silence_starr[['chromosome']], screen_silence_starr[['start']], screen_silence_starr[['end']], sep = '-')

# get the enhancers
screen_enhancer_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/ccre_by_class/REST-Enhancers.bed.gz'
screen_enhancer <- fread(screen_enhancer_loc, header = F, sep = '\t')
colnames(screen_enhancer) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
screen_enhancer[['name']] <- paste(screen_enhancer[['chromosome']], screen_enhancer[['start']], screen_enhancer[['end']], sep = '-')

# put all of this in a list
screen_all_list <- list(
  'screen_all' = screen_anno,
  'screen_els' = screen_els,
  'screen_pls' = screen_pls,
  'screen_ctcf' = screen_ctcf, 
  'screen_enhancer' = screen_enhancer,
  'screen_sil_cai' = screen_silence_cai, 
  'screen_sil_huan' = screen_silence_huan, 
  'screen_sil_jayavelu' = screen_silence_jayavelu, 
  'screen_sil_rest' = screen_silence_rest, 
  'screen_sil_starr' = screen_silence_starr
)

# add all the annotations
cpeaks_all_annos <- add_annotations_from_list(cpeaks_anno, screen_all_list)
# where to store the result
cpeaks_annotated_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_wscreenv4.tsv.gz'
# write the result
write.table(cpeaks_all_annos, gzfile(cpeaks_annotated_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# and make a checksum
mdfiver::create_md5_for_file(cpeaks_annotated_loc)

# check what we can find back of our annotations in screen
length(unique(cpeaks_overlaps_screen[!is.na(cpeaks_overlaps_screen[['overlapping_feature']]), ][['name']])) / nrow(cpeaks_anno)
# [1] 0.7326725
length(unique(screen_overlaps_cpeaks[!is.na(screen_overlaps_cpeaks[['overlapping_feature']]), ][['name']])) / nrow(screen_anno)
# [1] 0.7617672

# get the categories that are present in the original screen data
screen_category_numbers <- data.frame(table(screen_anno[['type']]))
# and which one we can find back
screen_category_numbers_filtered <- data.frame(table(screen_anno[screen_anno[['name']] %in% unique(screen_overlaps_cpeaks[!is.na(screen_overlaps_cpeaks[['overlapping_feature']]), ][['name']]), ][['type']]))
# make sure they are ordered the same
screen_category_numbers <- screen_category_numbers[order(screen_category_numbers[['Var1']]), ]
screen_category_numbers_filtered <- screen_category_numbers_filtered[order(screen_category_numbers_filtered[['Var1']]), ]
# create new table that has the numbers for what was found back or not
screen_found_not_found_numbers <- data.frame(
  'category' = rep(screen_category_numbers[['Var1']], times = 2), 
  'found' = c(rep('both', times = nrow(screen_category_numbers)), rep('screen only', times = nrow(screen_category_numbers))), 
  'number' = c(screen_category_numbers_filtered[['Freq']], (screen_category_numbers[['Freq']] - screen_category_numbers_filtered[['Freq']])), 
  'frac' = c((screen_category_numbers_filtered[['Freq']] / screen_category_numbers[['Freq']]), ((screen_category_numbers[['Freq']] - screen_category_numbers_filtered[['Freq']]) / screen_category_numbers[['Freq']]))
)
# add category + found
screen_found_not_found_numbers[['category_found']] <- paste(screen_found_not_found_numbers[['category']], screen_found_not_found_numbers[['found']])
# actual plot now
p_cpeaks_in_screen_numbers <- ggplot(data  = screen_found_not_found_numbers, mapping  = aes(x = category, y = number, fill = found)) + 
  # barplot
  geom_bar(position = 'stack', stat = 'identity') +
  # with manual colors
  scale_fill_manual(values = list('screen only' = 'darkblue', 'both' = 'purple')) +
  # labels
  xlab('Element category') + 
  ylab('Number of Elements') + 
  ggtitle('Categories of elements from ENCODE SCREEN v4\nfound in cPeaks reference') + 
  # add more whitespace
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# do the fractions as well
p_cpeaks_in_screen_fracs <- ggplot(data  = screen_found_not_found_numbers, mapping  = aes(x = category, y = frac, fill = found)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = list('screen only' = 'darkblue', 'both' = 'purple')) +
  xlab('Element category') + 
  ylab('Fraction of Elements') + 
  ggtitle('Categories of elements from ENCODE SCREEN v4\nfound in cPeaks reference') + 
  theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
  # rotate the x axis ticks
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# show the plots
p_cpeaks_in_screen_numbers
p_cpeaks_in_screen_fracs

# get the categories that are present in the original screen data
screen_category_numbers <- data.frame(table(screen_anno[['type']]))
# and which one we can find back
screen_category_numbers_filtered <- data.frame(table(screen_anno[screen_anno[['name']] %in% unique(screen_overlaps_cpeaks[!is.na(screen_overlaps_cpeaks[['overlapping_feature']]), ][['name']]), ][['type']]))
# make sure they are ordered the same
screen_category_numbers <- screen_category_numbers[order(screen_category_numbers[['Var1']]), ]
screen_category_numbers_filtered <- screen_category_numbers_filtered[order(screen_category_numbers_filtered[['Var1']]), ]
# create new table that has the numbers for what was found back or not
screen_found_not_found_numbers <- data.frame(
  'category' = rep(screen_category_numbers[['Var1']], times = 2), 
  'found' = c(rep('both', times = nrow(screen_category_numbers)), rep('screen only', times = nrow(screen_category_numbers))), 
  'number' = c(screen_category_numbers_filtered[['Freq']], (screen_category_numbers[['Freq']] - screen_category_numbers_filtered[['Freq']])), 
  'frac' = c((screen_category_numbers_filtered[['Freq']] / screen_category_numbers[['Freq']]), ((screen_category_numbers[['Freq']] - screen_category_numbers_filtered[['Freq']]) / screen_category_numbers[['Freq']]))
)

