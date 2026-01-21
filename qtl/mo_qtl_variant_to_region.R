#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_qtl_variant_to_region.R
# Function: get list of variants ever tested for QTLs
############################################################################################################################

####################
# libraries        #
####################

library(mdfiver)
library(data.table)
library(stringr)
library(IRanges)


####################
# Functions        #
####################

read_qtl_files_all <- function(unfiltered_loc, folders=NULL, unfiltered_file='qtl_results_all.txt.gz', variant_id_column='snp_id', variant_chromosome_column='snp_chromosome', variant_position_column='snp_position', verbose = T) {
  # get the folders in the directory, which should be the cell types
  cell_types <- list.dirs(unfiltered_loc, full.names = F, recursive = F)
  # check if overlaps with the folders that we want to take a look at
  if (!is.null(folders)) {
    cell_types <- intersect(cell_types, folders)
  }
  # initialize our final table
  all_variants <- NULL
  # check each cell type
  for (cell_type in cell_types) {
    # paste together the full path
    full_cell_type_path <- paste(unfiltered_loc, '/', cell_type, '/', unfiltered_file, sep = '')
    # check if the file exists
    if (file.exists(full_cell_type_path)) {
      # log if requested
      if (verbose) {
        print(paste('reading', full_cell_type_path))
      }
      # read the file
      cell_type_output <- fread(full_cell_type_path, sep = '\t', header = T)
      # subset to the variants
      cell_type_output <- cell_type_output[, c(..variant_id_column, ..variant_chromosome_column, ..variant_position_column)]
      # get only unique entries
      cell_type_output <- unique(cell_type_output)
      # and add to existing table
      if (!is.null(all_variants)) {
        all_variants <- rbind(all_variants, cell_type_output)
        # don't forget to make unique again
        all_variants <- unique(all_variants)
      }
      # if the table doesn't exist, we need to initialize it
      else {
        all_variants <- cell_type_output
      }
    }
    # warn
    else {
      warning(paste('directory for file exists, but file is not there, so skipped:', full_cell_type_path))
    }
  }
  # order the table
  all_variants <- all_variants[order(all_variants[[variant_chromosome_column]], all_variants[[variant_position_column]], all_variants[[variant_id_column]]), ]
  # and return the result
  return(all_variants)
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


####################
# Main Code        #
####################


###################
# mo eQTLs        #
###################

# location of the QTL outputs
eqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/output/L1/combined/'
# location of the caQTL
caqtl_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/output/L1/combined/'

# get all eQTL variants
eqtl_variants_all <- read_qtl_files_all(eqtl_output_loc)
caqtl_output_all <- read_qtl_files_all(caqtl_output_loc)

# get the variants where we do not have the chromosome or position (this happens due to long variants and how LIMIX handles these)
eqtl_variants_all_wmissings <- eqtl_variants_all[is.na(eqtl_variants_all[['snp_chromosome']]) | is.na(eqtl_variants_all[['snp_position']]), ]
# split up by chromosome, position, ref, alt
eqtl_variants_all_wmissings_info <- data.frame(str_split_fixed(eqtl_variants_all_wmissings[['snp_id']], ':', 4))
# make the first two columns numeric
eqtl_variants_all_wmissings_info[[1]] <- as.numeric(eqtl_variants_all_wmissings_info[[1]])
eqtl_variants_all_wmissings_info[[2]] <- as.numeric(eqtl_variants_all_wmissings_info[[2]])
# and update the values in the original table
eqtl_variants_all[is.na(eqtl_variants_all[['snp_chromosome']]) | is.na(eqtl_variants_all[['snp_position']]), c('snp_chromosome', 'snp_position')] <- eqtl_variants_all_wmissings_info[, c(1,2)]

# get the variants where we do not have the chromosome or position (this happens due to long variants and how LIMIX handles these)
caqtl_variants_all_wmissings <- caqtl_variants_all[is.na(caqtl_variants_all[['snp_chromosome']]) | is.na(caqtl_variants_all[['snp_position']]), ]
# split up by chromosome, position, ref, alt
caqtl_variants_all_wmissings_info <- data.frame(str_split_fixed(caqtl_variants_all_wmissings[['snp_id']], ':', 4))
# make the first two columns numeric
caqtl_variants_all_wmissings_info[[1]] <- as.numeric(caqtl_variants_all_wmissings_info[[1]])
caqtl_variants_all_wmissings_info[[2]] <- as.numeric(caqtl_variants_all_wmissings_info[[2]])
# and update the values in the original table
caqtl_variants_all[is.na(caqtl_variants_all[['snp_chromosome']]) | is.na(caqtl_variants_all[['snp_position']]), c('snp_chromosome', 'snp_position')] <- caqtl_variants_all_wmissings_info[, c(1,2)]

# sort now we have this new info
eqtl_variants_all <- eqtl_variants_all[order(eqtl_variants_all[['snp_chromosome']], eqtl_variants_all[['snp_position']], eqtl_variants_all[['snp_id']]), ]
caqtl_variants_all <- caqtl_variants_all[order(caqtl_variants_all[['snp_chromosome']], caqtl_variants_all[['snp_position']], caqtl_variants_all[['snp_id']]), ]

# where to store the result
eqtl_variants_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_eqtl_variants_tested.tsv.gz'
caqtl_variants_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_caqtl_variants_tested.tsv.gz'

# write these to files
write.table(eqtl_variants_all, gzfile(eqtl_variants_loc), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(caqtl_variants_all, gzfile(caqtl_variants_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# and make checksum
mdfiver::create_md5_for_file(eqtl_variants_loc)
mdfiver::create_md5_for_file(caqtl_variants_loc)

# now merge the eQTL and caQTLs
qtl_variants_all <- rbind(eqtl_variants_all, caqtl_variants_all)
# get the unique ones
qtl_variants_all <- unique(qtl_variants_all)
# sort them again
qtl_variants_all <- qtl_variants_all[order(qtl_variants_all[['snp_chromosome']], qtl_variants_all[['snp_position']], qtl_variants_all[['snp_id']]), ]
# and write these as well
qtl_variants_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested.tsv.gz'
write.table(qtl_variants_all, gzfile(qtl_variants_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_md5_for_file(qtl_variants_loc)

# get the numbers
nrow(qtl_variants_all)
# 7091383
nrow(eqtl_variants_all)
# 7068547
nrow(caqtl_variants_all)
# 6867159

# get the cpeaks annotations
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_info.tsv'
# read the cpeaks annotation
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = ' ')
# add signac style identifier
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')

# qtl variants that are multi-allelic, might be split up, causing duplicate positions. we'll make a mapping of the variants to the positions, so we can temporarily deduplicate the positions
qtl_variants_all[['chrompos']] <- paste(qtl_variants_all[['snp_chromosome']], qtl_variants_all[['snp_position']])
# by extracting just based on position, and making that unique
qtl_variant_positions <- unique(qtl_variants_all[, c('chrompos', 'snp_chromosome', 'snp_position')])
# and we'll add 'chr' to the chromosome, as that is what cpeaks does as well
qtl_variant_positions[['snp_chromosome']] <- paste0('chr', qtl_variant_positions[['snp_chromosome']])

# now we'll do the overlap with cpeaks
qtl_variant_pos_to_cpeaks <- overlap_all_regions(
  region_table1 = cpeaks_anno, 
  region_table2 = qtl_variant_positions, 
  feature_id_column1='signac_hg38', 
  feature_id_column2='chrompos', 
  chromosome_column1='chr_hg38', 
  chromosome_column2='snp_chromosome', 
  start_column1='start_hg38', 
  start_column2='snp_position', 
  end_column1='end_hg38', 
  end_column2='snp_position')

# now join that information back onto the origina position table
qtl_variants_all_cpeaks <- merge(qtl_variants_all, data.table(qtl_variant_pos_to_cpeaks[, c('chrompos', 'overlapping_feature')]), by = 'chrompos', all.x = T)

# write result to a file
qtl_variants_all_cpeaks_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_cpeaks_overlap.tsv.gz'
write.table(qtl_variants_all_cpeaks[, c('snp_id', 'snp_chromosome', 'snp_position', 'overlapping_feature')], gzfile(qtl_variants_all_cpeaks_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_md5_for_file(qtl_variants_all_cpeaks_loc)

# get the screen annotations
screen_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/encode_cres/v4/GRCh38-cCREs.bed'
screen_anno <- fread(screen_anno_loc, header = F, sep = '\t')
colnames(screen_anno) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
# add name based on the location
screen_anno[['signac_hg38']] <- paste(screen_anno[['chromosome']], screen_anno[['start']], screen_anno[['end']], sep = '-')

# and overlap with screen this time
qtl_variant_pos_to_screen <- overlap_all_regions(
  region_table1 = screen_anno, 
  region_table2 = qtl_variant_positions, 
  feature_id_column1='signac_hg38', 
  feature_id_column2='chrompos', 
  chromosome_column1='chromosome', 
  chromosome_column2='snp_chromosome', 
  start_column1='start', 
  start_column2='snp_position', 
  end_column1='end', 
  end_column2='snp_position')

# now join that information back onto the origina position table
qtl_variants_all_screen <- merge(qtl_variants_all, data.table(qtl_variant_pos_to_screen[, c('chrompos', 'overlapping_feature')]), by = 'chrompos', all.x = T)
# write result to a file
qtl_variants_all_screen_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/mo_qtl_variants_tested_screen_overlap.tsv.gz'
write.table(qtl_variants_all_screen[, c('snp_id', 'snp_chromosome', 'snp_position', 'overlapping_feature')], gzfile(qtl_variants_all_screen_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_md5_for_file(qtl_variants_all_screen_loc)
