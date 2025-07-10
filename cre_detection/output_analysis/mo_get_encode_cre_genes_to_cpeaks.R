#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_encode_cre_genes_to_cpeaks.R
# Function: use the encode ccre-to-gene links to annotate genes to the cpeaks reference
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(ggplot2)
library(IRanges)


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
# Main Code        #
####################

# location of the cpeaks annotations
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/cPeaks/cPeaks_info.tsv'
# read the table
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = ' ')

# read the screen annotation
screen_anno_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/encode_cres/v4/GRCh38-cCREs.bed'
screen_anno <- fread(screen_anno_loc, header = F, sep = '\t')
colnames(screen_anno) <- c('chromosome', 'start', 'end', 'id1', 'id2', 'type')
# add name based on the location
cpeaks_anno[['name']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
screen_anno[['name']] <- paste(screen_anno[['chromosome']], screen_anno[['start']], screen_anno[['end']], sep = '-')

# hi-c gene links
hic_gene_links_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/encode_cres/v4/human_gene_links/V4-hg38.Gene-Links.3D-Chromatin.txt.gz'
hic_gene_links <- fread(hic_gene_links_loc, header = F, sep = '\t')
# set the column names
colnames(hic_gene_links) <- c('cCRE ID','Gene ID','Common Gene Name','Gene Type','Assay Type','Experiment ID','Biosample','Score','P-value')

# CRISPR gene links
crispy_gene_links_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/encode_cres/v4/human_gene_links/V4-hg38.Gene-Links.CRISPR.txt.gz'
crispy_gene_links <- fread(crispy_gene_links_loc, header = F, sep = '\t')
# set the column names
colnames(crispy_gene_links) <- c('cCRE ID','Gene ID','Common Gene Name','Gene Type','gRNA ID','Assay Type','Experiment ID','Biosample','Effect size','P-value')

# keep only K562
hic_gene_links <- hic_gene_links[!is.na(hic_gene_links[['Biosample']]) & hic_gene_links[['Biosample']] == 'K562', ]
crispy_gene_links <- crispy_gene_links[!is.na(crispy_gene_links[['Biosample']]) & crispy_gene_links[['Biosample']] == 'K562', ]
# filter the crispr on p value
crispy_gene_links <- crispy_gene_links[crispy_gene_links[['P-value']] < 0.05, ]
# correlate the scores to the p values
ggplot(data = hic_gene_links[!is.na(hic_gene_links[['P-value']]), ], mapping = aes(x =`Score`, y = -log10(`P-value`))) + 
  geom_point()
# and distribution of scores
ggplot(data = hic_gene_links, mapping = aes(y = Score)) + 
  geom_density(fill = 'blue')

# read the names of the genes
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno <- fread(gene_anno_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')

# add gene symbols to the genes in hic and crispr
hic_gene_links[['gs']] <- gene_anno[match(hic_gene_links[['Gene ID']], gene_anno[['ens']]), ][['gs']]
crispy_gene_links[['gs']] <- gene_anno[match(crispy_gene_links[['Gene ID']], gene_anno[['ens']]), ][['gs']]

# remove any genes without gene symbols
hic_gene_links <- hic_gene_links[!is.na(hic_gene_links[['gs']]), ]
crispy_gene_links <- crispy_gene_links[!is.na(crispy_gene_links[['gs']]), ]

# add the genes to the screen data
screen_with_genes <- merge(screen_anno, hic_gene_links, by.x = 'id2', by.y = 'cCRE ID', all.x = F, all.y = F)
# subset to just the name and gene symbol, to aggregate genes for the same region
screen_with_genes <- screen_with_genes[, c('name', 'gs')][, .(gs = paste(unique(gs), collapse = ";")), by = name]
# finally add this back to the original screen data
screen_anno[['hic_gs']] <- screen_with_genes[match(screen_anno[['name']], screen_with_genes[['name']])][['gs']]
# remove what we don't have
screen_anno <- screen_anno[!is.na(screen_anno[['hic_gs']]), ]
# add this hi-c annotation to cpeaks
cpeaks_anno <- add_annotations_from_list(cpeaks_anno, list(
  # add the hi-C genes as a column under the name 'type', because that is what the method expects
  'hic_genes' = cbind(screen_anno[, -c('type'), drop = F], data.frame('type' = screen_anno[['hic_gs']]))
))

# do the same for the crispr links
screen_with_genes <- merge(screen_anno, crispy_gene_links, by.x = 'id2', by.y = 'cCRE ID', all.x = F, all.y = F)
# subset to just the name and gene symbol, to aggregate genes for the same region
screen_with_genes <- screen_with_genes[, c('name', 'gs')][, .(gs = paste(unique(gs), collapse = ";")), by = name]
# finally add this back to the original screen data
screen_anno[['crispr_gs']] <- screen_with_genes[match(screen_anno[['name']], screen_with_genes[['name']])][['gs']]
# remove what we don't have
screen_anno <- screen_anno[!is.na(screen_anno[['crispr_gs']]), ]
# add this crispr annotation to cpeaks
cpeaks_anno <- add_annotations_from_list(cpeaks_anno, list(
  # add the hi-C genes as a column under the name 'type', because that is what the method expects
  'crispr_genes' = cbind(screen_anno[, -c('type'), drop = F], data.frame('type' = screen_anno[['crispr_gs']]))
))

# split the region-gene pairs up as well
screen_cpeaks_crispr_links <- cpeaks_anno[!is.na(cpeaks_anno[['crispr_genes']]), c('name', 'crispr_genes')][, .(gene = unlist(strsplit(crispr_genes, ";"))), by = name]
screen_cpeaks_hic_links <- cpeaks_anno[!is.na(cpeaks_anno[['hic_genes']]), c('name', 'hic_genes')][, .(gene = unlist(strsplit(hic_genes, ";"))), by = name]
# rename the columns
colnames(screen_cpeaks_crispr_links) <- c('region', 'gene')
colnames(screen_cpeaks_hic_links) <- c('region', 'gene')
# and we'll put these somewhere
screen_cpeaks_crispr_links_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/cPeaks/cpeaks_to_screenv4_crispr.tsv.gz'
screen_cpeaks_crispr_hic_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/cPeaks/cpeaks_to_screenv4_hic.tsv.gz'
# save the tables
write.table(screen_cpeaks_crispr_links, gzfile(screen_cpeaks_crispr_links_loc), row.names = F, col.names = T, sep = '\t', quote = F)
write.table(screen_cpeaks_hic_links, gzfile(screen_cpeaks_crispr_hic_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# make some checksums as well
mdfiver::create_md5_for_file(screen_cpeaks_crispr_links_loc)
mdfiver::create_md5_for_file(screen_cpeaks_crispr_hic_loc)

# we'll also save the cpeaks annotations to the genes
cpeaks_anno_wscreen_genes_loc <- '/groups/umcg-franke-scrna/tmp02/external_datasets/cPeaks/cPeaks_wk562_genes_screenv4.tsv.gz'
write.table(cpeaks_anno, gzfile(cpeaks_anno_wscreen_genes_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# and a checksum to top it off
mdfiver::create_md5_for_file(cpeaks_anno_wscreen_genes_loc)
