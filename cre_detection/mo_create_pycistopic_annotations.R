#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_pycistopic_annotations.R
# Function: create region annotation file for pycistopic
############################################################################################################################

####################
# libraries        #
####################

library(Signac)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(mdfiver)

####################
# Main Code        #
####################

# location of the annotations
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_info.tsv'
# read the annotations
cpeaks_anno <- read.table(cpeaks_anno_loc, header = T, sep = ' ')

# select the regions we included
features_file_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/merged_major_celltypes/features.tsv.gz'
# read the features
features <- read.table(features_file_loc, header = F)$V1

# set annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# load a dummy object
dc_object <- readRDS('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_dc_wstatus_1_80_20240709.rds')

# get annotations for all of the features
region_to_gene <- ClosestFeature(object = dc_object, regions = StringToGRanges(regions = features), annotation = annotations)
# rename some of the columns
colnames(region_to_gene) <- c('closest_gene_transcript_id', 'closest_gene_name', 'closest_gene_id', 'closest_gene_biotype', 'closest_gene_type', 'gene_closest_region', 'query_region', 'closest_gene_distance')
# reorder them
region_to_gene <- region_to_gene[, c('closest_gene_transcript_id', 'closest_gene_name', 'closest_gene_id', 'closest_gene_biotype', 'closest_gene_type', 'gene_closest_region', 'closest_gene_distance', 'query_region' )]
# reorder the rows to match the features, just to make sure
region_to_gene <- region_to_gene[match(features, region_to_gene[['query_region']]), ]

# remove any NA positions from the cpeaks annotations
cpeaks_anno <- cpeaks_anno[!is.na(cpeaks_anno[['chr_hg38']]) & !is.na(cpeaks_anno[['start_hg38']]) & !is.na(cpeaks_anno[['end_hg38']]), ]
# add that as a name
cpeaks_anno <- cbind(data.frame(name = paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']],cpeaks_anno[['end_hg38']], sep = '-')), cpeaks_anno)
# also order that to be the same as our features
cpeaks_anno <- cpeaks_anno[match(features, cpeaks_anno[['name']]), ]

# remove the query region from the region-to-gene, as we don't need that any more
region_to_gene[['query_region']] <- NULL

# merge the two tables
regions_anno_full <- cbind(cpeaks_anno, region_to_gene)

# put the result somewhere
regions_anno_full_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/merged_major_celltypes/region_annotations.tsv.gz'
write.table(regions_anno_full, gzfile(regions_anno_full_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# checksum as well
mdfiver::create_md5_for_file(regions_anno_full_loc)