#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_signac_peaks_celltype_lane.R
# Function: 
############################################################################################################################


####################
# libraries        #
####################

# it's a seurat object
library(Seurat)
# with ATAC data as well
library(Signac)
# this database is needed for annotations
library(EnsDb.Hsapiens.v86)
# this is for parsing the parameters
library(optparse) # NOT IN CONTAINER

####################
# Functions        #
####################

add_lower_classification <- function(seurat_object, reclassification_mapping, mapping_original_column, mapping_reclass_column, metadata_original_column, metadata_reclassification_column){
  # add the new column
  seurat_object@meta.data[[metadata_reclassification_column]] <- NA
  # get each cell type in the data
  metadata_original_cts <- unique(seurat_object@meta.data[[metadata_original_column]])
  # and the originals in the mapping
  reclassification_original_cts <- unique(reclassification_mapping[[mapping_original_column]])
  # we can only map what is present in both
  originals_both <- intersect(metadata_original_cts, reclassification_original_cts)
  # check what is missing
  only_metadata <- setdiff(metadata_original_cts, reclassification_original_cts)
  only_mapping <- setdiff(reclassification_original_cts, metadata_original_cts)
  # warn what is missing
  if(length(only_metadata) > 0){
    print('some celltypes only in metadata')
    print(only_metadata)
  }
  if(length(only_mapping) > 0){
    print('some celltypes only in remapping ')
    print(only_mapping)
  }
  # check each cell type
  for(celltype_original in originals_both){
    # get the appropriate remapping
    celltype_remapped <- reclassification_mapping[reclassification_mapping[[mapping_original_column]] == celltype_original, mapping_reclass_column]
    # now remap in the metadata
    seurat_object@meta.data[seurat_object@meta.data[[metadata_original_column]] == celltype_original, metadata_reclassification_column] <- celltype_remapped
  }
  return(seurat_object)
}


do_processing <- function(options, annotations) {
  # read lane object
  lane <- readRDS(opt[['seurat']])
  # set the assay
  DefaultAssay(lane) <- 'peaks'
  # set the annotations to the object now
  Annotation(lane) <- annotations
  # hack the path to be vaxtron compatible
  #lane$peaks@fragments[[1]]@path <- '/groups/umcg-franke-scrna/tmp03/projects/multiome/processed/joint/alignment/b38/230302_lane6/outs/atac_fragments.tsv.gz'
  lane$peaks@fragments[[1]]@path <- gsub('tmp02', 'tmp03', lane$peaks@fragments[[1]]@path)
  # compute nucleosome signal score per cell
  lane <- NucleosomeSignal(object = lane)
  # compute TSS enrichment score per cell
  lane <- TSSEnrichment(object = lane, fast = FALSE)
  # show TSS
  lane$high.tss <- ifelse(lane$TSS.enrichment > 3, 'High', 'Low')
  
  # initially set the celltype column
  celltype_column <- opt[['cell_type_column']]
  # add new celltype column if a path to one was supplied
  if (!is.null(opt[['reclassification_table']])) {
    # read that table
    reclassification <- read.table(opt[['reclassification_table']], header = T, sep = '\t')
    # get the column names of the original and the reclassification
    original_column <- colnames(reclassification)[1]
    reclassification_column <- colnames(reclassification)[2]
    # do the reclassification
    lane <- add_lower_classification(lane, reclassification, original_column, reclassification_column, celltype_column, reclassification_column)
    # set the new column as the cell type now
    celltype_column <- reclassification_column
  }
  
  # calculate gene activity
  gene_activities <- GeneActivity(lane)
  
  # subset the object to interested locations if we have a confinement
  if (!is.null(opt[['location_confinement']])) {
    
  }
  
  # set the cell type as the identity
  Idents(lane) <- celltype_column
  # calculate the peaks
  da_peaks <- FindMarkers(
    object = lane,
    test.use = 'LR',
    latent.vars = 'nCount_peaks',
    ident.1 = opt[['cell_type']]
  )
  # now get the closest genes
  da_closest_open_features <- ClosestFeature(lane, regions = rownames(da_peaks))
  # add the open features
  da_peaks <- cbind(da_peaks, 
                    da_closest_open_features[match(rownames(da_peaks), da_closest_open_features[['query_region']]), c('tx_id','gene_name','gene_id','gene_biotype','type','closest_region', 'distance')])
  
  # add the gene and the celltype as explicit columns
  da_peaks <- cbind(data.frame('cell_type' = rep(opt[['cell_type']], times = nrow(da_peaks)), 'region' = rownames(da_peaks)),
                    da_peaks)
  # write all results
  peaks_output_loc <- paste(opt[['out']], '.peaks.tsv', sep = '')
  write.table(da_peaks, peaks_output_loc, sep = '\t', row.names = F, col.names = T, quote = T)
}


do_debug <- function(annotations) {
  options <- list(
    'seurat' = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_230302_lane6multimodal_azi_mapped.rds',
    'out' = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/DC',
    'cell_type_column' = 'predicted.mo_10x_cell_type',
    'cell_type' = 'DC',
    'reclassification_table' = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cell_type_assignment/azimuth/celltypes_10x_ref_to_lowerres.tsv'
  )
  do_processing(options, annotations)
}

####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-s", "--seurat"), type="character", default=NULL,
              help="seurat object input", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file base path without extention [default= %default]", metavar="character"),
  make_option(c("-c", "--cell_type_column"), type="character", default='predicted.mo_10x_cell_type',
              help="column describing the celltype in the Seurat metadata [default= %default]", metavar="character"),
  make_option(c("-t", "--cell_type"), type="character", default=NULL,
              help="the cell type to do the analysis for [default= %default]", metavar="character"),
  make_option(c("-r", "--reclassification_table"), type="character", default=NULL,
              help="celltype reclassification table [default= %default]", metavar="character"),
  make_option(c("-l", "--location_confinement"), type="character", default=NULL,
              help="table with genomic locations to confine the search to [default= %default]", metavar="character")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# get annotations from ensemble database
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# we use UCSC gencode
seqlevelsStyle(annotations) <- "UCSC"
# and this was aligned on build 38
genome(annotations) <- "hg38"

# do debug
#do_debug(annotations)
# do the analysis
do_processing(opt, annotations)
