#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_hybrid_cre_inputs.R
# Function: 
#
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(Signac)
library(data.table)


####################
# Functions        #
####################


#' Write Data Sample
#'
#' This function writes expression data, chromatin data, and metadata to specified folders, ensuring that the sample name is POSIX-safe. It also creates MD5 checksums for the written files.
#'
#' @param expression_data A matrix containing the expression data.
#' @param chromatin_data A matrix containing the chromatin data.
#' @param metadata A data frame containing the metadata.
#' @param sample_name A character string representing the sample name.
#' @param output_folder A character string specifying the output folder path.
#' @param expression_name A character string for the expression data folder name. Default is 'RNA'.
#' @param chromatin_name A character string for the chromatin data folder name. Default is 'peaks'.
#' @param binarize_chromatin_matrix A logical value indicating whether to binarize the chromatin matrix. Default is TRUE.
#' @return An integer value of 0 upon successful completion.
#' @examples
#' \dontrun{
#' write_data_sample(expression_data, chromatin_data, metadata, "sample1", "/path/to/output")
#' }
#' @export
write_data_sample <- function(expression_data, chromatin_data, metadata, sample_name, output_folder, expression_name='RNA', chromatin_name='peaks', binarize_chromatin_matrix=T) {
  # make the sample name posix safe
  sample_name_posix <- make_posix_safe(sample_name)
  # warn if this make the name different
  if (sample_name != sample_name_posix) {
    warning(paste('sample name was not POSIX safe,', sample_name, 'was renamed to', sample_name_posix))
  }
  # create folder for the expression data
  expression_data_folder <- paste0(output_folder, '/', sample_name_posix, '/', expression_name, '/')
  # the chromatin folder as well
  chromatin_data_folder <- paste0(output_folder, '/', sample_name_posix, '/', chromatin_name, '/')
  # create the directories
  dir.create(expression_data_folder, recursive = T)
  dir.create(chromatin_data_folder, recursive = T)
  
  # add barcode as explicit column
  metadata <- cbind(data.frame('barcode' = rownames(metadata)), metadata)
  # write the metadata
  metadata_loc <- paste0(output_folder, '/', sample_name_posix, '/metadata.tsv.gz')
  # write it
  write.table(metadata, gzfile(metadata_loc), row.names = F, col.names = T, quote = F, sep = '\t')
  # and make an md5
  mdfiver::create_md5_for_file(metadata_loc)
  
  # write the expression data
  expression_barcodes_loc <- paste0(expression_data_folder, 'barcodes.tsv.gz')
  write.table(data.frame(x = colnames(expression_data)), gzfile(expression_barcodes_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(expression_barcodes_loc)
  expression_features_loc <- paste0(expression_data_folder, 'features.tsv.gz')
  write.table(data.frame(x = rownames(expression_data)), gzfile(expression_features_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(expression_features_loc)
  # write matrix
  expression_matrix_loc <- paste0(expression_data_folder, 'matrix.mtx')
  writeMM(expression_data, expression_matrix_loc)
  # Compress the MTX file using gzip-0  
  expression_matrix_gz_loc <- paste0(expression_data_folder, 'matrix.mtx.gz')
  gzip(expression_matrix_loc, expression_matrix_gz_loc, overwrite = TRUE)
  # and make md5 checksum
  mdfiver::create_md5_for_file(expression_matrix_gz_loc)
  
  # binarize the chromatin data if requested
  if (binarize_chromatin_matrix) {
    chromatin_data@x[chromatin_data@x > 1] <- 1
  }
  
  # and write the chromatin data
  chromatin_barcodes_loc <- paste0(chromatin_data_folder, 'barcodes.tsv.gz')
  write.table(data.frame(x = colnames(chromatin_data)), gzfile(chromatin_barcodes_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(chromatin_barcodes_loc)
  chromatin_features_loc <- paste0(chromatin_data_folder, 'features.tsv.gz')
  write.table(data.frame(x = rownames(chromatin_data)), gzfile(chromatin_features_loc), row.names = F, col.names = F, quote = F)
  mdfiver::create_md5_for_file(chromatin_features_loc)
  # write matrix
  chromatin_matrix_loc <- paste0(chromatin_data_folder, 'matrix.mtx')
  writeMM(chromatin_data, chromatin_matrix_loc)
  # Compress the MTX file using gzip
  chromatin_matrix_gz_loc <- paste0(chromatin_data_folder, 'matrix.mtx.gz')
  gzip(chromatin_matrix_loc, chromatin_matrix_gz_loc, overwrite = TRUE)
  # and make md5 checksum
  mdfiver::create_md5_for_file(chromatin_matrix_gz_loc)
  return(0)
}


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
    seurat_object[['data']] <- CreateAssayObject(data = norm_count_matrix)
    
  } else {
    print('using Seurat v3/4 style \'slot\'')
    seurat_object[['data']] <- CreateAssayObject(data = norm_count_matrix)
  }
  return(seurat_object)
}


write_matrix_slices <- function(multimodal_object, output_loc, chromosome, gene_anno, region_anno, confinement, gene_anno_gene_column='gs', gene_anno_chrom_column='chrom', gene_anno_start_column='start', gene_anno_end_column='end', confinement_gene_column='feature_id', confinement_region_column='snp_id', gene_chunk_size=100, assay_expression='data', layer_expression='data', assay_accessibility='peaks', layer_accessibility='counts', binarize_chromatin_matrix=T) {
  # get expression assay
  expression <- GetAssayData(multimodal_object, layer = layer_expression, assay = assay_expression)
  accessibility <- GetAssayData(multimodal_object, layer = layer_accessibility, assay = assay_accessibility)
  # binarize the chromatin data if requested
  if (binarize_chromatin_matrix) {
    accessibility@x[accessibility@x > 1] <- 1
  }
  # subset the gene annotations to that chromosome
  gene_annos_chromosome <- gene_anno[!is.na(gene_anno[[gene_anno_chrom_column]]) &
                                       gene_anno[[gene_anno_chrom_column]] == chromosome, ]
  # and make sure it is ordered as we expect
  gene_annos_chromosome <- gene_annos_chromosome[
    order(gene_annos_chromosome[[gene_anno_start_column]], gene_annos_chromosome[[gene_anno_end_column]])
  ]
  # get all unique genes
  chromosome_genes <- unique(gene_annos_chromosome[[gene_anno_gene_column]])
  # go through the chunks
  chunk_start <- 1
  # message process
  message(paste('processing', as.character(length(chromosome_genes)), 'features in chunks of', as.character(gene_chunk_size)))
  # keep checking each chunk
  while(chunk_start < length(chromosome_genes)) {
    # the end of the chunk
    chunk_end <- chunk_start + gene_chunk_size - 1
    # unless we don't have a full chunk left
    if (chunk_end > length(chromosome_genes)) {
      chunk_end <- length(chromosome_genes)
    }
    # message
    message(paste('processing chunk', chunk_start, 'to', chunk_end))
    # grab the genes of this chunk
    chunk_genes <- chromosome_genes[chunk_start : chunk_end]

    # now get the chromatin regions for those genes
    chunk_regions <- unique(confinement[!is.na(confinement[[confinement_gene_column]]) &
                                   !is.na(confinement[[confinement_region_column]]) &
                                   confinement[[confinement_gene_column]] %in% chunk_genes , ][[confinement_region_column]])
    
    # subset the annotation to these genes
    gene_annos_chunk <- gene_annos_chromosome[gene_annos_chromosome[[gene_anno_gene_column]] %in% chunk_genes, ]
    
    
    # get the smallest number from the starts and ends
    chunk_min <- min(c(gene_annos_chunk[[gene_anno_start_column]], gene_annos_chunk[[gene_anno_end_column]]))
    # and the largest
    chunk_max <- max(c(gene_annos_chunk[[gene_anno_start_column]], gene_annos_chunk[[gene_anno_end_column]]))
    
    # subset to these genes
    assay_genes <- expression[chunk_genes[chunk_genes %in% rownames(expression)], ]
    # and regions
    assay_regions <- accessibility[chunk_regions[chunk_regions %in% rownames(accessibility)], ]
    
    # convert the gene assay to a non-sparse format
    genes_ns <- as.matrix(assay_genes)
    # and the accessibility
    regions_ns <- as.matrix(assay_regions)
    
    # create full path of the matrix
    output_directory_chunk <- paste0(output_loc, '/', chromosome, '-', chunk_min, '-', chunk_max, '/')
    # create that directory
    dir.create(output_directory_chunk, recursive = T)
    
    # create the output file location
    expression_chunk_output_loc <- paste0(output_directory_chunk, '/', 'expression.tsv.gz')
    accessibility_chunk_output_loc <- paste0(output_directory_chunk, '/', 'accessibility.tsv.gz')
    
    # actually write the files now
    write.table(genes_ns, gzfile(expression_chunk_output_loc), row.names = T, col.name = T, sep = '\t', quote = F)
    write.table(regions_ns, gzfile(accessibility_chunk_output_loc), row.names = T, col.name = T, sep = '\t', quote = F)
    
    # and finally make some checksums
    mdfiver::create_md5_for_file(expression_chunk_output_loc)
    mdfiver::create_md5_for_file(accessibility_chunk_output_loc)
    
    # update chunk
    chunk_start <- chunk_start + gene_chunk_size
  }
  return(0)
}



####################
# Main code        #
####################

# location of the region-to-gene files
region_to_gene_confinement_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/featureVariantFile.w150k.filtered0.0001_cts.txt'
# read that file
region_to_gene_confinement <- fread(region_to_gene_confinement_loc, header = T, sep = '\t')

# read the cpeaks annotation
cpeaks_anno_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/cPeaks/cPeaks_info.tsv'
cpeaks_anno <- fread(cpeaks_anno_loc, header = T, sep = ' ')
# add the Signac style name
cpeaks_anno[['signac_hg38']] <- paste(cpeaks_anno[['chr_hg38']], cpeaks_anno[['start_hg38']], cpeaks_anno[['end_hg38']], sep = '-')
# as well as the SCENIC+ style name
cpeaks_anno[['scenic_hg38']] <- paste0(cpeaks_anno[['chr_hg38']], ':', cpeaks_anno[['start_hg38']], '-', cpeaks_anno[['end_hg38']])

# read the location of the genes
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno <- fread(gene_anno_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')

# location of the cell type object
mono_object_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_monocyte_1_80_20240521.rds'

# location of the output
output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/L1/monocyte/'

# donor annotation column
donor_annotation_column <- 'sample_final'

# test monocyte for now
mono_object <- readRDS(mono_object_loc)

# add pflogpf
mono_object <- normalize_mj(mono_object)

# write the slices
write_matrix_slices(mono_object, output_loc, 'chr17', gene_anno, cpeaks_anno, region_to_gene_confinement)

# get a cell barcode to sample
smf <- data.frame(genotype_id = mono_object@meta.data[[donor_annotation_column]], phenotype_id = rownames(mono_object@meta.data))
# location of the smf
smf_output_loc <- paste0(output_loc, '/smf.tsv.gz')
# write the result
write.table(smf, gzfile(smf_output_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# make checksum
mdfiver::create_md5_for_file(smf_output_loc)
