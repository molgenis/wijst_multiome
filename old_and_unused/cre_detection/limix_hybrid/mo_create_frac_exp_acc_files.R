#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_frac_exp_acc_files.R
# Function: create files that list the fraction of cells expressing a gene, or having a read in a chromatin region
############################################################################################################################

####################
# libraries        #
####################

# for reading the data
library(Seurat)
library(Signac)

####################
# Functions        #
####################


add_binarized_assay <- function(seurat_object, assay_to_binarize='peaks', layer_to_binarize='counts', assay_to_add='binpeaks') {
  # set the default assay
  DefaultAssay(seurat_object) <- assay_to_binarize
  # extract the assay
  assay_data <- Seurat::GetAssayData(seurat_object, layer = layer_to_binarize, assay = assay_to_binarize)
  # binarize the chromatin data if requested
  assay_data@x[assay_data@x > 0] <- 1
  # add the assay data back
  if (layer_to_binarize == 'counts') {
    seurat_object[[assay_to_add]] <- CreateAssayObject(counts = assay_data, cells = colnames(assay_data))
  }
  else {
    seurat_object[[assay_to_add]] <- CreateAssayObject(data = assay_data, cells = colnames(assay_data))
  }
  return(seurat_object)
}



####################
# settings         #
####################

set.seed(7777)


####################
# debug code      #
####################


####################
# Main Code        #
####################

# location of the Seurat objects
objects_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
# prepend and append of objects
object_prepend <- 'mo_multimodal_'
object_append <- '_1_80_20240521.rds'

# set the location of the output file for expression
fracs_exp_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/frac_exp/'
# the prepend of this file
fracs_exp_prepend <- ''
# the append of this file
fracs_exp_append <- '.tsv.gz'

# same for the fraction of accessible
fracs_acc_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/frac_acc/'
# the prepend of this file
fracs_acc_prepend <- ''
# the append of this file
fracs_acc_append <- '.tsv.gz'

# these are the cell types to use
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')

# check each cell type
for (cell_type in cell_types) {
  # past the path together
  object_full_loc <- paste0(objects_loc, '/', object_prepend, tolower(cell_type), object_append)
  # read the object
  object <- readRDS(object_full_loc)
  # add the binarized chromatin assay
  object <- add_binarized_assay(object)
  # add the binarized expression assay
  object <- add_binarized_assay(object, assay_to_binarize = 'RNA', layer_to_binarize = 'counts', assay_to_add = 'binexp')
  
  
  # we can now calculate the average expression, which will be the fraction expressed, because we made every value >1 into 1
  frac_exp_nonzero <- Matrix::rowMeans(object@assays$binexp@counts)
  # make into the format we expect
  frac_exp_nonzero <- data.frame('feature' = names(frac_exp_nonzero), 'frac_exp' = as.vector(frac_exp_nonzero))
  # order by feature
  frac_exp_nonzero <- frac_exp_nonzero[order(frac_exp_nonzero[['feature']]), ]
  # get the location of where to store
  fracs_exp_loc_full <- paste0(fracs_exp_loc, '/', fracs_exp_prepend, cell_type, fracs_exp_append)
  # write the file
  write.table(frac_exp_nonzero, gzfile(fracs_exp_loc_full), row.names = F, col.names = T, sep = '\t', quote = F)
  # and make a checksum
  mdfiver::create_md5_for_file(fracs_exp_loc_full)
  
  # we can now calculate the average accessibility, which will be the fraction expressed, because we made every value >1 into 1
  frac_acc_nonzero <- Matrix::rowMeans(object@assays$binpeaks@counts)
  # make into the format we accect
  frac_acc_nonzero <- data.frame('feature' = names(frac_acc_nonzero), 'frac_acc' = as.vector(frac_acc_nonzero))
  # order by feature
  frac_acc_nonzero <- frac_acc_nonzero[order(frac_acc_nonzero[['feature']]), ]
  # get the location of where to store
  fracs_acc_loc_full <- paste0(fracs_acc_loc, '/', fracs_acc_prepend, cell_type, fracs_acc_append)
  # write the file
  write.table(frac_acc_nonzero, gzfile(fracs_acc_loc_full), row.names = F, col.names = T, sep = '\t', quote = F)
  # and make a checksum
  mdfiver::create_md5_for_file(fracs_acc_loc_full)
  
  # clear memory
  rm(object)
}

# check each cell type
for (cell_type in cell_types) {
  # past the path together
  object_full_loc <- paste0(objects_loc, '/', object_prepend, tolower(cell_type), object_append)
  # read the object
  object <- readRDS(object_full_loc)
  # add the binarized chromatin assay
  object <- add_binarized_assay(object)
  # add the binarized expression assay
  object <- add_binarized_assay(object, assay_to_binarize = 'RNA', layer_to_binarize = 'counts', assay_to_add = 'binexp')
  # extract the expression matrix
  object_exp <- object@assays$binexp@counts
  # extract the accessiblity matrix
  object_acc <- object@assays$binpeaks@counts
  # add sample and lane
  object@meta.data[['sample_lane']] <- paste(object@meta.data[['sample_final']], object@meta.data[['lane']], sep = ';;')
  # extract the samples 
  samples <- unique(object@meta.data[['sample_lane']])
  # extract the features
  exp_features <- rownames(object_exp)
  # extract the features
  acc_features <- rownames(object_acc)
  # save binary matrix per sample
  frac_exp_nonzero_psample <- list()
  frac_acc_nonzero_psample <- list()
  # add those features as a column
  frac_exp_nonzero_psample[['feature']] <- data.frame('feature' = exp_features)
  frac_acc_nonzero_psample[['feature']] <- data.frame('feature' = acc_features)
  # check each sample
  for (sample_mtdt in samples) {
    # get indices where it is this sample
    sample_this_sample <- object@meta.data[['sample_lane']] == sample_mtdt
    
    # subset the matrix
    object_exp_sample <- object_exp[, sample_this_sample]
    # now the fraction of nonzero
    frac_exp_nonzero <- NULL
    if (!is.null(dim(object_exp_sample))) {
      frac_exp_nonzero <- Matrix::rowMeans(object_exp_sample)
    }
    else {
      frac_exp_nonzero <- object_exp_sample
    }
    # make into df
    frac_exp_nonzero <- data.frame(x = as.vector(unlist(frac_exp_nonzero)))
    # set the column to be the sample
    colnames(frac_exp_nonzero) <- c(sample_mtdt)
    # and add to the list
    frac_exp_nonzero_psample[[sample_mtdt]] <- frac_exp_nonzero
    
    # subset the matrix
    object_acc_sample <- object_acc[, sample_this_sample]
    frac_acc_nonzero <- NULL
    if (!is.null(dim(object_acc_sample))) {
      frac_acc_nonzero <- Matrix::rowMeans(object_acc_sample)
    }
    else {
      frac_acc_nonzero <- object_acc_sample
    }
    # make into df
    frac_acc_nonzero <- data.frame(x = as.vector(unlist(frac_acc_nonzero)))
    # set the column to be the sample
    colnames(frac_acc_nonzero) <- c(sample_mtdt)
    # and add to the list
    frac_acc_nonzero_psample[[sample_mtdt]] <- frac_acc_nonzero
  }
  
  # merge all the dfs in the list
  frac_exp_nonzero_asample <- do.call('cbind', frac_exp_nonzero_psample)
  frac_acc_nonzero_asample <- do.call('cbind', frac_acc_nonzero_psample)
  
  # get the location of where to store
  fracs_exp_asample_loc_full <- paste0(fracs_exp_loc, '/', fracs_exp_prepend, cell_type, '_persample', fracs_exp_append)
  # write the file
  write.table(frac_exp_nonzero_asample, gzfile(fracs_exp_asample_loc_full), row.names = F, col.names = T, sep = '\t', quote = F)
  # and make a checksum
  mdfiver::create_md5_for_file(fracs_exp_asample_loc_full)
  
  # get the location of where to store
  fracs_acc_asample_loc_full <- paste0(fracs_acc_loc, '/', fracs_acc_prepend, cell_type, '_persample', fracs_acc_append)
  # write the file
  write.table(frac_acc_nonzero_asample, gzfile(fracs_acc_asample_loc_full), row.names = F, col.names = T, sep = '\t', quote = F)
  # and make a checksum
  mdfiver::create_md5_for_file(fracs_acc_asample_loc_full)
}
