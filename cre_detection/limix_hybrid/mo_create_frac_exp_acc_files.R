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
