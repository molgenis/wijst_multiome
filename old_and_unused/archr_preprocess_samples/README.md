### ArchR ATAC data processing
'*archr_preprocess_samples/build-arrowfiles.R*' create Arrow files to use in archR\
'*archr_preprocess_samples/build-project.R*'  create ArchR project from Arrow chunk_files\
'*archr_preprocess_samples/batch_correction.R*'  do batch correction over all lanes (unfortunately names 'sample' in ArchR)\
'*archr_preprocess_samples/preprocess.R*'  perform preprocesing to remove low-quality nuclei and doublets\
'*archr_preprocess_samples/iterative-LSI.R*'  perform dimension reduction and clustering\
'*archr_preprocess_samples/iterative-LSI.R*'  perform dimension reduction and clustering\
'*archr_preprocess_samples/addMetadata_subset_celltypes.R*' add celltype annotation from Seurat to ArchR project\
'*archr_preprocess_samples/imputeCelltypes.R*' impute missing celltypes for nuclei, based on Seurat annotation clusters (majority vote)\
'*archr_preprocess_samples/peakCalling.R*' perform peak calling based on celltypes assigned
