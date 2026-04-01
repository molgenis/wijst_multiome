### cpeaks-based ATAC data processing
'*cpeaks_preprocess_samples/mo_get_credible_atac_barcodes.R*'   get the barcodes from archR that passed QC\
'*cpeaks_preprocess_samples/mo_get_arc_metadata.R*' get the cellranger-arc metadata from the summary csvs\
'*cpeaks_preprocess_samples/mo_create_cpeaks_jobs.sh*'  use the fragments from cellranger, and transform them to cpeaks reference mapped fragments\
'*cpeaks_preprocess_samples/mo_cpeaks_to_signac.R*' read cpeaks reference-based fragments and read them into Signac\
'*cpeaks_preprocess_samples/mo_cpeaks_to_celltypes.R*'  create multimodal object per cell type, and write beds of the peak data\
'*cpeaks_preprocess_samples/mo_cpeaks_to_conditions.R*'  write beds of the peak data, seperated by condition and cell type\
'*cpeaks_preprocess_samples/mo_subset_fragment_files.ipynb*'    subset a fragments file based on the annotations of the barcodes, so check for biases\
'*cpeaks_preprocess_samples/mo_create_multimodal_object_per_celltype.R*'    merge the Seurat countdata and the Signac ATAC data into multimodal objects per cell type\
'*cpeaks_preprocess_samples/mo_multimodal_clustering.R*'    perform PCA, knn-clustering and 2d UMAP using both modalities separately, and together\
'*cpeaks_preprocess_samples/mo_merge_cpeaks_with_screenv4.R*'    merge cpeaks defined regions with screen v4 annotations of regions

