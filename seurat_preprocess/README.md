### seurat preprocess
'*seurat_preprocess/mo_lane_to_seurat.R*'   read the barcodes/features/matrix files into Seurat objects for the RNA counts\
'*seurat_preprocess/mo_lane_to_seuratv4.R*'   read the barcodes/features/matrix files into Seurat objects for the RNA counts, using Seurat v4 format

### merge seurat objects
'*seurat_preprocess/mo_merge_seurat_objects.R*'     merge the celltyped Seurat files for each lane\
'*seurat_preprocess/mo_filter_qc.R*'    do QC on the count data\
'*seurat_preprocess/mo_normalize_and_cluster.R*'    perform SCT normalization, PCA, knn-clustering and 2d UMAP\
'*seurat_preprocess/mo_add_rna_metadata.R*'     add covid assignments to seurat metadata\
'*seurat_preprocess/mo_create_l2_objects.R*'    create a Seurat object per L2 Azimuth covid19 annotation

### LONG-COVID
'*seurat_preprocess/lc_create_per_celltype_objects.R*'  create a Seurat object for each cell type, with only UT or both UT and 24hCA\
'*seurat_preprocess/lc_export_rna_objects.R*'  export Seurat objects for the long-covid lifelines data
