# Wijst Multi-ome study

This is the github repository of the Multiome study containing CA stimulated PBMCs

## Software

Here we will list the software used to generate the data

Cellranger 7.1.0\
Cellranger-arc 2.0.0\
R 4.3.1\
plink2-20230707\
CellBender v3\
Seurat 4.9.9.9059\
Souporcell v2\
Demuxlet v2

## Custom code

Below we will outline the steps taken to process and analyse the data

### alignment to reference genome

'*alignment/mo_create_cellranger_jobs_batch1.sh*' create SLURM jobs to align the RNA-seq data to the b38 human reference\
'*alignment/mo_create_multiome_csvs_batch1.sh*' create SLURM jobs to do a joint align of the RNA-seq and ATAC-seq data to the b38 human reference\
'*alignment/mo_create_cellranger_arc_jobs_batch1.sh*' create SLURM jobs to align the RNA+ATAC data to the b38 human reference\
'*alignment/mo_download_websummaries.sh*' download the CellRanger web summary files\
'*alignment/mo_get_cellranger_summaries.R*' created an Excel sheet describing the CellRanger results, from the CellRanger output directories

### genotypeing

'*genotyping/mo_collect_previously_genotyped_individuals.R*' collect data on which individuals were previously genotyped\
'*genotyping/mo_get_previous_genotypes.R*' create list of previously genotyped individuals, create sample sheets, create age/sex metadata file\
'*genotyping/mo_preprocess_previous_genotypes.sh*' subset previously genotyped individuals genotype data, and convert to format for imputation\
'*genotyping/mo_create_ugli_psam.R*' create psam file required for imputation, for the previously genotyped individuals\
'*genotyping/mo_create_batch2_psam.R*' create psam file required for imputation, for the second batch of individuals\
'*genotyping/PreImputation_ugli.yaml*' the imputation config used for imputing the previously genotyped UGLI individuals\
'*genotyping/PreImputation_mo_batch2.yaml*' the imputation config used for imputing the second batch of new individuals

### ambient RNA correction
'*ambient_rna_correction/mo_cellbend_rna.sh*' create sbatch jobs to do ambient RNA correction using CellBender on the RNA data

### demultiplexing and doublet detection
'*demultiplexing/mo_run_scrublet.py*' python script to run Scrublet on a CellBender corrected 10x lane\
'*demultiplexing/mo_test_scrublet.ipynb*' jupyter notebook to test Scrublet on individual 10x lanes\
'*demultiplexing/mo_subset_for_genotype_correlations.sh*' script to subset genotypes for each 10x experiment, using the individuals in present in that experiment\
'*demultiplexing/mo_subset_gex_bams_snp_regions.sh*' script to create jobs to subset the gene expression alignment files, to only have reads overlapping variants we genotyped\
'*demultiplexing/mo_create_popscle_sort_vcfs_jobs.sh*' script to create jobs to sort the per-lane VCF files by the order of chromosomes in the alignment files\
'*demultiplexing/mo_create_demuxlet_jobs.sh*' script to create jobs to run Demuxlet on the gene expression alignments\
'*demultiplexing/mo_create_souping_samples_jobs_externalbarcodes.sh*' script to create jobs that perform Souporcell on each 10x lane while allowing external barcode files\
'*demultiplexing/mo_correlate_genotypes.R*' correlate the Souporcell cluster genotypes to the genotypes generated for the individuals, to do sample assignment\
'*demultiplexing/mo_plot_demultiplexing_assignments.Rmd*' plot the demultiplexing assignments\
'*demultiplexing/mo_sample_assignment.R*' do sample assignment on the Seurat object, based on the correlated souporcell genotypes\
'*demultiplexing/mo_get_rna_qced_barcodes.R*' extract the valid barcodes in the Seurat object\

### scanpy preprocess
'*scanpy_preprocess/mo_h5_to_scanpy_objects.ipynb*' read the Seurat-incompatible h5 CellBender outputs into Scanpy H5AD objects\
'*scanpy_preprocess/mo_deconstruct_scanpy_objects.ipynb*' deconstruct the H5AD Scanpy objects into Seurat compatible barcodes/features/matrix files

### seurat preprocess
'*seurat_preprocess/mo_lane_to_seurat.R*' read the barcodes/features/matrix files into Seurat objects for the RNA counts

### cell type assignment
'*cell_type_annotation/mo_create_10x_mo_reference.R*' create a reference Azimuth file from the 10x PBMC multi-ome reference dataset, which will also add the ATAC data using signac\
'*cell_type_annotation/mo_azimuth_reference_map_mo.R*' do an Azimuth reference mapping of each 10x lane onto the 10x PMBC reference\
'*cell_type_annotation/mo_azimuth_reference_map_pbmc.R*' do an Azimuth reference mapping of each 10x lane onto the Azimuth COVID-19 PBMC dataset\
'*cell_type_annotation/mo_check_celltype_assignments.Rmd*' check the cell type assignments that were done

### merge seurat objects
'*seurat_preprocess/mo_merge_seurat_objects.R*' merge the celltyped Seurat files\
'*seurat_preprocess/mo_filter_qc.R*' do QC on the count data\
'*seurat_preprocess/mo_normalize_and_cluster.R*' perform SCT normalization, PCA, knn-clustering and 2d UMAP\
'*seurat_preprocess/mo_add_rna_metadata.R*' add covid assignments to seurat metadata\
'*seurat_preprocess/lc_create_per_celltype_objects.R*' create a Seurat object for each cell type, with only UT or both UT and 24hCA]
'*seurat_preprocess/mo_create_l2_objects.R*' create a Seurat object per L2 Azimuth covid19 annotation

### differential celltype abundance analysis
'*differential_abundance/mo_plot_celltype_proportions.Rmd*' plot the cell type proportions present in the dataset\
'*differential_proportion/mo_differential_proportion_analysis.R*' use Speckle to do differential abundance analysis

### differential gene expression analysis
'*differential_expression/1m_create_per_celltype_objects.R*' create per-celltype Seurat objects from the NC2022 study which also used Candida stimulation\
'*differential_expression/1m_differential_expression_limma_parameterised.R*' perform differential gene expression using limma dream on the stimulation status in the NC2022 data\
'*differential_expression/lc_create_per_celltype_objects.R*' create per-celltype Seurat objects to use for running limma dream on the long-covid status\
'*differential_expression/lc_differential_expression_limma_parameterised.R*' perform differential gene expression using limma dream on the long-covid status in the multiome gene expression data\
'*differential_expression/lc_create_limma_jobs.sh*' create sbatch jobs for running limma per celltype on the long-covid status in the multiome gene expression data\
'*differential_expression/mo_differential_expression_limma_parameterised.R*' perform differential gene expression using limma dream on the stimulation status in the multiome gene expression data\
'*differential_expression/mo_plot_de.R*' plot differential gene expression results

### gene set enrichment
'*differential_expression/mo_gse_enrichr.R*' create per-celltype Seurat objects from the NC2022 study which also used Candida stimulation\

### differential protein expression analysis
'*gene_set_enrichment/lc_covid_olink.R*' perform differential protein analysis with a linear regression model on the long-covid status\

### cpeaks-based ATAC data processing
'*cpeaks_preprocess_samples/mo_create_cpeaks_jobs.sh*' use the fragments from cellranger, and transform them to cpeaks reference mapped fragments\
'*cpeaks_preprocess_samples/mo_get_credible_atac_barcodes.R*' get the barcodes from archR that passed QC\
'*cpeaks_preprocess_samples/mo_get_arc_metadata.R*' get the cellranger-arc metadata from the summary csvs\
'*cpeaks_preprocess_samples/mo_cpeaks_to_signac.R*' read cpeaks reference-based fragments and read them into Signac\
'*cpeaks_preprocess_samples/mo_cpeaks_to_celltypes.R*' create multimodal object per cell type, and write beds of the peak data\
'*cpeaks_preprocess_samples/mo_cpeaks_to_celltypes.R*' write beds of the peak data, seperated by condition and cell type\
'*cpeaks_preprocess_samples/mo_subset_fragment_files.ipynb*' subset a fragments file based on the annotations of the barcodes\
'*cpeaks_preprocess_samples/mo_create_multimodal_object_per_celltype.R*' merge the Seurat countdata and the Signac ATAC data into multimodal objects per cell type\
'*cpeaks_preprocess_samples/mo_multimodal_clustering.R*' perform PCA, knn-clustering and 2d UMAP using both modalities separately, and together

### ATAC data versus other datasets
'*atac_replication/mo_download_ihec.sh*' download data from ihec to compare ATAC data from this dataset against other ones\
'*atac_replication/mo_ihec_bigbed_to_bed.sh*' convert bigbed files to bed files from ihec

### differentially accessible region detection
'*differential_accessibility/mo_create_per_celltype_signac_objects.R*' create a Signac object with ATAC data per cell type\
'*differential_accessibility/mo_anonymize_signac_objects.R*' remap sample IDs to other identifiers so data can be moved to other clusters\
'*differential_accessibility/lc_differential_accessibility_limma_parameterised.R*' check differential accessible regions using limma on the long-covid status\
'*differential_accessibility/mo_differential_accessibility_limma_parameterised.R*' check differential accessible regions using limma on the stimulation status\
'*differential_accessibility/mo_create_limma_dar_jobs.sh*' create sbatch jobs to perform limma for DAR identification\
'*differential_accessibility/mo_differential_accessibility_kimma_parameterised.R*' check differential accessible regions using kimma on the stimulation status\
'*differential_accessibility/mo_differential_accessibility_limma_add_perm_fdr.R*' add permutation-based FDR to the DAR identification step\
'*differential_accessibility/mo_annotate_limma_dar_output.R*' add closest gene annotation to DAR output\
'*differential_accessibility/mo_analyse_dar_output.Rmd*' analyse the DAR output from limma

### eQTL mapping
'*qtl/eqtl/LIMIX/mo_create_limix_qtl_input.R*' create input for LIMIX eQTL mapping

### eQTL mapping
'*qtl/caqtl/mo_create_limix_chromatin_input.R*' create input for LIMIX caQTL mapping



