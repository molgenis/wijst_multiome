# Wijst Multi-ome study

This is the github repository of the Multiome (snRNA+snATAC) study containing untreated (UT) and Candida (24hCA) stimulated PBMCs

## Software

Here we will list the software used to generate the data

Cellranger 7.1.0\
Cellranger-arc 2.0.0\
R 4.4.2\
plink2-20230707\
CellBender v3\
Seurat 5.2.1\
Signac 1.14.0\
ArchR 1.0.3\
Souporcell v2\
Demuxlet v2

## Custom code

python/jupyter scripts have environment files supplied
R software was run through this singularity container: https://github.com/royoelen/single-cell-container-server/releases/tag/v2.0.1

Below we will outline the steps taken to process and analyse the data


### alignment to reference genome
'*alignment/mo_create_cellranger_jobs_batch1.sh*'   create SLURM jobs to align the RNA-seq data to the b38 human reference\
'*alignment/mo_create_multiome_csvs_batch1.sh*'     create SLURM jobs to do a joint align of the RNA-seq and ATAC-seq data to the b38 human reference\
'*alignment/mo_create_cellranger_arc_jobs_batch1.sh*'   create SLURM jobs to align the RNA+ATAC data to the b38 human reference\
'*alignment/mo_download_websummaries.sh*' download the  CellRanger web summary files\
'*alignment/mo_get_cellranger_summaries.R*'     created an Excel sheet describing the CellRanger results, from the CellRanger output directories


### genotypeing
'*genotyping/mo_collect_previously_genotyped_individuals.R*'    collect data on which individuals were previously genotyped\
'*genotyping/mo_get_previous_genotypes.R*'  create list of previously genotyped individuals, create sample sheets, create age/sex metadata file\
'*genotyping/mo_preprocess_previous_genotypes.sh*'  subset previously genotyped individuals genotype data, and convert to format for imputation\
'*genotyping/mo_create_ugli_psam.R*'    create psam file required for imputation, for the previously genotyped individuals\
'*genotyping/mo_create_batch2_psam.R*'  create psam file required for imputation, for the second batch of individuals\
'*genotyping/PreImputation_ugli.yaml*'  the imputation config used for imputing the previously genotyped UGLI individuals\
'*genotyping/PreImputation_mo_batch2.yaml*'     the imputation config used for imputing the second batch of new individuals


### ambient RNA correction
'*ambient_rna_correction/mo_cellbend_rna.sh*' create sbatch jobs to do ambient RNA correction using CellBender on the RNA data


### demultiplexing and doublet detection
'*demultiplexing/mo_subset_for_genotype_correlations.sh*'   script to subset genotypes for each 10x experiment, using the individuals in present in that experiment\
'*demultiplexing/mo_subset_gex_bams_snp_regions.sh*'    script to create jobs to subset the gene expression alignment files, to only have reads overlapping variants we genotyped\
'*demultiplexing/mo_create_popscle_sort_vcfs_jobs.sh*' s    cript to create jobs to sort the per-lane VCF files by the order of chromosomes in the alignment files\
'*demultiplexing/mo_create_demuxlet_jobs.sh*'   script to create jobs to run Demuxlet on the gene expression alignments\
'*demultiplexing/mo_create_souping_samples_jobs_externalbarcodes.sh*'   script to create jobs that perform Souporcell on each 10x lane while allowing external barcode files\
'*demultiplexing/mo_correlate_genotypes.R*'     correlate the Souporcell cluster genotypes to the genotypes generated for the individuals, to do sample assignment\
'*demultiplexing/mo_plot_demultiplexing_assignments.Rmd*'   plot the demultiplexing assignments\
'*demultiplexing/mo_sample_assignment.R*'   do sample assignment on the Seurat object, based on the correlated souporcell genotypes\
'*demultiplexing/mo_get_rna_qced_barcodes.R*'   extract the valid barcodes in the Seurat object


### scanpy preprocess
'*scanpy_preprocess/mo_h5_to_scanpy_objects.ipynb*'     read the Seurat-incompatible h5 CellBender outputs into Scanpy H5AD objects\
'*scanpy_preprocess/mo_deconstruct_scanpy_objects.ipynb*'   deconstruct the H5AD Scanpy objects into Seurat compatible barcodes/features/matrix files


### seurat preprocess
'*seurat_preprocess/mo_lane_to_seurat.R*'   read the barcodes/features/matrix files into Seurat objects for the RNA counts


### cell type assignment
'*cell_type_annotation/mo_create_10x_mo_reference.R*'   create a reference Azimuth file from the 10x PBMC multi-ome reference dataset, which will also add the ATAC data using signac\
'*cell_type_annotation/mo_azimuth_reference_map_mo.R*'  do an Azimuth reference mapping of each 10x lane onto the 10x PMBC reference\
'*cell_type_annotation/mo_azimuth_reference_map_pbmc.R*'    do an Azimuth reference mapping of each 10x lane onto the Azimuth COVID-19 PBMC dataset\
'*cell_type_annotation/mo_check_celltype_assignments.Rmd*'  check the cell type assignments that were done


### merge seurat objects
'*seurat_preprocess/mo_merge_seurat_objects.R*'     merge the celltyped Seurat files for each lane\
'*seurat_preprocess/mo_filter_qc.R*'    do QC on the count data\
'*seurat_preprocess/mo_normalize_and_cluster.R*'    perform SCT normalization, PCA, knn-clustering and 2d UMAP\
'*seurat_preprocess/mo_add_rna_metadata.R*'     add covid assignments to seurat metadata\
'*seurat_preprocess/lc_create_per_celltype_objects.R*'  create a Seurat object for each cell type, with only UT or both UT and 24hCA\
'*seurat_preprocess/mo_create_l2_objects.R*'    create a Seurat object per L2 Azimuth covid19 annotation


### differential celltype abundance analysis
'*differential_abundance/mo_plot_celltype_proportions.Rmd*'     plot the cell type proportions present in the dataset\
'*differential_proportion/mo_differential_proportion_analysis.R*'   use Speckle to do differential abundance analysis


### differential gene expression analysis
'*differential_expression/1m_create_per_celltype_objects.R*'    create per-celltype Seurat objects from the NC2022 study which also used Candida stimulation\
'*differential_expression/1m_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the stimulation status in the NC2022 data\
'*differential_expression/lc_create_per_celltype_objects.R*'    create per-celltype Seurat objects to use for running limma dream on the long-covid status\
'*differential_expression/lc_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the long-covid status in the multiome gene expression data\
'*differential_expression/lc_create_limma_jobs.sh*'     create sbatch jobs for running limma per celltype on the long-covid status in the multiome gene expression data\
'*differential_expression/mo_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the stimulation status in the multiome gene expression data\
'*differential_expression/mo_plot_de.R*'    plot differential gene expression results


### gene set enrichment
'*gene_set_enrichment/mo_gse_enrichr.R*'    perform gene set enrichment on the differential gene expression results


### differential protein expression analysis
'*differential_protein/lc_covid_olink.R*'   perform differential protein analysis with a linear regression model on the long-covid status


### cpeaks-based ATAC data processing
'*cpeaks_preprocess_samples/mo_get_credible_atac_barcodes.R*'   get the barcodes from archR that passed QC\
'*cpeaks_preprocess_samples/mo_get_arc_metadata.R*' get the cellranger-arc metadata from the summary csvs\
'*cpeaks_preprocess_samples/mo_create_cpeaks_jobs.sh*'  use the fragments from cellranger, and transform them to cpeaks reference mapped fragments\
'*cpeaks_preprocess_samples/mo_cpeaks_to_signac.R*' read cpeaks reference-based fragments and read them into Signac\
'*cpeaks_preprocess_samples/mo_cpeaks_to_celltypes.R*'  create multimodal object per cell type, and write beds of the peak data\
'*cpeaks_preprocess_samples/mo_cpeaks_to_celltypes.R*'  write beds of the peak data, seperated by condition and cell type\
'*cpeaks_preprocess_samples/mo_subset_fragment_files.ipynb*'    subset a fragments file based on the annotations of the barcodes, so check for biases\
'*cpeaks_preprocess_samples/mo_create_multimodal_object_per_celltype.R*'    merge the Seurat countdata and the Signac ATAC data into multimodal objects per cell type\
'*cpeaks_preprocess_samples/mo_multimodal_clustering.R*'    perform PCA, knn-clustering and 2d UMAP using both modalities separately, and together


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


### ATAC data versus other datasets
'*atac_replication/mo_download_ihec.sh*'    download data from ihec to compare ATAC data from this dataset against other ones\
'*atac_replication/mo_ihec_bigbed_to_bed.sh*'   convert bigbed files to bed files from ihec


### differentially accessible region detection
'*differential_accessibility/mo_create_per_celltype_signac_objects.R*'  create a Signac object with ATAC data per cell type\
'*differential_accessibility/mo_anonymize_signac_objects.R*'    remap sample IDs to other identifiers so data can be moved to other clusters\
'*differential_accessibility/lc_differential_accessibility_limma_parameterised.R*'  check differential accessible regions using limma on the long-covid status\
'*differential_accessibility/mo_differential_accessibility_limma_parameterised.R*'  check differential accessible regions using limma on the stimulation status\
'*differential_accessibility/mo_create_limma_dar_jobs.sh*'  create sbatch jobs to perform limma for DAR identification\
'*differential_accessibility/mo_differential_accessibility_limma_add_perm_fdr.R*'   add permutation-based FDR to the DAR identification step\
'*differential_accessibility/mo_annotate_limma_dar_output.R*'   add closest gene annotation to DAR output\
'*differential_accessibility/mo_analyse_dar_output.Rmd*'    analyse the DAR output from limma
'*differential_accessibility/mo_compare_imputed_nonimputed.ipynb*'  compare the imputed vs the non-imputed count matrices\
'*differential_accessibility/mo_check_pycistopic_imputations.R*'  compare the imputed vs the non-imputed count matrices
'*differential_accessibility/mo_export_imputed_pycistopic_matrices.ipynb*'  perform imputation in pycistopic and extract the imputed count matrices\
'*differential_accessibility/mo_extract_topic_memberships.ipynb*'  extract the topic membership of each cell from pycistopic and export that to a table\
'*differential_accessibility/mo_differential_accessibility_topics.R*'  use limma to identify DARs across topics\
'*differential_accessibility/mo_differential_accessibility_imputed.R*'  use limma to identify DARs across topics, with the imputed count matrices\
'*differential_accessibility/mo_limma_to_bed.R*'  export the DAR output from limma to bed format\
'*differential_accessibility/mo_plot_dars_topics.Rmd*'  plot the sharedness of the DARs across the topics\
'*differential_accessibility/mo_check_dar_methods.Rmd*'  plot the sharedness of the DARs across different methods


### CRE detection
'*cre_detection/mo_scenicplus_env.yml*' environment with packages used for CRE detection

#### data peparation
'*cre_detection/data_preparation/mo_deconstruct_signac_objects.R'* deconstruct the Signac objects into matrix/features/barcodes, to merge for pycistopic
'*cre_detection/data_preparation/mo_deconstruct_rna_objects.R'* deconstruct the Seurat objects into matrix/features/barcodes, to use in the SCENIC+ step
'*cre_detection/data_preparation/mo_parts_to_scanpy.py'*    create scanpy object from deconstructed Seurat object
'*cre_detection/data_preparation/mo_get_outer_join_filtered_chromatin_regions.R'*   get all the regions that are represented in .1% cells in any cell type
'*cre_detection/data_preparation/mo_create_pycistopic_annotations.R'*   create annotations for regions to use in pycistopic
'*cre_detection/data_preparation/mo_merge_chunked_mtx_files.py'*    merge chunked mtx files into one file

#### topic modelling
'*cre_detection/pycistopic/mo_create_pycistopic_object.ipynb'*  create pycistopic object\
'*cre_detection/pycistopic/mo_pycistopic_model_topics.ipynb'*   perform topic modelling\
'*cre_detection/pycistopic/mo_pycistopic_binarize_topics.ipynb'*    binarize topics\
'*cre_detection/pycistopic/mo_pycistopic_impute_accessibility.ipynb'*   impute accessibility matrices\
'*cre_detection/pycistopic/mo_pycistopic_normalize_accessibility.ipynb'*   normalize accessibility matrices\
'*cre_detection/pycistopic/mo_pycistopic_normalize_accessibilityfind_variable_features.ipynb'*   find DARs that vary\
'*cre_detection/pycistopic/mo_pycistopic_export_topic_contributions.ipynb'*   export contributions of cells and regions to topics to tables\
'*cre_detection/pycistopic/mo_region_topics_to_beds.R'*   export topic membership of regions to bed files

#### DAR identification
'*cre_detection/dar_identification/mo_identify_dars_wilcoxon.py'*  perform DAR detection between topics in pycistopic using the wilcoxon rank sum test\
'*cre_detection/dar_identification/mo_identify_dars_celltypes_wilcoxon.py'*  perform DAR detection between cell type in pycistopic using the wilcoxon rank sum test\
'*cre_detection/dar_identification/mo_extract_dar_outputs.ipynb'*  convert binary DAR detection output into tsv format\
'*cre_detection/dar_identification/mo_check_dar_numbers.ipynb'*  check DAR and topic membership against available metadata

#### eRegulon detection
'*cre_detection/pycistarget/mo_cistarget.ipynb*'  run pycistarget and DEM to identify overrepresented motifs in the DARs

'*cre_detection/scenicplus/scenicplus_config.yaml*'  config for running scenic+ pipeline after setting up all inputs

#### output analysis
'*cre_detection/output_analysis/mo_calculate_eregulon_enrichments.ipynb*'  calculate enriched eRegulons for specific metadata variables\
'*cre_detection/output_analysis/mo_plot_scenic_output.Rmd*'  plot eRegulons, genes and regions found in SCENIC, and compare to smaller 10x dataset\
'*cre_detection/output_analysis/mo_scenic_add_region_info_to_cres.R*'  add information regarding QTLs to SCENIC output\
'*cre_detection/output_analysis/mo_scenic_vs_granie_claringbould.Rmd*'  compare SCENIC output to macrophage dataset\
'*cre_detection/output_analysis/mo_cre_method_outputs_comparison.R*'  compare SCENIC to pseudobulk and regression-based models


### eQTL mapping
'*qtl/eqtl/LIMIX/mo_create_limix_qtl_input.R*'  create input for LIMIX eQTL mapping\
'*qtl/eqtl/mo_create_n_cellss_expressed_tables.R*'  create table of number of non-zero nuclei per gene and donor\
'*qtl/eqtl/LIMIX/mo_annotation_to_chunking_file.R*' create chunking file for eQTL mapping in LIMIX\
'*qtl/eqtl/LIMIX/limix_qtl.smk*'    LIMIX snakemake file to do eQTL mapping\
'*qtl/eqtl/LIMIX/mo_qtl_template.yaml*' LIMIX configuration file to do eQTL mapping


### caQTL mapping
'*qtl/caqtl/mo_create_limix_chromatin_input.R*'     create input for LIMIX caQTL mapping\
'*qtl/caqtl/mo_create_n_cellss_accessible_tables.R*'  create table of number of non-zero nuclei per region and donor\
'*qtl/caqtl/mo_create_caqtl_feature_filters.R*'   create lists of features to test for caQTL mapping\
'*qtl/caqtl/mo_create_caqtl_feature_filters.R*'   create lists of features to test for caQTL mapping\
'*qtl/caqtl/mo_caqtl_lcl_replication.Rmd*'  plot replication of LCL caQTLs in our B caQTL output


### QTL results
'*qtl/mo_filter_down_significant_results.R*'    subset QTL output by significance\
'*qtl/mo_regress_pcs_qtlinput.ipynb*'   regress principal components out of QTL input matrices so they plot more like they are modelled\
'*qtl/mo_plot_qtls.ipynb*'    plot how the QTLs and CREs look\
'*qtl/mo_qtl_caqtl_eqtl_overlaps.ipynb*'    plot overlapping caQTLs and eQTLs\
'*qtl/mo_plot_qtl_numbers.Rmd*'     plot the number of QTLs\
'*qtl/mo_get_finemapped_variants.R*'     extract finemapped and non-finemapped eQTLs from sc-eQTLgen to compare the variants\
'*qtl/mo_plot_independent_qtls.Rmd*'     extract and plot overlapping and colocalizing QTLs


### interaction-QTL mapping
'*qtl/interaction_eqtl/mo_get_significant_variant_feature_pairs.R*'    get significant variant-feature pairs from the QTL mappings\
'*qtl/interaction_eqtl/mo_create_limix_interaction_qtl_input.R*'    create interaction-eQTL input files\
'*qtl/interaction_eqtl/limix_interactions.smk*'    LIMIX-QTL interaction snakemake file\
'*qtl/interaction_eqtl/mo_interaction_template.yaml*'    LIMIX-QTL interaction yaml file for interaction-eQTLs\
'*qtl/interaction_caqtl/mo_create_limix_chromatin_interaction_input.R*'    create interaction-caQTL input files\
'*qtl/interaction_caqtl/mo_interaction_caqtls.yaml*'    LIMIX-QTL interaction yaml file for interaction-caQTL


### QTL mediation
'*qtl/mediation/mo_create_mediation_confinements.R*'    create confinement files for eQTL-by-caQTL mediation analyses\
'*qtl/mediation/mo_create_mediation_metadata.R*'    create metadata for eQTL-by-caQTL mediation analyses\
'*qtl/mediation/mo_perform_qtl_mediation_analysis.R*'    run eQTL-by-caQTL mediaton analysis script\
'*qtl/mediation/mo_create_mediation_jobs.sh*'    create run eQTL-by-caQTL mediaton analysis jobs


### QTL finemapping
'*qtl/finemapping/mo_finemap_qtls.R*'   perform finemapping on QTL summary statistics coming from LIMIX-QTL\
'*qtl/finemapping/lpmcv2_format_finemapping*'   convert binary .rds finemapping results into tsv files\
'*qtl/finemapping/mo_create_finemap_jobs.sh*'   create sbatch jobs to do finemapping


### QTL utility scripts
'*qtl/mo_regress_qtlinputs.py*'    regress PCs out of QTL input files\
'*qtl/mo_eigenmt_correct_limix_qtls.R*'    perform eigenMT MTC on QTL outputs


### QTL CRE replication
'*qtl/cre_qtl/mo_creqtl_env.yml*'  yaml for creating conda environment to do CRE replication
'*qtl/cre_qtl/mo_split_sample_and_celltype.R*'  split Seurat object into separate matrices/features/barcodes for each sample
'*qtl/cre_qtl/mo_split_sample_and_celltype_each.sh*'    split each Seurat object
'*qtl/cre_qtl/mo_calculate_atac_rna_betas.py*'  calculate the scaled beta+se between accessibility and expression using a binomial model for a combination of an accessiblity and an expression matrix
'*qtl/cre_qtl/mo_create_beta_calculation_jobs.sh*'  create jobs to calculate beta+se of accessibility and expression for each sample
'*qtl/cre_qtl/mo_aggregated_creqtl_inputs.sh*'  aggregate per-sample outputs of beta calculation jobs
'*qtl/cre_qtl/mo_meta_analyse_creqtl_cres*' meta-analyse betas and ses calculated and aggregated in previous steps to get to significant region-gene pairs
'*qtl/cre_qtl/mo_plot_replicating_cres.ipynb*'  plot properties of region-gene pairs that were overlapping eQTLs/caQTLs and replicate as CREs
'*qtl/cre_qtl/mo_plot_creqtls.Rmd*' plot proportion of region-gene pairs that were overlapping eQTLs/caQTLs and replicate as CREs


### cell type composition GWAS
'*ctc_gwas/mo_ctc_gwas_wg2_step2.sh*'    for the cell type composition GWAS, perform step 2 of sc-eQTLgen WG2 to get consortium-compatible cell types\
'*ctc_gwas/mo_ctc_gwas_wg2_step3.sh*'    for the cell type composition GWAS, perform step 3 of sc-eQTLgen WG2 to get consortium-compatible cell types\
'*ctc_gwas/mo_ctc_gwas_wg2_step4.sh*'    for the cell type composition GWAS, perform step 4 of sc-eQTLgen WG2 to get consortium-compatible cell types


## scripts no longer used
'*demultiplexing/mo_run_scrublet.py*'   python script to run Scrublet on a CellBender corrected 10x lane\
'*demultiplexing/mo_test_scrublet.ipynb*'   jupyter notebook to test Scrublet on individual 10x lanes\
'*differential_accessibility/mo_differential_accessibility_kimma_parameterised.R*'  check differential accessible regions using kimma on the stimulation status



## License

The code availabe in this repository is available under GPLv2 License:
https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
