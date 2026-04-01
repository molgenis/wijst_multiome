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
'*alignment/mo_get_cellranger_summaries.R*'     created an Excel sheet describing the CellRanger results, from the CellRanger output directories\
'*alignment/mo_copy_trailing_cellranger_data.sh*'     copy cellranger folder back to tmp while keeping folder structure


### genotyping
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
'*demultiplexing/mo_subset_for_genotype_correlations.sh*'   script to subset genotypes for each 10x experiment, using the individuals in present in that experiment as job\
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
'*seurat_preprocess/mo_lane_to_seurat.R*'   read the barcodes/features/matrix files into Seurat objects for the RNA counts\
'*seurat_preprocess/mo_lane_to_seuratv4.R*'   read the barcodes/features/matrix files into Seurat objects for the RNA counts, using Seurat v4 format


### cell type assignment
'*cell_type_annotation/mo_create_10x_mo_reference.R*'   create a reference Azimuth file from the 10x PBMC multi-ome reference dataset, which will also add the ATAC data using signac\
'*cell_type_annotation/mo_azimuth_reference_map_mo.R*'  do an Azimuth reference mapping of each 10x lane onto the 10x PMBC reference\
'*cell_type_annotation/mo_azimuth_reference_map_pbmc.R*'    do an Azimuth reference mapping of each 10x lane onto the Azimuth COVID-19 PBMC dataset\
'*cell_type_annotation/mo_check_celltype_assignments.Rmd*'  check the cell type assignments that were done\
'*cell_type_annotation/mo_plot_celltype_markers.R*'    check marker expression for cells and the cell types they were assigned to\
'*cell_type_annotation/mo_plot_umap_celltypes.Rmd*'  plot the assignments of cells onto the UMAP


### merge seurat objects
'*seurat_preprocess/mo_merge_seurat_objects.R*'     merge the celltyped Seurat files for each lane\
'*seurat_preprocess/mo_filter_qc.R*'    do QC on the count data\
'*seurat_preprocess/mo_normalize_and_cluster.R*'    perform SCT normalization, PCA, knn-clustering and 2d UMAP\
'*seurat_preprocess/mo_add_rna_metadata.R*'     add covid assignments to seurat metadata\
'*seurat_preprocess/mo_create_l2_objects.R*'    create a Seurat object per L2 Azimuth covid19 annotation


### differential celltype abundance analysis
'*differential_proportion/mo_plot_celltype_proportions.Rmd*'     plot the cell type proportions present in the dataset\
'*differential_proportion/mo_differential_proportion_analysis.R*'   use Speckle to do differential abundance analysis\
'*differential_proportion/mo_plot_cell_numbers.Rmd*'   plot cell numbers recovered


### differential gene expression analysis
'*differential_expression/1m_create_per_celltype_objects.R*'    create per-celltype Seurat objects from the NC2022 study which also used Candida stimulation\
'*differential_expression/1m_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the stimulation status in the NC2022 data\
'*differential_expression/mo_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the stimulation status in the multiome gene expression data\
'*differential_expression/mo_plot_de.R*'    plot differential gene expression results\
'*differential_expression/mo_create_limix_qtl_input_l2.R*'    create per-celltype Seurat objects at L2 resolution\
'*differential_expression/mo_merge_dar_de_tables.R.R*'    merge DAR and DE tables


### gene set enrichment
'*gene_set_enrichment/mo_gse_enrichr.R*'    perform gene set enrichment on the differential gene expression results


### round ATAC fragments
'*fragment_rounding/round_fragments.py*'   round ATAC fragments


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
'*differential_accessibility/mo_differential_accessibility_limma_parameterised.R*'  check differential accessible regions using limma on the stimulation status\
'*differential_accessibility/mo_create_limma_dar_jobs.sh*'  create sbatch jobs to perform limma for DAR identification\
'*differential_accessibility/mo_create_limma_dar_topic_jobs.sh*'  create sbatch jobs to perform limma for DAR identification across topics\
'*differential_accessibility/mo_differential_accessibility_limma_add_perm_fdr.R*'   add permutation-based FDR to the DAR identification step\
'*differential_accessibility/mo_annotate_limma_dar_output.R*'   add closest gene annotation to DAR output\
'*differential_accessibility/mo_analyse_dar_output.Rmd*'    analyse the DAR output from limma
'*differential_accessibility/mo_compare_imputed_nonimputed.ipynb*'  compare the imputed vs the non-imputed count matrices\
'*differential_accessibility/mo_check_pycistopic_imputations.R*'  compare the imputed vs the non-imputed count matrices
'*differential_accessibility/mo_export_imputed_pycistopic_matrices.ipynb*'  perform imputation in pycistopic and extract the imputed count matrices\
'*differential_accessibility/mo_extract_topic_memberships.ipynb*'  extract the topic membership of each cell from pycistopic and export that to a table\
'*differential_accessibility/mo_differential_accessibility_topics.R*'  use limma to identify DARs across topics\
'*differential_accessibility/mo_differential_accessibility_limma_add_perm_fdr.R*'  add permutation based FDR to DARs\
'*differential_accessibility/mo_differential_accessibility_imputed.R*'  use limma to identify DARs across topics, with the imputed count matrices\
'*differential_accessibility/mo_limma_to_bed.R*'  export the DAR output from limma to bed format\
'*differential_accessibility/mo_plot_dars_topics.Rmd*'  plot the sharedness of the DARs across the topics\
'*differential_accessibility/mo_check_dar_methods.Rmd*'  plot the sharedness of the DARs across different methods\
'*differential_accessibility/mo_pycistopic_dars.ipynb*'  extract DARs from monocyte pycistopic run\
'*differential_accessibility/mo_pycistopic_topic_modeling_all_mono.py*'  perform pycistopic on only monocytes


### gene score calculation
'*/gene_score_calculation/mo_check_signac_peaks_celltype_lane.R*'   add gene annotations to signac peaks\
'*/gene_score_calculation/mo_check_signac_peaks_celltype_lane_jobs.sh*'   run jobs to add gene annotations to signac peaks


### CRE detection

#### comparisons
'*cre_detection/comparisons/mo_check_opposite_cre_effects.R'* check CRE-gene pairs that have opposite effects in different methods\
'*cre_detection/comparisons/mo_cre_method_outputs_comparison.R*'  compare SCENIC to pseudobulk and regression-based models\
'*cre_detection/comparisons/mo_create_creqtl_plots.R*'  plot CRE-i-eQTL effects\
'*cre_detection/comparisons/mo_do_full_cre_overlap_check.R*'  check overlap and concordance of TF-CRE-gene sets across CRE detection methods\
'*cre_detection/comparisons/mo_get_encode_cre_genes_to_cpeaks.R*'  add encode region-gene link information to cpeaks\
'*cre_detection/comparisons/mo_merge_cpeaks_with_screenv4.R*'  add encode screen v4 information to cpeaks defined regions\
'*cre_detection/comparisons/mo_overlap_qtl_with_cres.R*'  add information on cpeaks regions to where variants might be located in


#### LIMIX/single-cell method
'*cre_detection/limix_hybrid/mo_create_frac_exp_acc_files.R*'  create annotations for how prevalent regions are accessible across samples\
'*cre_detection/limix_hybrid/mo_create_frac_exp_files.R*'  create annotations for how prevalent genes are expressed are accessible across samples\
'*cre_detection/limix_hybrid/mo_create_hybrid_cre_inputs_atac.R*'  create chunked accessiblity tables for all data at once\
'*cre_detection/limix_hybrid/mo_create_hybrid_cre_inputs_rna.R*'  create chunked expression tables for all data at once\
'*cre_detection/limix_hybrid/mo_create_hybrid_cre_cov_matrix.R*'  create binary covariate and kinship data per cell\
'*cre_detection/limix_hybrid/mo_merge_hybrid_cres.R*'  merge the chunked CRE mappings\
'*cre_detection/limix_hybrid/mo_sample_hybrid_cre_inputs.R*'  randomly sample CRE inputs to check stability


##### LIMIX output analysis
'*cre_detection/limix_hybrid/output_analysis/mo_hybrid_vs_hic_comparison.R*'  overlap LIMIX CRE-gene links with Hi-C data from encode\
'*cre_detection/limix_hybrid/output_analysis/mo_limix_cre_celltype_replication.Rmd.R*'  plot replication of LIMIX CRE-gene links across cell types\
'*cre_detection/limix_hybrid/output_analysis/mo_limix_hybrid_vs_reunion.R.R*'  overlap LIMIX CRE-gene links with CRE-gene pairs in REUNION paper\
'*cre_detection/limix_hybrid/output_analysis/mo_plot_sccres.R*'  plot LIMIX CRE-gene links


##### single-cell CRE interaction analysis
'*cre_detection/limix_hybrid/interaction_analysis/mo_hybrid_cre_interaction_overlaps.R*'  plot interaction-eQTL at single-cell level with TF or ATAC as interaction terms overlaps across methods
'*cre_detection/limix_hybrid/interaction_analysis/mo_hybrid_cre_plot_interaction.R*'  plot specific interaction-eQTLs from chunk\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/mo_hybrid_cre_interaction.R*'  perform single-cell interaction-eQTL analysis\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/mo_hybrid_interactions.smk*'  merge single-cell interaction-eQTL analysis chunks\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/mo_hybrid_cre_interaction_merge.R*'  snakemake to run single-cell interaction-eQTL analysis chunks\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/mo_plot_qtl_method_overlap.Rmd*'  overlap interaction-eQTLs with eQTL/caQTL/SCENIC+\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/region_interaction/mo_cre_sccre_interaction_confinement_r2g.R*'  create CRE-i-eQTL confinement file\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/region_interaction/mo_hybrid_interactions_regions_template.yaml*'  yaml configuration for CRE-i-eQTL run\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/mo_add_pseudobulked_tf_activities.R*'  calculate pseudobulked TF activities\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/mo_calc_auc_eregulons.R*'  calculate gene AUC based TF activities\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/mo_compare_tf_to_region_ieqtls.Rmd*'  compare TF-i-eQTLs to CRE-i-eQTLs\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/mo_cre_sccre_interaction_confinement.R*'  create TF-i-eQTLs confinement\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/mo_hybrid_interations_template.yaml*'  yaml configuration for TF-i-eQTLs run\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/mo_recalc_tf_i_eqtl_auc_eregulons.R*'  recalculate TF-i-eQTL TF activities excluding TF-i-eGenes\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/mo_hybrid_interations_template_replication.yaml*'  yaml configuration for TF-i-eQTLs validation run that removed TF-i-eGenes from TF activities\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/viacheslav_replication/mo_hybrid_interations_template_viacheslavrep.yaml*'  yaml configuration for TF-i-eQTLs replication from Viacheslav paper\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/coeqtl_replication/mo_sccre_coeqtl_tf_interaction_confinement.R*'  create confinement for replicating co-eQTL TF-gene pairs\
'*cre_detection/limix_hybrid/interaction_analysis/interaction_like/tf_interaction/coeqtl_replication/mo_hybrid_interations_template_coeqtlrep.yaml*'  yaml configuration for TF-i-eQTLs replication from coeqtl paper\
'*cre_detection/limix_hybrid/interaction_analysis/coeqtl_like/mo_hybrid_cre_coeqtl_replication.R*'  perform co-eQTL style TF-i-eQTL analysis for chunk\
'*cre_detection/limix_hybrid/interaction_analysis/coeqtl_like/mo_hybrid_interactions_coeqtl.smk*'  perform co-eQTL style TF-i-eQTL analysis across chunks\
'*cre_detection/limix_hybrid/interaction_analysis/coeqtl_like/mo_hybrid_interations_coeqtl_template.yaml*'  configuration to perform co-eQTL style TF-i-eQTL analysis across chunks

#### PyCistop/SCENIC+ method

'*cre_detection/pycistopic_scenic/mo_scenicplus_env.yml'* yaml to create conda environment for pyscistopic and SCENIC+

##### data peparation
'*cre_detection/pycistopic_scenic/data_preparation/mo_deconstruct_signac_objects.R'* deconstruct the Signac objects into matrix/features/barcodes, to merge for pycistopic
'*cre_detection/pycistopic_scenic/data_preparation/mo_deconstruct_rna_objects.R'* deconstruct the Seurat objects into matrix/features/barcodes, to use in the SCENIC+ step
'*cre_detection/pycistopic_scenic/data_preparation/mo_parts_to_scanpy.py'*    create scanpy object from deconstructed Seurat object
'*cre_detection/pycistopic_scenic/data_preparation/mo_get_outer_join_filtered_chromatin_regions.R'*   get all the regions that are represented in .1% cells in any cell type
'*cre_detection/pycistopic_scenic/data_preparation/mo_create_pycistopic_annotations.R'*   create annotations for regions to use in pycistopic
'*cre_detection/pycistopic_scenic/data_preparation/mo_merge_chunked_mtx_files.py'*    merge chunked mtx files into one file


##### topic modelling
'*cre_detection/pycistopic_scenic/pycistopic/mo_create_pycistopic_object.ipynb'*  create pycistopic object\
'*cre_detection/pycistopic_scenic/pycistopic/mo_pycistopic_model_topics.ipynb'*   perform topic modelling\
'*cre_detection/pycistopic_scenic/pycistopic/mo_pycistopic_binarize_topics.ipynb'*    binarize topics\
'*cre_detection/pycistopic_scenic/pycistopic/mo_pycistopic_impute_accessibility.ipynb'*   impute accessibility matrices\
'*cre_detection/pycistopic_scenic/pycistopic/mo_pycistopic_impute_accessibility.py'*   impute accessibility matrices\
'*cre_detection/pycistopic_scenic/pycistopic/mo_pycistopic_normalize_accessibility.ipynb'*   normalize accessibility matrices\
'*cre_detection/pycistopic_scenic/pycistopic/mo_pycistopic_find_variable_features.ipynb'*   find DARs that vary\
'*cre_detection/pycistopic_scenic/pycistopic/mo_pycistopic_export_topic_contributions.ipynb'*   export contributions of cells and regions to topics to tables\
'*cre_detection/pycistopic_scenic/pycistopic/mo_region_topics_to_beds.R'*   export topic membership of regions to bed files\
'*cre_detection/pycistopic_scenic/pycistopic/mo_check_imputation_amounts.ipynb'*   check how much data is imputed


##### DAR identification
'*cre_detection/pycistopic_scenic/dar_identification/mo_identify_dars_wilcoxon.py'*  perform DAR detection between topics in pycistopic using the wilcoxon rank sum test\
'*cre_detection/pycistopic_scenic/dar_identification/mo_identify_dars_celltypes_wilcoxon.py'*  perform DAR detection between cell type in pycistopic using the wilcoxon rank sum test\
'*cre_detection/pycistopic_scenic/dar_identification/mo_extract_dar_outputs.ipynb'*  convert binary DAR detection output into tsv format\
'*cre_detection/pycistopic_scenic/dar_identification/mo_check_dar_numbers.ipynb'*  check DAR and topic membership against available metadata


##### cisTarget
'*cre_detection/pycistopic_scenic/pycistarget/mo_cistarget.ipynb*'  run pycistarget and DEM to identify overrepresented motifs in the DARs


##### SCENIC+
'*cre_detection/pycistopic_scenic/scenicplus/scenicplus_config.yaml*'  config for running scenic+ pipeline after setting up all inputs
'*cre_detection/pycistopic_scenic/scenicplus/mo_extract_auc_signature_sets.py*'  config for running scenic+ pipeline after setting up all inputs


##### output analysis
'*cre_detection/pycistopic_scenic/output_analysis/mo_calculate_eregulon_enrichments.ipynb*'  calculate enriched eRegulons for specific metadata variables\
'*cre_detection/pycistopic_scenic/output_analysis/mo_plot_scenic_output.Rmd*'  plot eRegulons, genes and regions found in SCENIC, and compare to smaller 10x dataset\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic_add_region_info_to_cres.R*'  add information regarding QTLs to SCENIC output\
'*cre_detection/pycistopic_scenic/output_analysis/mo_calc_auc_eregulons.R*'  recalculate TF activity based AUC as done in SCENIC+\
'*cre_detection/pycistopic_scenic/output_analysis/mo_check_scenic_celltype_specificity.R*'  check specificity of SCENIC+ TFs/eRegulons based on RSS\
'*cre_detection/pycistopic_scenic/output_analysis/mo_check_scenic_eregulon_vs_iegenes.R*'  check overlap of i-eGenes and TF associated genes\
'*cre_detection/pycistopic_scenic/output_analysis/mo_differential_tf_activity.R*'  check difference of TF activity between conditions and cell types\
'*cre_detection/pycistopic_scenic/output_analysis/mo_export_tf_activity_matrix.ipynb*'  export AUC gene TF activity matrix from SCENIC+\
'*cre_detection/pycistopic_scenic/output_analysis/mo_plot_celltype_tf_activities.Rmd*'  plot TF activity levels across cell types\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic10x_vs_reunion.R*'  check overlap and concordance of TF-gene combinations in SCENIC+ and REUNION, using the public 10x multiome PBMC dataset\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic_vs_reunion.R*'  check overlap and concordance of TF-gene combinations in SCENIC+ and REUNION\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic_vs_granie_claringbould.Rmd*'  compare SCENIC output to macrophage dataset\
'*cre_detection/pycistopic_scenic/output_analysis/mo_check_scenic_vs_eqtlgen.R*'  check overlap and concordance of TF-gene combinations in SCENIC+ and eQTLgen cis-trans pairs\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic_to_sceqtlgen_overlap.R*'  check overlap and concordance of TF-gene combinations in SCENIC+ and sc-eQTLgen colocalizing cis/trans genes\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic_vs_granie_claringbould.Rmd*'  check overlap and concordance of TF-region-gene combinations in SCENIC+ and GRANIE method\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic_vs_hic_comparison.R*'  check overlap and concordance of region-gene combinations in SCENIC+ and HiC data from ENCODE\
'*cre_detection/pycistopic_scenic/output_analysis/mo_scenic_vs_string_comparison.R*'  check overlap and concordance of TF-gene combinations in SCENIC+ and STRING database


### eQTL mapping
'*qtl/eqtl/LIMIX/mo_create_limix_qtl_input.R*'  create input for LIMIX eQTL mapping\
'*qtl/eqtl/LIMIX/mo_annotation_to_chunking_file.R*' create chunking file for eQTL mapping in LIMIX\
'*qtl/eqtl/LIMIX/limix_qtl.smk*'    LIMIX snakemake file to do eQTL mapping\
'*qtl/eqtl/LIMIX/mo_qtl_template.yaml*' LIMIX configuration file to do eQTL mapping\
'*qtl/eqtl/LIMIX/mo_qtl_template.yaml*' LIMIX configuration file to do eQTL mapping in UT\
'*qtl/eqtl/mo_create_n_cellss_expressed_tables.R*'  create table of number of non-zero nuclei per gene and donor\
'*qtl/eqtl/mo_check_eqtl_correlations_celltypes.Rmd*'  plot eQTL replication across cell types\
'*qtl/eqtl/mo_sceqtl_replication.Rmd*'  plot eQTL replication in sc-eQTLgen\
'*qtl/eqtl/mo_get_esnp_enrichment.R*'  check characteristics of eSNPs


### caQTL mapping
'*qtl/caqtl/limix_qtls_confined.smk*'     run confined caQTL mapping using snakemake file\
'*qtl/caqtl/mo_qtl_confined_template.yaml*'     run confined caQTL mapping using yaml file\
'*qtl/caqtl/mo_create_limix_chromatin_input.R*'     create input for LIMIX caQTL mapping\
'*qtl/caqtl/mo_create_n_cellss_accessible_tables.R*'  create table of number of non-zero nuclei per region and donor\
'*qtl/caqtl/mo_create_caqtl_chunking_and_annotations.R*'  create chunking and annotation files for caQTL mapping\
'*qtl/caqtl/mo_create_caqtl_feature_filters.R*'   create lists of features to test for caQTL mapping\
'*qtl/caqtl/mo_check_caqtl_correlations_celltypes.Rmd*'   plot replication of caQTLs across cell types\
'*qtl/caqtl/mo_caqtl_lcl_replication.Rmd*'  plot replication of LCL caQTLs in our B caQTL output


### interaction-QTL mapping
'*qtl/interaction_eqtl/mo_get_significant_variant_feature_pairs.R*'    get significant variant-feature pairs from the QTL mappings\
'*qtl/interaction_eqtl/mo_create_limix_interaction_qtl_input.R*'    create interaction-eQTL input files\
'*qtl/interaction_eqtl/limix_interactions.smk*'    LIMIX-QTL interaction snakemake file\
'*qtl/interaction_eqtl/mo_interaction_template.yaml*'    LIMIX-QTL interaction yaml file for interaction-eQTLs\
'*qtl/interaction_eqtl/mo_compare_ieqtls_vs_non_ieqtls.R*'    compar characteristics of interacting vs not interacting eQTLs\
'*qtl/interaction_eqtl/mo_compare_interaction_qtls_vs_de_or_dar.R*'    compare characteristics of interaction-QTLs vs DE/DAR numbers\
'*qtl/interaction_eqtl/mo_plot_interaction_qtls.ipynb*'    plot interaction-QTLs\
'*qtl/interaction_eqtl/mo_plot_interaction_vs_condition_qtls.Rmd*'    plot interaction-QTLs versus NC2022 data\
'*qtl/interaction_caqtl/mo_create_limix_chromatin_interaction_input.R*'    create interaction-caQTL input files\
'*qtl/interaction_caqtl/mo_interaction_caqtls.yaml*'    LIMIX-QTL interaction yaml file for interaction-caQTL\
'*qtl/interaction_caqtl/mo_compare_icaqtls_vs_non_icaqtls.R*'    compare interacting vs non interacting caQTLs


### QTL mediation
'*qtl/mediation/mo_create_mediation_confinements.R*'    create confinement files for eQTL-by-caQTL mediation analyses\
'*qtl/mediation/mo_create_mediation_metadata.R*'    create metadata for eQTL-by-caQTL mediation analyses\
'*qtl/mediation/mo_perform_qtl_mediation_analysis.R*'    run eQTL-by-caQTL mediaton analysis script\
'*qtl/mediation/mo_create_mediation_jobs.sh*'    create run eQTL-by-caQTL mediaton analysis jobs\
'*qtl/mediation/mo_merge_mediation_results.R*'    merge mediation results


### QTL finemapping
'*qtl/finemapping/mo_finemap_qtls.R*'   perform finemapping on QTL summary statistics coming from LIMIX-QTL\
'*qtl/finemapping/lpmcv2_format_finemapping*'   convert binary .rds finemapping results into tsv files\
'*qtl/finemapping/mo_create_finemap_jobs.sh*'   create sbatch jobs to do finemapping\
'*qtl/finemapping/mo_add_fm_info_to_qtls.R*'   add finemapping info to QTL summary statistics\
'*qtl/finemapping/mo_get_finemapped_variants.R*'   merge and filter finemapped results\
'*qtl/finemapping/mo_plot_independent_qtls.Rmd*'   plot independent effect numbers per feature, based on finemapping


### QTL colocalization
'*qtl/colocalization/mo_annotate_gwas_variants.R*'   add variant annotations as chrom:var:ref:alt format in GWAS\
'*qtl/colocalization/mo_coloc_qtls.R*'   colocalize eQTLs/caQTLs per cell type\
'*qtl/colocalization/mo_coloc_traits_eqtlgen.R*'   colocalize QTLs with eQTLgen GWAS sumstats\
'*qtl/colocalization/mo_do_gwas_colocs.sh*'   colocalize each cell type with GWAS sumstats


### eQTL/caQTL overlap
'*qtl/eqtl_caqtl_overlap/mo_annotate_overlapping_qtls_with_cres.R*'   add CRE info for overlapping caQTL/eQTL pairs\
'*qtl/eqtl_caqtl_overlap/mo_check_overlapping_qtl_effects_per_cs.R*'   add overlapping effects for caQTL/eQTL pairs per credible set\
'*qtl/eqtl_caqtl_overlap/mo_finemapped_eqtl_to_caqtl.R*'   merge finemapping and overlapping caQTL/eQTLs\
'*qtl/eqtl_caqtl_overlap/mo_qtl_caqtl_eqtl_overlaps.ipynb*'   plot overlapping caQTL/eQTL pairs\
'*qtl/eqtl_caqtl_overlap/mo_plot_dual_qtl_examples.R*'   plot examples of overlapping caQTL/eQTL pairs


### QTL CRE replication
'*qtl/cre_qtl/mo_creqtl_env.yml*'  yaml for creating conda environment to do CRE replication\
'*qtl/cre_qtl/mo_split_sample_and_celltype.R*'  split Seurat object into separate matrices/features/barcodes for each sample\
'*qtl/cre_qtl/mo_split_sample_and_celltype_each.sh*'    split each Seurat object\
'*qtl/cre_qtl/mo_calculate_atac_rna_betas.py*'  calculate the scaled beta+se between accessibility and expression using a binomial model for a combination of an accessiblity and an expression matrix\
'*qtl/cre_qtl/mo_create_beta_calculation_jobs.sh*'  create jobs to calculate beta+se of accessibility and expression for each sample\
'*qtl/cre_qtl/mo_create_beta_calculation_jobs.sh*'  create jobs to calculate beta+se of accessibility and expression for each sample and run on CPU\
'*qtl/cre_qtl/mo_aggregated_creqtl_inputs.sh*'  aggregate per-sample outputs of beta calculation jobs\
'*qtl/cre_qtl/mo_meta_analyse_creqtl_cres*' meta-analyse betas and ses calculated and aggregated in previous steps to get to significant region-gene pairs\
'*qtl/cre_qtl/mo_plot_replicating_cres.ipynb*'  plot properties of region-gene pairs that were overlapping eQTLs/caQTLs and replicate as CREs\
'*qtl/cre_qtl/mo_plot_creqtls.Rmd*' plot proportion of region-gene pairs that were overlapping eQTLs/caQTLs and replicate as CREs\
'*qtl/cre_qtl/mo_compare_naive_cres.Rmd*' compare LIMIX single-cell CRE mapping results across cell types


### Transcript Factor QTL mapping
'*qtl/tfqtl/mo_sc_tfqtl.R*' run single-cell TF-QTL analysis\
'*qtl/tfqtl/mo_tf_qtl.smk*' run single-cell TF-QTL analysis using snakemake file\
'*qtl/tfqtl/mo_tf_qtl_template.yaml*' run single-cell TF-QTL analysis using yaml file\


### LD
'*qtl/ld/mo_calc_lds.R*'    calculate LD with sparse matrices between all variants or two sets of variants\
'*qtl/ld/mo_get_ld_pairs.R*'    get the variants that are in LD with one another


### QTL utility scripts
'*qtl/mo_regress_qtlinputs.py*'    regress PCs out of QTL input files\
'*qtl/mo_regress_pcs_qtlinput.ipynb*'    regress PCs out of QTL input files\
'*qtl/mo_eigenmt_correct_limix_qtls.R*'    perform eigenMT MTC on QTL outputs\
'*qtl/mo_get_ncells_analysis.R*'    get number of cells used for generating each pseudobulk\
'*qtl/mo_merge_qtl_tables.R*'    merge QTL result tables for supplements\
'*qtl/mo_plot_pseudobulk_correlations.R*'    plot correlation of pseudobulk gene values against one another\
'*qtl/mo_qtl_variant_to_region.R*'    match cpeaks or screen regions to variants


### QTL results
'*qtl/mo_filter_down_significant_results.R*'    subset QTL output by significance\
'*qtl/mo_regress_pcs_qtlinput.ipynb*'   regress principal components out of QTL input matrices so they plot more like they are modelled\
'*qtl/mo_plot_qtls.ipynb*'    plot how the QTLs and CREs look\
'*qtl/mo_plot_qtls.Rmd*'    plot how the QTLs and CREs look\
'*qtl/mo_qtl_caqtl_eqtl_overlaps.ipynb*'    plot overlapping caQTLs and eQTLs\
'*qtl/mo_plot_qtl_numbers.Rmd*'     plot the number of QTLs\
'*qtl/mo_get_finemapped_variants.R*'     extract finemapped and non-finemapped eQTLs from sc-eQTLgen to compare the variants\
'*qtl/mo_plot_independent_qtls.Rmd*'     extract and plot overlapping and colocalizing QTLs\
'*qtl/mo_annotate_qtls_with_cres.R*'     check QTL tables for overlap with DAR/SCENIC+/openness data\
'*qtl/mo_annotate_overlapping_qtls_with_cres.R*'     check overlapping eQTLs/caQTLs if they are present in SCENIC+ output


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


### EGA
'*ega/mo_create_multiome_tmp_copies_batch1.sh*' copy sequence data to tmp\
'*ega/mo_create_multiome_egacryptor_files_batch1.sh*' create encrypted copies of sequence data\
'*ega/mo_create_multiome_ega_pairing_file_batch1.sh*' create pairing file for encrypted and uncrypted files\
'*ega/mo_create_multiome_aspera_jobs_batch1.sh*' upload encryted data to EGA\


## scripts for other papers

### aging
'*aging_subset_atac_objects.R*' subset ATAC data object to samples to try in aging analysis
'*aging_merge_atac_objects.R*' merge ATAC data objects that were subsetted for samples to try in aging analysis
'*aging_subset_rna_object.R*' subset RNA data object to samples to try in aging analysis


### differential protein expression analysis
'*differential_protein/lc_covid_olink.R*'   perform differential protein analysis with a linear regression model on the long-covid status


### merge seurat objects
'*seurat_preprocess/lc_create_per_celltype_objects.R*'  create a Seurat object for each cell type, with only UT or both UT and 24hCA\
'*seurat_preprocess/lc_export_rna_objects.R*'  export Seurat objects for the long-covid lifelines data


### differentially accessible region detection
'*differential_accessibility/lc_differential_accessibility_limma_parameterised.R*'  check differential accessible regions using limma on the long-covid status


### differential gene expression analysis
'*differential_expression/lc_create_per_celltype_objects.R*'    create per-celltype Seurat objects to use for running limma dream on the long-covid status\
'*differential_expression/lc_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the long-covid status in the multiome gene expression data\
'*differential_expression/lc_create_limma_jobs.sh*'     create sbatch jobs for running limma per celltype on the long-covid status in the multiome gene expression data\
'*differential_expression/lc_filter_de.R*'     filter DE output for ones that we tried to replicate


### cell type composition GWAS
'*ctc_gwas/mo_ctc_gwas_wg2_step2.sh*'    for the cell type composition GWAS, perform step 2 of sc-eQTLgen WG2 to get consortium-compatible cell types\
'*ctc_gwas/mo_ctc_gwas_wg2_step3.sh*'    for the cell type composition GWAS, perform step 3 of sc-eQTLgen WG2 to get consortium-compatible cell types\
'*ctc_gwas/mo_ctc_gwas_wg2_step4.sh*'    for the cell type composition GWAS, perform step 4 of sc-eQTLgen WG2 to get consortium-compatible cell types


### cis- vs trans-acting QTLs
'*ctc_gwas/mo_annotate_qtls_with_caqtls.R*'    annotate QTLs on whether they overlap with caQTLs\
'*ctc_gwas/mo_check_scenic_vs_eqtlgen.R*'    check overlap between eQTLgen and SCENIC+\
'*ctc_gwas/mo_create_creqtl_confinement_cistrans_eqtlgen.R*'    create CRE-i-eQTL confinement using eQTLgen trans-genes


### fungus meta-analysis
'*fungus_metaanalysis/funmeta_lookup_funmeta_variants.R*'   look up fungal meta-analysis variants in our QTL data

## scripts no longer used
'*demultiplexing/mo_run_scrublet.py*'   python script to run Scrublet on a CellBender corrected 10x lane\
'*demultiplexing/mo_test_scrublet.ipynb*'   jupyter notebook to test Scrublet on individual 10x lanes\
'*differential_accessibility/mo_differential_accessibility_kimma_parameterised.R*'  check differential accessible regions using kimma on the stimulation status



## License

The code availabe in this repository is available under GPLv2 License:
https://www.gnu.org/licenses/old-licenses/gpl-2.0.html
