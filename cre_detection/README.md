### CRE detection
'*cre_detection/mo_scenicplus_env.yml*' environment with packages used for CRE detection


#### comparisons
'*cre_detection/comparisons/mo_check_opposite_cre_effects.R'* check CRE-gene pairs that have opposite effects in different methods\
'*cre_detection/comparisons/mo_cre_method_outputs_comparison.R*'  compare SCENIC to pseudobulk and regression-based models\
'*cre_detection/comparisons/mo_create_creqtl_plots.R*'  plot CRE-i-eQTL effects\
'*cre_detection/comparisons/mo_do_full_cre_overlap_check.R*'  check overlap and concordance of TF-CRE-gene sets across CRE detection methods\
'*cre_detection/comparisons/mo_get_encode_cre_genes_to_cpeaks.R*'  add encode region-gene link information to cpeaks\
'*cre_detection/comparisons/mo_merge_cpeaks_with_screenv4.R*'  add encode screen v4 information to cpeaks defined regions\
'*cre_detection/comparisons/mo_overlap_qtl_with_cres.R*'  add information on cpeaks regions to where variants might be located in


##### single-cell CRE interaction analysis
'*cre_detection//interaction_analysis/mo_hybrid_cre_interaction_overlaps.R*'  plot interaction-eQTL at single-cell level with TF or ATAC as interaction terms overlaps across methods
'*cre_detection//single_cell_interaction/mo_hybrid_cre_plot_interaction.R*'  plot specific interaction-eQTLs from chunk\
'*cre_detection//single_cell_interaction/interaction_like/mo_hybrid_cre_interaction.R*'  perform single-cell interaction-eQTL analysis\
'*cre_detection//single_cell_interaction/interaction_like/mo_hybrid_interactions.smk*'  merge single-cell interaction-eQTL analysis chunks\
'*cre_detection//single_cell_interaction/interaction_like/mo_hybrid_cre_interaction_merge.R*'  snakemake to run single-cell interaction-eQTL analysis chunks\
'*cre_detection//single_cell_interaction/interaction_like/mo_plot_qtl_method_overlap.Rmd*'  overlap interaction-eQTLs with eQTL/caQTL/SCENIC+\
'*cre_detection//single_cell_interaction/interaction_like/region_interaction/mo_cre_sccre_interaction_confinement_r2g.R*'  create CRE-i-eQTL confinement file\
'*cre_detection//single_cell_interaction/interaction_like/region_interaction/mo_hybrid_interactions_regions_template.yaml*'  yaml configuration for CRE-i-eQTL run\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/mo_add_pseudobulked_tf_activities.R*'  calculate pseudobulked TF activities\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/mo_calc_auc_eregulons.R*'  calculate gene AUC based TF activities\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/mo_compare_tf_to_region_ieqtls.Rmd*'  compare TF-i-eQTLs to CRE-i-eQTLs\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/mo_cre_sccre_interaction_confinement.R*'  create TF-i-eQTLs confinement\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/mo_hybrid_interations_template.yaml*'  yaml configuration for TF-i-eQTLs run\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/mo_recalc_tf_i_eqtl_auc_eregulons.R*'  recalculate TF-i-eQTL TF activities excluding TF-i-eGenes\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/mo_hybrid_interations_template_replication.yaml*'  yaml configuration for TF-i-eQTLs validation run that removed TF-i-eGenes from TF activities\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/viacheslav_replication/mo_hybrid_interations_template_viacheslavrep.yaml*'  yaml configuration for TF-i-eQTLs replication from Viacheslav paper\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/coeqtl_replication/mo_sccre_coeqtl_tf_interaction_confinement.R*'  create confinement for replicating co-eQTL TF-gene pairs\
'*cre_detection//single_cell_interaction/interaction_like/tf_interaction/coeqtl_replication/mo_hybrid_interations_template_coeqtlrep.yaml*'  yaml configuration for TF-i-eQTLs replication from coeqtl paper\
'*cre_detection//single_cell_interaction/coeqtl_like/mo_hybrid_cre_coeqtl_replication.R*'  perform co-eQTL style TF-i-eQTL analysis for chunk\
'*cre_detection//single_cell_interaction/coeqtl_like/mo_hybrid_interactions_coeqtl.smk*'  perform co-eQTL style TF-i-eQTL analysis across chunks\
'*cre_detection//single_cell_interaction/coeqtl_like/mo_hybrid_interations_coeqtl_template.yaml*'  configuration to perform co-eQTL style TF-i-eQTL analysis across chunks

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

