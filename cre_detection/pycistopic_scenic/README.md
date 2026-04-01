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
