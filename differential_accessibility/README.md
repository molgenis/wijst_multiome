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

### LONG-COVID
'*differential_accessibility/lc_differential_accessibility_limma_parameterised.R*'  check differential accessible regions using limma on the long-covid status
