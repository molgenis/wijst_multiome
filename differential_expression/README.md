### differential gene expression analysis
'*differential_expression/1m_create_per_celltype_objects.R*'    create per-celltype Seurat objects from the NC2022 study which also used Candida stimulation\
'*differential_expression/1m_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the stimulation status in the NC2022 data\
'*differential_expression/mo_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the stimulation status in the multiome gene expression data\
'*differential_expression/mo_plot_de.R*'    plot differential gene expression results\
'*differential_expression/mo_create_limix_qtl_input_l2.R*'    create per-celltype Seurat objects at L2 resolution\
'*differential_expression/mo_merge_dar_de_tables.R.R*'    merge DAR and DE tables

### LONG COVID
'*differential_expression/lc_create_per_celltype_objects.R*'    create per-celltype Seurat objects to use for running limma dream on the long-covid status\
'*differential_expression/lc_differential_expression_limma_parameterised.R*'    perform differential gene expression using limma dream on the long-covid status in the multiome gene expression data\
'*differential_expression/lc_create_limma_jobs.sh*'     create sbatch jobs for running limma per celltype on the long-covid status in the multiome gene expression data\
'*differential_expression/lc_filter_de.R*'     filter DE output for ones that we tried to replicate
