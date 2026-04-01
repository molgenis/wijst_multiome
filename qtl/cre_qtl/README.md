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
