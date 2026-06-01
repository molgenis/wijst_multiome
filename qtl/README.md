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
