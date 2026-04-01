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
