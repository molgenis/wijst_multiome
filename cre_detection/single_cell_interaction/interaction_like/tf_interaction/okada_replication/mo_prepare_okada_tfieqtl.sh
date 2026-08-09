#!/bin/bash

###################################################################
#Script Name	  : mo_prepare_okada_tfieqtl.sh
#Description	  : take the Okada data and convert to format required for TFa-i-eQTL analysis
#Args           : 
#Author       	: Roy Oelen
###################################################################

# location of seurat file
OKADA_RDS_LOC='/groups/umcg-franke-scrna/tmp04/external_datasets/onek1k/seurat_objects/okada.rds'
# location of chunks
OKADA_CHUNK_LOC='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/okada/L1/all/'
# confinements
CONFINEMENT_LOC='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/featureVariantFile.w150k.filtered0.0001_cts.txt'
GENE_ANNO_LOC='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'

# convert object into chunks
#Rscript mo_create_hybrid_cre_inputs_rna.R \
#    --in ${OKADA_RDS_LOC} \
#    --out ${OKADA_CHUNK_LOC} \
#    --confinement ${CONFINEMENT_LOC} \
#    --gene_annotation_file ${GENE_ANNO_LOC} \
#    --donor_annotation_column Assignment

# location of the genotypes
GENOTYPE_DIR='/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_okada/genotype_input/'
# location of plink2 binary
PLINK_LOC='/groups/umcg-franke-scrna/tmp04/software/plink2_amd/plink2'

# convert genotype data from bgen to bed
${PLINK_LOC} \
  --bgen ${GENOTYPE_DIR}/EUR_imputed_hg38_varFiltered.bgen ref-last \
  --make-bed \
  --out ${GENOTYPE_DIR}/EUR_imputed_hg38_varFiltered

# split genotype data into chromosome-specific files
for chr in {1..22} X Y; do
   ${PLINK_LOC} --bfile ${GENOTYPE_DIR}/EUR_imputed_hg38_varFiltered --chr $chr --make-bed --out ${GENOTYPE_DIR}/EUR_imputed_hg38_varFiltered$chr
done