#!/bin/bash

###################################################################
#Script Name	  : mo_do_gwas_colocs.sh
#Description	  : do GWAS coloc with SUSie
#Args           : name of the dataset, prepend of gwas files, append of gwas files
#Author       	: Roy Oelen
#example        : ./mo_do_gwas_colocs.sh \
# White_blood_cell_count \
# White_blood_cell_count__ \
# ___gwas.txt.gz
###################################################################

# the parameters
# DS2_NAME='White_blood_cell_count'
# DS2_PREPEND='White_blood_cell_count__'
# DS2_APPEND='___gwas.txt.gz'
DS2_NAME=$1
DS2_PREPEND=$2
DS2_APPEND=$3
# the cell types
CTS=('B' 'CD4T' 'CD8T' 'DC' 'monocyte' 'NK')

for CT in ${CTS}
  do
    ~/start_Rscript.sh ~/mo_coloc_traits_eqtlgen.R \
    --dataset1_in /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_eqtl/sc-eqtlgen/output/ut_and_24hca_significant/L1/${CT}_finemapped.tsv.gz \
    --dataset2_in_directory /groups/umcg-franke-scrna/tmp04/external_datasets/GWAS/eqtlgen_phase2_processed/${DS2_NAME}/ \
    --dataset2_in_prepend ${DS2_PREPEND} \
    --dataset2_in_append ${DS2_APPEND} \
    --dataset1_name ${CT} \
    --dataset2_name ${DS2_NAME} \
    --output_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/colocalization/eqtl_gwas/eqtlgen_processed/${DS2_NAME}/${CT}.tsv.gz \
    --binary_rds_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/colocalization/eqtl_gwas/eqtlgen_processed/${DS2_NAME}/${CT}.rds \
    --variant_mapping_loc /groups/umcg-franke-scrna/tmp04/external_datasets/GWAS/eqtlgen_phase2_processed/1000G-30x_index.parquet
done
