#!/bin/bash

###################################################################
#Script Name	  : mo_prepare_nc2022_genotypes.sh
#Description	  : prepare the genotype files for the interaction-eQTL replication
#Args           : 
#Author       	: Roy Oelen
#
###################################################################

# convert each chromosome to bgen format
CHROMS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22')
for chrom in ${CHROMS[*]}
    do
    /groups/umcg-franke-scrna/tmp02/software/plink2_amd/plink2 \
        --vcf /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/imputed/vcf_all_merged/imputed_hg38_info_filled.vcf.gz \
        --export bgen-1.2 \
        --out /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/imputed/bgen_split/imputed_hg38_info_filled_chr${chrom} \
        --chr ${chrom}
done

# also make a bed file for KING
/groups/umcg-franke-scrna/tmp02/software/plink2_amd/plink2 \
    --vcf /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/imputed/vcf_all_merged/imputed_hg38_info_filled.vcf.gz \
    --make-bed \
    --out  /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/imputed/plink_converted/imputed_hg38_info_filled

# and calculate the kinship
/groups/umcg-franke-scrna/tmp02/software/king-2.3.2/king \
    -b /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/imputed/plink_converted/imputed_hg38_info_filled.bed \
    --kinship \
    --prefix /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/imputed/king_kinship/imputed_hg38_info_filled_king_kinship

# index each genotype file by opening it once
python /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_wijst2018/coeqtl_redo_test/software/snakemake-to-map-coeqtls/scripts/coeqtl_open_bgens.py \
    --path /groups/umcg-franke-scrna/tmp02/projects/sc-eqtlgen-consortium-pipeline/wg1-preprocessing/wg1_wijst2020/genotype/imputed/bgen_split/ \
    --regex imputed_hg38_info_filled_chr*.bgen
