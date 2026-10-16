#!/bin/bash

###################################################################
#Script Name	  : mo_create_maf_files.sh
#Description	  : take the Okada genotype data and create MAF files
#Args           : 
#Author       	: Roy Oelen
###################################################################

# location of the genotypes
GENOTYPE_DIR='./input/genotype/'
# location of plink2 binary
PLINK_LOC='/groups/umcg-franke-scrna/tmp04/software/plink2_amd/plink2'

# split genotype data into chromosome-specific files
for chr in {1..22}; do
   ${PLINK_LOC} --bfile ${GENOTYPE_DIR}/AES_imputed_hg38_varFiltered_chr${chr} --freq cols=chrom,pos,ref,alt,altfreq,nobs --out ${GENOTYPE_DIR}/AES_imputed_hg38_varFiltered_chr${chr}
   sha256sum ${GENOTYPE_DIR}/AES_imputed_hg38_varFiltered_chr${chr}.afreq > ${GENOTYPE_DIR}/AES_imputed_hg38_varFiltered_chr${chr}.afreq.sha256
done
