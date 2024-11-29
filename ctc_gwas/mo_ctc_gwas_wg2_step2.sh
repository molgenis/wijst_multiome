#!/bin/bash

###################################################################
#Script Name	  : mo_ctc_gwas_wg2_step2.sh
#Description	  : for the celltype composition gwas, do step2 of wg2
#Args           :
#Author       	: Roy Oelen
###################################################################

# go to the directory that the data is at
cd /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/
# create the directory
mkdir step2_azimuth
# check each rds file
for i in $(ls step1_split);
do
  out=$(echo $i | awk 'gsub(".rds", "")') # Use same base filename as output
  singularity run -B /groups/umcg-franke-scrna/tmp01/ /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/cell_classification.sif Rscript /map_azimuth.R --file /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/step1_split/${i} --path /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/step2_azimuth --out ${out} --batch lane
done