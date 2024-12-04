#!/bin/bash

###################################################################
#Script Name	  : mo_ctc_gwas_wg2_step3.sh
#Description	  : for the celltype composition gwas, do step3 of wg2
#Args           :
#Author       	: Roy Oelen
###################################################################

# go to the directory that the data is at
cd /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/
# create the directory
mkdir step3_hierscpred
# check each rds file
for i in $(ls step2_azimuth);
do
  out=$(echo $i | awk 'gsub(".rds", "")') # Use same base filename as output
  singularity run -B /groups/umcg-franke-scrna/tmp01/ /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/cell_classification.sif Rscript /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/Depracated/wg2-cell-type-classification/wg2_wijst2020/map_hierscpred.R --file /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/step2_azimuth/${i} --path /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/step3_hierscpred --out ${out} --batch lane
done