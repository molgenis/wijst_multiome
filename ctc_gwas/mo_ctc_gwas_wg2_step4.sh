#!/bin/bash

###################################################################
#Script Name	  : mo_ctc_gwas_wg2_step4.sh
#Description	  : for the celltype composition gwas, do step4 of wg2
#Args           :
#Author       	: Roy Oelen
###################################################################

singularity exec -B /groups/umcg-franke-scrna/tmp01/ /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/cell_classification.sif Rscript /reduce.R --file /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/step3_hierscpred --out /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/reduced_data --path /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/ctc_gwas/wg2/step4_reduce
