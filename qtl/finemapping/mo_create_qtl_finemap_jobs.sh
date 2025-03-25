#!/usr/bin/env bash
###################################################################
#Script Name	  : mo_create_qtl_finemap_jobs.sh
#Description	  : create jobs to run finemapping on the QTL output files
#Args           : location of the QTL output folders, location to store the finemapping results, location to store the sbatch job files
#Author       	: Roy Oelen
#example        : 
# mo_create_qtl_finemap_jobs.sh \
#   /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/sc-eqtlgen/output/nominal_condition/L1/ \
#   /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_eqtl/sc-eqtlgen/output/nominal_condition/L1/ \
#   /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/finemapping/interaction_eqtl/sc-eqtlgen/jobs/nominal_condition/L1/ \
#   /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/genotype/imputed_hg38_all_anc
###################################################################


# location where the folders with QTL output are
QTL_DIR=$1
# location to place the finemapped output
FINEMAPPED_OUT_DIR=$2
# location to place the job files
FINEMAP_JOB_DIR=$3
# location of the genotype file used for the finemapping
GENOTYPE_LOC=$4

# these we always have
CORES='2'
MEMORY_GB='128'
TMP_SIZE='512mb'
RUNTIME='23:59:59'

# location to the R script
FM_SCRIPT_LOC='/groups/umcg-franke-scrna/tmp04/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_finemap_qtls.R'
FORMAT_SCRIPT_LOC='/groups/umcg-franke-scrna/tmp04/users/umcg-roelen/singularity/rstudio-server/simulated_home/lpmcv2_format_finemapping.R'
# and the R command
RSCRIPT_COMMAND='~/start_Rscript.sh'

# prepends and appends of results
FINEMAP_RESULT_PREPEND=''
FINEMAP_RESULT_APPEND='_finemapped.rds'
FINEMAP_FORMATTED_PREPEND=''
FINEMAP_FORMATTED_APPEND='_finemapped.tsv.gz'
JOB_PREPEND='finemap_'
JOB_APPEND=''
# and for the input
INPUT_FILE='iqtl_results_all.txt.gz'

# and how the columns are named
SIGNIFICANCE_COLUMN='empirical_feature_p_value'
SIGNIFICANCE_CUTOFF='0'
VARIANT_COLUMN='snp_id'
FEATURE_COLUMN='feature_id'
SLOPE_COLUMN='beta_SNP'
SE_COLUMN='beta_se_SNP'
N_SAMPLE='318'

# list all files that have the search
dirlist=(${QTL_DIR}*)
# loop the files and directories
for e in "${dirlist[@]}"; do
    # if is is a directory, make the script
    if [ -d "$e" ];
    then
        # get the basename of that folder
        e_basename=$(basename ${e})
        # get the full input path
        INPUT_PATH=${QTL_DIR}'/'${e_basename}'/'${INPUT_FILE}
        # construct the finemapped rds file
        FINEMAP_RESULT=${FINEMAPPED_OUT_DIR}'/'${FINEMAP_RESULT_PREPEND}${e_basename}${FINEMAP_RESULT_APPEND}
        # and table
        FINEMAP_TABLE=${FINEMAPPED_OUT_DIR}'/'${FINEMAP_FORMATTED_PREPEND}${e_basename}${FINEMAP_FORMATTED_APPEND}
        
        # construct the first command
        FM_COMMAND=${RSCRIPT_COMMAND}' '${FM_SCRIPT_LOC}' 
    --qtl_file '${INPUT_PATH}' 
    --genotype_file '${GENOTYPE_LOC}' 
    --significance_column '${SIGNIFICANCE_COLUMN}' 
    --significance_cutoff '${SIGNIFICANCE_CUTOFF}' 
    --output_rds '${FINEMAP_RESULT}' 
    --num_threads '${CORES}' 
    --variant_column '${VARIANT_COLUMN}' 
    --feature_column '${FEATURE_COLUMN}' 
    --slope_column '${SLOPE_COLUMN}' 
    --se_column '${SE_COLUMN}' 
    --n_sample '${N_SAMPLE}

        # and the second command
        FORMAT_COMMAND=${RSCRIPT_COMMAND}' '${FORMAT_SCRIPT_LOC}' --in '${FINEMAP_RESULT}' --out '${FINEMAP_TABLE}

        # now create a job name
        JOB_NAME=${JOB_PREPEND}${e_basename}${JOB_APPEND}
        JOB_LOC=${FINEMAP_JOB_DIR}'/'${JOB_NAME}'_SBATCH.sh'
        JOB_OUT=${FINEMAP_JOB_DIR}'/'${JOB_NAME}'.out'
        JOB_ERR=${FINEMAP_JOB_DIR}'/'${JOB_NAME}'.err'

        # echo the file together
        echo '#!/bin/bash
#SBATCH --job-name='${JOB_NAME}'
#SBATCH --output='${JOB_OUT}'
#SBATCH --error='${JOB_ERR}'
#SBATCH --time='${RUNTIME}'
#SBATCH --cpus-per-task='${CORES}'
#SBATCH --mem='${MEMORY_GB}'GB
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp='${TMP_SIZE}'

'${FM_COMMAND}'

'${FORMAT_COMMAND}'

' > ${JOB_LOC}

    fi
done