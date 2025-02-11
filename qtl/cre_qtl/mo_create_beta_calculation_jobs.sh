#!/bin/bash

###################################################################
#Script Name	  : mo_create_beta_calculation_jobs.sh
#Description	  : create beta and se file using multiomics data
#Author       	: Roy Oelen
#Example
# ./mo_create_beta_calculation_jobs.sh \
#  /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/pseudobulk_replication/matrices/monocyte/ \
#  /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/pseudobulk_replication/betas_ps/monocyte/ \
#  /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/coeqtl/trial_run/cre_lists/mono_cre_scenic_and_pseudo.tsv.gz \
#  /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/cre_eqtl/pseudobulk_replication/jobs/monocyte/ \
#  RNA \
#  peaks
###################################################################

# name parameters
PER_SAMPLE_MATRICES=$1
PER_SAMPLE_OUTPUT=$2
CRE_LOC=$3
JOB_OUT_LOC=$4
EXPRESSION_ASSAY=$5
ACCESSIBILITY_ASSAY=$6

# some defaults
RUNTIME='05:59:59'
CORES='2'
MEMORY_GB='16'
TMP_SIZE='512MB'
SCRIPT_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/qtl/coeqtl/trial_run/scripts/mo_calculate_atac_rna_betas.py'
PYTHON_BIN='~/miniconda3/envs/gpu_env/bin/python'

# list all files in directory
dirlist=(${PER_SAMPLE_MATRICES}*)

# check each run
for dir_full in ${dirlist[*]}
    do
    # check if it is actually a directory
    if [ -d "$dir_full" ];
        then
        # extract the name of the directory
        dir="$(basename ${dir_full})"
        # paste together the paths
        expression_folder_loc=${PER_SAMPLE_MATRICES}'/'${dir}/${EXPRESSION_ASSAY}'/'
        accessibility_folder_loc=${PER_SAMPLE_MATRICES}'/'${dir}/${ACCESSIBILITY_ASSAY}'/'
        # check if both of the folders exist
        if [ -d "$expression_folder_loc" ];
            then
            if [ -d "$accessibility_folder_loc" ];
                then
                # paste together the output location
                output_loc_full=${PER_SAMPLE_OUTPUT}'/'${dir}'/'
                # get the name for the job
                job_name='betacalc_'${dir}
                # create the output for the job
                output_job_full=${JOB_OUT_LOC}'/'${job_name}'_SBATCH.sh'
                # and out/err
                output_err_full=${JOB_OUT_LOC}'/'${job_name}'.err'
                output_out_full=${JOB_OUT_LOC}'/'${job_name}'.out'
                # make the actual file
                echo '#!/bin/bash
#SBATCH --job-name='${job_name}'
#SBATCH --output='${output_out_full}'
#SBATCH --error='${output_err_full}'
#SBATCH --time='${RUNTIME}'
#SBATCH --cpus-per-task='${CORES}'
#SBATCH --gres=gpu:a40:1
#SBATCH --mem='${MEMORY_GB}'GB
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp='${TMP_SIZE}'
'> ${output_job_full}
                # do the prerequisites
                echo 'mkdir -p '${output_loc_full}'/' >> ${output_job_full}
                echo 'conda activate gpu_env' >> ${output_job_full}
                # also add the actual work
                echo ${PYTHON_BIN}' '${SCRIPT_LOC}' \
    --expression_folder '${expression_folder_loc}' \
    --chromatin_folder '${accessibility_folder_loc}' \
    --output_folder '${output_loc_full}' \
    --use_gpu \
    --cre_loc '${CRE_LOC}'
'>> ${output_job_full}
            fi
        fi
    fi
done