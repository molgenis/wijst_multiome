#!/bin/bash

###################################################################
#Script Name	  : mo_create_limma_dar_topic_jobs.sh
#Description	  : create SBATCH jobs scripts to do DAR identification per cell type object
#Args           : location of seurat object, output location of DE, column in the metadata describing cell type, location to put job files
#Author       	: Roy Oelen
#example        : ./mo_create_limma_dar_topic_jobs.sh \
#/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_ \
#/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/topics20_aucell/pct01/ \
#cell_type \
#/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/jobs/topics20_aucell/pct01/ \
#/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/topic_annotations/mo_topic_20_aucell.tsv \
#0
###################################################################

# standard parameters
CORES='2'
MEMORY_GB='200'
TMP_SIZE='512mb'
RUNTIME='23:59:59'

# whether or not to permute
PERMUTE=0
# and if we permute, how many times
N_PERMUTATIONS=10
# if we are doing topics, how many
N_TOPICS=20

# location to the R script
script_loc='/groups/umcg-franke-scrna/tmp04/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_differential_accessibility_topics.R'

# the seurat object
seurat_objects_loc=$1 # like /groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240619_seuratv5_annotated_agesexcovid_
output_loc=$2 # like /groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/LONG_COVID/
cell_type_column=$3 # celltype_imputed_lowerres
jobs_loc=$4 # /groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/jobs/LONG_COVID/
topic_ann_loc=$5 # like /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/topic_annotations/mo_topic_20_aucell.tsv

# whether or not to permute
PERMUTE=$6
# and if we permute, how many times
N_PERMUTATIONS=$7

# set the number of iteractions as one if we are doing a true run
iterations=(1)
# however if we are permuting, we will do multiple iterations
if [ ${PERMUTE} -eq 1 ]
  then
    iterations=($(seq 1 1 ${N_PERMUTATIONS}))
  fi

# also we need multiple topics
topics=($(seq 1 1 ${N_TOPICS}))

# create the regex for the seurat objects
REGEX_SEURAT_OBJECTS=${seurat_objects_loc}'*.rds'

# check each RDS file
for object in ${REGEX_SEURAT_OBJECTS}
  do
  # create directories if they didn't exist
  mkdir -p ${output_loc}
  mkdir -p ${jobs_loc}
  # extract the file basename to use in the job name
  file_basename=$(basename "${object}" .rds)
  # do each iteration
  for i in "${iterations[@]}"
    do
        for t in "${topics[@]}"
            do
                # check if we are permuting
                if [ ${PERMUTE} -eq 0 ]
                    then
                        JOB_NAME='limma_topic'${t}${file_basename}
                    fi
                # we add which permutation it is, if we are permuting
                if [ ${PERMUTE} -eq 1 ]
                    then
                        JOB_NAME='limma_topic'${t}${file_basename}'_permutation_'${i}
                    fi
                JOB_LOC=${jobs_loc}'/'${JOB_NAME}'_SBATCH.sh'
                JOB_OUT=${jobs_loc}'/'${JOB_NAME}'.out'
                JOB_ERR=${jobs_loc}'/'${JOB_NAME}'.err'
                # echo the header
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
~/start_Rscript.sh '${script_loc}' \
-f '${object}' \
-o '${output_loc}' \
-c '${cell_type_column}' \
-m 10 \
-u 200 \
-l 2000 \' > ${JOB_LOC}
                # if we are permuting, add that as true
                if [ ${PERMUTE} -eq 1 ]
                    then
                        echo '--permute T \' >> ${JOB_LOC}
                    fi
                # and false if we are not
                if [ ${PERMUTE} -eq 0 ]
                    then
                        echo '--permute F \' >> ${JOB_LOC}
                    fi
                echo '--topic_ann '${topic_ann_loc}' \' >> ${JOB_LOC}
                echo '--topic Topic'${t}'
' >> ${JOB_LOC}
            done
    done
done


# /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/scripts/mo_create_limma_dar_topic_jobs.sh \
# /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_ \
# /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/output/topics20_aucell/pct01/ \
# cell_type \
# /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/limma_dream/jobs/topics20_aucell/pct01/ \
# /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/differential_accessibility/topic_annotations/mo_topic_20_aucell.tsv \
# 0