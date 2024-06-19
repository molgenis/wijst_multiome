#!/bin/bash

###################################################################
#Script Name	  : lc_create_limma_jobs.sh
#Description	  : create SBATCH jobs scripts to do DE per cell type object
#Args           : location of seurat object, output location of DE, column in the metadata describing cell type, location to put job files
#Author       	: Roy Oelen
#example        : ./lc_create_limma_jobs.sh \
#/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240619_seuratv5_annotated_agesexcovid_ \
#/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/LONG_COVID/ \
#celltype_imputed_lowerres \
#/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/jobs/LONG_COVID/
###################################################################

# standard parameters
CORES='2'
MEMORY_GB='200'
TMP_SIZE='512mb'
RUNTIME='23:59:59'

# location to the R script
script_loc='/groups/umcg-franke-scrna/tmp03/users/umcg-roelen/singularity/rstudio-server/simulated_home/lc_differential_expression_limma_parameterised.R'

# the seurat object
seurat_objects_loc=$1 # like /groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240619_seuratv5_annotated_agesexcovid_
output_loc=$2 # like /groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/output/LONG_COVID/
cell_type_column=$3 # celltype_imputed_lowerres
jobs_loc=$4 # /groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/differential_expression/limma_dream/jobs/LONG_COVID/

# create the regex for the seurat objects
REGEX_SEURAT_OBJECTS=${seurat_objects_loc}'*.rds'

for object in ${REGEX_SEURAT_OBJECTS}
  do
  mkdir -p ${output_loc}
  mkdir -p ${jobs_loc}
  file_basename=$(basename "${object}" .rds)
    JOB_NAME='limma_'${file_basename}
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
-l 2000
' > ${JOB_LOC}
done