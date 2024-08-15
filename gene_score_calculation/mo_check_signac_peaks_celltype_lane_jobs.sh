#!/bin/bash

###################################################################
#Script Name	  : mo_check_signac_peaks_celltype_lane_jobs.sh
#Description	  : 
#Args           :
#Author       	: Roy Oelen
###################################################################


# seurat related things
SEURAT_OBJECTS_LOC='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/'
SEURAT_PREPEND='mo_'
SEURAT_APPEND='multimodal_azi_mapped.rds'
# where to place the output
BASE_OUT_LOC='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/'
# what the column is of the cell types
CELL_TYPE_COLUMN='predicted.mo_10x_cell_type'
# if there is a classification
RECLASSIFICATION_TABLE='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cell_type_assignment/azimuth/celltypes_10x_ref_to_lowerres.tsv'
# where we will place the jobs
JOBS_LOC='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/jobs/'
# location of the R executable
R_EXEC='~/start_Rscript.sh'
# location of the R script
SCRIPT_LOC='/groups/umcg-franke-scrna/tmp03/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_check_signac_peaks_celltype_lane.R'

# settings for job
TIME='6:59:59'
NCPU='8'
MEMGB='64'

# we will run through all of the lanes
LANES=('230105_lane1' '230105_lane2' '230105_lane3' '230105_lane4' \
'230105_lane5' '230105_lane6' '230105_lane7' '230105_lane8' \
'230112_lane1' '230112_lane2' '230112_lane3' '230112_lane4' \
'230112_lane5' '230112_lane6' '230112_lane7' '230112_lane8' \
'230120_lane1' '230120_lane2' '230120_lane3' '230120_lane4' \
'230120_lane5' '230120_lane6' '230120_lane7' '230120_lane8' \
'230127_lane1' '230127_lane2' '230127_lane3' '230127_lane4' \
'230127_lane5' '230127_lane6' '230127_lane7' '230127_lane8' \
'230202_lane1' '230202_lane2' '230202_lane3' '230202_lane4' \
'230202_lane5' '230202_lane6' '230202_lane7' '230202_lane8' \
'230209_lane1' '230209_lane2' '230209_lane3' '230209_lane4' \
'230209_lane5' '230209_lane6' '230209_lane7' '230209_lane8' \
'230216_lane1' '230216_lane2' '230216_lane3' '230216_lane4' \
'230216_lane5' '230216_lane6' '230216_lane7' '230216_lane8' \
'230223_lane1' '230223_lane2' '230223_lane3' '230223_lane4' \
'230223_lane5' '230223_lane6' '230223_lane7' '230223_lane8' \
'230302_lane1' '230302_lane2' '230302_lane3' '230302_lane4' \
'230302_lane5' '230302_lane6' '230302_lane7' '230302_lane8' \
'230316_lane1' '230316_lane2' '230316_lane3' '230316_lane4' \
'230316_lane5' '230316_lane6' '230316_lane7' '230316_lane8' \
)

# we will only use these cell types
CELL_TYPES=( \
'B' 'CD4T' 'CD8T' 'DC' 'monocyte' 'NK' \
)

# check each lane
for lane in ${LANES[*]} ; do
    # check each cell type
    for cell_type in ${CELL_TYPES[*]} ; do
        # paste together the paths that we want
        seurat=${SEURAT_OBJECTS_LOC}'/'${SEURAT_PREPEND}${lane}${SEURAT_APPEND}
        out=${BASE_OUT_LOC}'/'${lane}'_'${cell_type}
        # as well as the actual job location
        job_name='calc_peak_'${lane}'_'${cell_type}
        job_loc_full=${JOBS_LOC}'/calc_peak_'${lane}'_'${cell_type}'_SBATCH.sh'
        log_out_full=${JOBS_LOC}'/calc_peak_'${lane}'_'${cell_type}'.out'
        log_err_full=${JOBS_LOC}'/calc_peak_'${lane}'_'${cell_type}'.err'
        # create the job
        echo -e "#!/usr/bin/env bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${log_out_full}
#SBATCH --error=${log_err_full}
#SBATCH --time=${TIME}
#SBATCH --cpus-per-task=${NCPU}
#SBATCH --mem=${MEMGB}gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE

${R_EXEC} ${SCRIPT_LOC} \\
    --seurat ${seurat} \\
    --out ${out} \\
    --cell_type_column ${CELL_TYPE_COLUMN} \\
    --cell_type ${cell_type} \\
    --reclassification_table ${RECLASSIFICATION_TABLE} \\

" > ${job_loc_full}

    done
done