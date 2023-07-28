#!/bin/bash

###################################################################
#Script Name	  : mo_create_cellranger_arc_jobs_batch1.sh
#Description	  : create SBATCH jobs scripts to do alignment RNA and ATAC data
#Author       	: Roy Oelen
###################################################################

CELLRANGER_LOC='/groups/umcg-franke-scrna/tmp02/software/cellranger-arc-2.0.2/cellranger-arc'
REFDATA_LOC='/groups/umcg-franke-scrna/tmp02/external_datasets/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/'
SAMPLE_SHEET_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_sample_sheet_batch1.tsv'
OUTPUT_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/'
JOB_DIR='/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/jobs/b38/'
CSVS_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/combined_library_csvs/'

# these are the lanes to run through cellranger
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


CORES='22'
MEMORY_GB='128'
TMP_SIZE='512mb'
RUNTIME='23:59:59'

# check each run
for lane in ${LANES[*]}
  do
    JOB_NAME='align_'${lane}
    JOB_LOC=${JOB_DIR}'/'${JOB_NAME}'_SBATCH.sh'
    JOB_OUT=${JOB_DIR}'/'${JOB_NAME}'.out'
    JOB_ERR=${JOB_DIR}'/'${JOB_NAME}'.err'

    # temporary and permanent storage for cellranger output per lane
    OUTPUT_LOC_FULL=${OUTPUT_LOC}'/'${lane}'/'

    # grab the data for this lane
    lane_data=$(grep ${lane} ${SAMPLE_SHEET_LOC})

    # the samples are the second entry
    samples=$(echo ${lane_data} | awk '{print $2}')

    # and the location of the csvs
    library_csv=${CSVS_LOC}'/'${lane}'.csv'

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
'> ${JOB_LOC}

    # do the prerequisites
    echo 'mkdir -p '${OUTPUT_LOC_FULL}'/' >> ${JOB_LOC}

    # go to that directory
    echo 'cd '${OUTPUT_LOC}'' >> ${JOB_LOC}

    # for compatibility with older clusters, we need this
    echo 'export TENX_IGNORE_DEPRECATED_OS=1' >> ${JOB_LOC}

    # build the job
    echo ${CELLRANGER_LOC}' count \
--id='${lane}' \
--reference='${REFDATA_LOC}' \
--libraries='${library_csv}' \
--localcores='${CORES}' \
--localmem='${MEMORY_GB}'  
' >> ${JOB_LOC}

done