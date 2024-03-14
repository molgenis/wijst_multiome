#!/bin/bash

###################################################################
#Script Name	  : mo_create_cpeaks_jobs.sh
#Description	  : create SBATCH jobs scripts to do peak calling with the cpeaks reference
#Author       	: Roy Oelen
###################################################################

#LANE_DIR='/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/'
LANE_DIR='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/rounded_fragments/'
OUTPUT_DIR='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/output/rounded/'
JOB_DIR='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/jobs/rounded/'
#FRAGMENTS_APPEND='outs/atac_fragments.tsv.gz'
FRAGMENTS_APPEND='_rounded_fragments.tsv.gz'
CPEAKS_DIR='/groups/umcg-franke-scrna/tmp03/software/cPeaks/'
BARCODES_DIR='/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/credible_barcodes/'
BARCODES_PREPEND=''
BARCODES_APPEND='.txt'

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
'230316_lane5' '230316_lane6' '230316_lane7' '230316_lane8'\
)


CORES='8'
MEMORY_GB='86'
TMP_SIZE='512mb'
RUNTIME='23:59:59'
NR_OF_GPUS='1'

# check each run
for lane in ${LANES[*]}
  do
    JOB_NAME='cpeak_'${lane}
    JOB_LOC=${JOB_DIR}'/'${JOB_NAME}'_SBATCH.sh'
    JOB_OUT=${JOB_DIR}'/'${JOB_NAME}'.out'
    JOB_ERR=${JOB_DIR}'/'${JOB_NAME}'.err'

    # output directory
    output_loc=${OUTPUT_DIR}'/'${lane}'/'

    # location of fragments
    #fragment_loc=${LANE_DIR}'/'${lane}'/'${FRAGMENTS_APPEND}
    fragment_loc=${LANE_DIR}'/'${lane}''${FRAGMENTS_APPEND}

    # location of barcodes
    barcodes_loc=${BARCODES_DIR}'/'${BARCODES_PREPEND}${lane}${BARCODES_APPEND}

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
' > ${JOB_LOC}

# load environment
echo '~/miniconda3/bin/activate cpeaks_env' >> ${JOB_LOC}

# go to the directory that has the cpeaks script
echo 'cd '${CPEAKS_DIR} >> ${JOB_LOC}

# execute script
echo '~/miniconda3/envs/cpeaks_env/bin/python '${CPEAKS_DIR}'main.py \
    --fragment_path '${fragment_loc}' \
    --barcode_path '${barcodes_loc}' \
    --output '${output_loc}' \
    --num_cores '${CORES}'
' >> ${JOB_LOC}

done


# conda environment was created with
# conda create -n cpeaks_env python==3.9
# conda activate cpeaks_env
# pip install anndata scanpy numpy tqdm joblib