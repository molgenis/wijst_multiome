#!/bin/bash

###################################################################
#Script Name	  : mo_cellbend_rna.sh
#Description	  : create SBATCH jobs scripts to do ambient RNA correction
#Author       	: Roy Oelen, Martijn Vochteloo
###################################################################

LANE_DIR='/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/'
JOB_DIR='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/ambient_rna_correction/CellBender/jobs/joint/b38/include_introns/expectcells/'
OUTPUT_DIR='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/ambient_rna_correction/CellBender/output/joint/b38/include_introns/expectcells/'
OUTPUT_FILE_APPPEND='cellbent_feature_bc_matrix.h5'
INPUT_FILE_APPEND='outs/raw_feature_bc_matrix.h5'
CELL_NUMBERS_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_cell_numbers.tsv'

# option for setting expect_cells
SET_EXPECTED_CELLS=1

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
'230316_lane5' '230316_lane6' '230316_lane7' \
)


CORES='4'
MEMORY_GB='8'
TMP_SIZE='512mb'
RUNTIME='23:59:59'
NR_OF_GPUS='2'

# check each run
for lane in ${LANES[*]}
  do
    JOB_NAME='cellbend_'${lane}
    JOB_LOC=${JOB_DIR}'/'${JOB_NAME}'_SBATCH.sh'
    JOB_OUT=${JOB_DIR}'/'${JOB_NAME}'.out'
    JOB_ERR=${JOB_DIR}'/'${JOB_NAME}'.err'

    # output directory
    OUTPUT_LOC_FULL=${OUTPUT_DIR}'/'${lane}'/'

    # create that directory
    mkdir -p ${OUTPUT_LOC_FULL}

    # and the output file
    OUTPUT_FILE_FULL=${OUTPUT_LOC_FULL}'/'${OUTPUT_FILE_APPPEND}
    
    # the input file as well
    INPUT_FILE_FULL=${LANE_DIR}'/'${lane}'/'${INPUT_FILE_APPEND}

    # echo the header
    echo '#!/bin/bash
#SBATCH --job-name='${JOB_NAME}'
#SBATCH --output='${JOB_OUT}'
#SBATCH --error='${JOB_ERR}'
#SBATCH --time='${RUNTIME}'
#SBATCH --cpus-per-task='${CORES}'
#SBATCH --mem='${MEMORY_GB}'GB
#SBATCH --gres=gpu:a40:'${NR_OF_GPUS}'
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp='${TMP_SIZE}'
'> ${JOB_LOC}

    # load environment
    echo '~/miniconda3/bin/activate cellbender_env' >> ${JOB_LOC}

    # and CUDA
    echo 'ml CUDA/11.7.0' >> ${JOB_LOC}

    # build the job
    echo '~/miniconda3/envs/cellbender_env/bin/cellbender \
 remove-background \
  --input='${INPUT_FILE_FULL}' \
  --output='${OUTPUT_FILE_FULL}' \' >> ${JOB_LOC}

    # add the expected cells as a parameter if requested
    if [ ${SET_EXPECTED_CELLS} -eq 1 ]
    then
      # grab the data for this lane
      lane_data=$(grep ${lane} ${CELL_NUMBERS_LOC})

      # the cells are the second entry
      nr_cells=$(echo ${lane_data} | awk '{print $2}')

      # add that parameter to the bash script
      echo '  --expected-cells='${nr_cells}' \' >> ${JOB_LOC}
    fi

    echo '  --cuda' >> ${JOB_LOC}

done


# required environment set up like this:
# conda create -n cellbender_env python=3.10
# conda activate cellbender_env
# ml CUDA/11.7.0
# pip install tables torch torchvision
# git clone https://github.com/broadinstitute/CellBender
# pip install -e CellBender