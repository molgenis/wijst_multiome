#!/bin/bash

###################################################################
#Script Name	  : mo_copy_trailing_cellranger_data.sh
#Description	  : copy any specific data from the cellranger output
#Args           : file or directory in the cellranger per-lane output that you want to copy
#Author       	: Roy Oelen
#example        : ./mo_copy_trailing_cellranger_data.sh \
#                 outs/per_barcode_metrics.csv \
###################################################################

# these we know
PRM_CELLRANGER_DIR='/groups/umcg-franke-scrna/prm02/projects/multiome/processed/joint/alignment/b38/'
TMP_CELLRANGER_DIR='/groups/umcg-franke-scrna/tmp04/projects/multiome/processed/joint/alignment/b38/'

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

# this is the relative path of what we want to copy from the cellranger lanes
DATA_TO_CP=$1

# check each lane
for lane in ${LANES[*]}
    do
    # paste together the full destination
    data_to_copy_to=${TMP_CELLRANGER_DIR}'/'${lane}'/'${DATA_TO_CP}
    # get the directory of this destination
    data_dir=${data_to_copy_to%/*}
    # make this directory
    mkdir -p ${data_dir}
    # get the full source path
    data_to_copy_from=${PRM_CELLRANGER_DIR}'/'${lane}'/'${DATA_TO_CP}
    # copy the data
    cp -r ${data_to_copy_from} ${data_dir}
done