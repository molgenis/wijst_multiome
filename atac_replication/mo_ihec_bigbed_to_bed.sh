#!/bin/bash

###################################################################
#Script Name	  : mo_ihec_bigbed_to_bed.sh
#Description	  : convert all the bigbed files to bed files
#Author       	: Roy Oelen
###################################################################

# input directory
IHEC_IN_DIR='/groups/umcg-franke-scrna/tmp02/external_datasets/ihec/chromatin/'
# output directory
IHEC_OUT_DIR='/groups/umcg-franke-scrna/tmp02/external_datasets/ihec/chromatin/'
# bigbed file extention
BIGBED_APPEND='.bigBed'
# bed extention
BED_APPEND='.bed'

# location of the executable for bigbed to bed
BIGBED_TO_BED_BIN_LOC='/groups/umcg-franke-scrna/tmp02/software/ucsc/bigBedToBed'

# regex we'll use for the search variable
FILEREGEX='.+'${BIGBED_APPEND}

# list the files and dirs
dirlist=(${IHEC_IN_DIR}*)

# loop the files and directories
for e in "${dirlist[@]}"
    do
    # if it is a file
    if [ -f "$e" ];
        then
        # now check if it matches with our search
        if [[ ${e##*/} =~ ${FILEREGEX} ]];
            then
            # get the filename without the path
            base_e="$(basename -- $e)"
            # do a search and replace
            e_replaced=$(echo $base_e | sed "s/$BIGBED_APPEND/$BED_APPEND/g")
            # actually do the conversion
            ${BIGBED_TO_BED_BIN_LOC} \
                ${e} \
                ${IHEC_OUT_DIR}'/'${e_replaced}
        fi
    fi
done