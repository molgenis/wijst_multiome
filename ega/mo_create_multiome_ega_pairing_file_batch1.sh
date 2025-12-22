#!/bin/bash

############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_multiome_ega_pairing_file_batch1.sh
# Function: create pairing file of checksums for EGA
# Arguments: location of pairing file to create
# Example: ./mo_create_multiome_ega_pairing_file_batch1.sh /groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/mo_pairing_file_lanes.csv
############################################################################################################################


# location of the encrypted md5s
RNA_ENCRYPTED='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/rna/sequence_data_encrypted/'
ATAC_ENCRYPTED='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/atac/sequence_data_encrypted/'
# location of the sample mapping file
#PAIRING_FILE_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/mo_pairing_file_lanes.csv'
PAIRING_FILE_LOC=$1

# the lanes we want to do
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

# define function for writing md5s to file so we can use if for RNA and ATAC
add_to_pairing_file() {
  # get sequence data parameter
  SEQUENCE_DATA_LOC=$1
  # check each lane
  for lane in ${LANES[*]}
    do
    # paste together the location of the sequence file
    SEQUENCE_R1_PATTERN=${SEQUENCE_DATA_LOC}${lane}'/*R1*.fastq.gz.md5'
    # match
    R1_FILES_MATCH=( ${SEQUENCE_R1_PATTERN} )
    # check if the glob expanded or not, if there is no expansion, there is no match
    if [[ ${R1_FILES_MATCH[0]} == "$SEQUENCE_R1_PATTERN" ]]; then
        echo 'skipping '${lane}' due to no file match'
    else
        # now check each file
        for R1 in ${R1_FILES_MATCH[*]}
            do
            # get the corresponding R2
            R2=${R1/R1/R2}
            # get the basenames of the R1 and R2
            R1_NOPATH=$(basename ${R1} '.md5')
            R2_NOPATH=$(basename ${R2} '.md5')
            # get the location of the corresponding md5s
            R1_MD5_UNENCRYPTED_LOC=${SEQUENCE_DATA_LOC}${lane}/${R1_NOPATH}'.md5'
            R2_MD5_UNENCRYPTED_LOC=${SEQUENCE_DATA_LOC}${lane}/${R2_NOPATH}'.md5'
            R1_MD5_ENCRYPTED_LOC=${SEQUENCE_DATA_LOC}${lane}/${R1_NOPATH}'.gpg.md5'
            R2_MD5_ENCRYPTED_LOC=${SEQUENCE_DATA_LOC}${lane}/${R2_NOPATH}'.gpg.md5'
            # get the md5s from the files
            R1_UNENCRYPTED_MD5=$(sed -n '1p' ${R1_MD5_UNENCRYPTED_LOC} | awk '{print $1}' | tr -d '\n') # this file is a bit messy, we need the first row and column value [0,0] and then need to remove the newline
            R2_UNENCRYPTED_MD5=$(sed -n '1p' ${R2_MD5_UNENCRYPTED_LOC} | awk '{print $1}' | tr -d '\n')
            R1_ENCRYPTED_MD5=$(cat ${R1_MD5_ENCRYPTED_LOC})
            R2_ENCRYPTED_MD5=$(cat ${R2_MD5_ENCRYPTED_LOC})
            # now strip extention and read from the filenames so we get just the sample names, then remove the newline
            SAMPLE=$(echo "${R1_NOPATH}" | cut -d'_' -f3)
            # start writing to file
            # Sample Alias,First Fastq File,First Checksum,First Unencrypted checksum,Second Fastq File,Second Checksum,Second Unencrypted checksum
            LINE_TO_WRITE=${SAMPLE}','${R1_NOPATH}','${R1_ENCRYPTED_MD5}','${R1_UNENCRYPTED_MD5}','${R2_NOPATH}','${R2_ENCRYPTED_MD5}','${R2_UNENCRYPTED_MD5}
            echo ${LINE_TO_WRITE} >> ${PAIRING_FILE_LOC}
        done
    fi
  done
}

# write the header
echo 'Sample Alias,First Fastq File,First Checksum,First Unencrypted checksum,Second Fastq File,Second Checksum,Second Unencrypted checksum' > ${PAIRING_FILE_LOC}
# add to the file the RNA data
add_to_pairing_file ${RNA_ENCRYPTED}
# and the ATAC data
add_to_pairing_file ${ATAC_ENCRYPTED}
