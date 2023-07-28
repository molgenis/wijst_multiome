#!/bin/bash

###################################################################
#Script Name	  : mo_create_multiome_csvs_batch1.sh
#Description	  : create csv files containing libraries needed for joint ATAC and RNA alignment
#Author       	: Roy Oelen
###################################################################

RNA_LANES_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/rna/sequence_data/'
ATAC_LANES_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/atac/sequence_data/'
SAMPLE_SHEET_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/mo_sample_sheet_batch1.tsv'
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

# check each run
for lane in ${LANES[*]}
  do
    # grab the data for this lane
    lane_data=$(grep ${lane} ${SAMPLE_SHEET_LOC})

    # the samples are the second entry
    samples=$(echo ${lane_data} | awk '{print $2}')

    # split the samples
    IFS=',' read -ra samples_split <<< "${samples}"

    # get the full RNA location
    full_rna_loc=${RNA_LANES_LOC}'/'${lane}'/'
    # and full ATAC location
    full_atac_loc=${ATAC_LANES_LOC}'/'${lane}'/'

    # create the full output location of the csv
    full_csv_loc=${CSVS_LOC}'/'${lane}'.csv'

    # write the header
    echo 'fastqs,sample,library_type' > ${full_csv_loc}

    # write each sample
    for sample in "${samples_split[@]}"
    do
        echo ${full_rna_loc}','${sample}',Gene Expression' >>  ${full_csv_loc}
        echo ${full_atac_loc}','${sample}',Chromatin Accessibility' >>  ${full_csv_loc}
    done
done