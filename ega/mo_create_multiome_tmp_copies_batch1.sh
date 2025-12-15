#!/bin/bash

###################################################################
#Script Name	  : mo_create_multiome_tmp_copies_batch1.sh
#Description	  : create scripts that copy the relevant sequence files from prm to tmp with correct names
#Author       	: Roy Oelen
###################################################################

# the original locations of the sequence files
RNA_LANES_LOC='/groups/umcg-franke-scrna/prm02/projects/multiome/raw/rna/sequence_data/'
ATAC_LANES_LOC='/groups/umcg-franke-scrna/prm02/projects/multiome/raw/atac/sequence_data/'
# where to place the copied files
RNA_LANES_COPIED_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/rna/sequence_data_unencrypted/'
ATAC_LANES_COPIED_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/atac/sequence_data_unencrypted/'
# to be able to distinguish the ATAC and RNA files, we have to add a prepend to them
RNA_LANES_PREPEND='RNA_'
ATAC_LANES_PREPEND='ATAC'
# add we'll make a script per lane so we can do the copies in batches
COPY_SCRIPTS_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/scripts/copy_scripts/'
# these are the lanes with sequence data
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

# check each lane
for lane in ${LANES[*]}
    do
    # get the full RNA location
    full_rna_loc=${RNA_LANES_LOC}'/'${lane}'/'
    # and full ATAC location
    full_atac_loc=${ATAC_LANES_LOC}'/'${lane}'/'
    # get the full copied RNA location
    full_rna_copied_loc=${RNA_LANES_COPIED_LOC}'/'${lane}'/'
    # and full copied ATAC location
    full_atac_copied_loc=${ATAC_LANES_COPIED_LOC}'/'${lane}'/'
    # make these directories
    mkdir -p ${full_rna_copied_loc}
    mkdir -p ${full_atac_copied_loc}
    # create the copy script
    copy_script_loc=${COPY_SCRIPTS_LOC}'/mo_copy_'${lane}'.sh'
    # write the header for the script
    echo '#!/bin/bash
' > ${copy_script_loc}
    # create the regex for the files
    REGEX_RNA=${full_rna_loc}'*'
    # check each file
    for object in ${REGEX_RNA}
        do
        # get the basename
        object_basename=$(basename "${object}")
        # make the new full path
        object_new_loc=${full_rna_copied_loc}${RNA_LANES_PREPEND}${object_basename}
        # paste the copy operation
        copy_operation_rna='cp '${object}' '${object_new_loc}
        # and echo this into the script
        echo ${copy_operation_rna} >> ${copy_script_loc}
    done
    # echo a newline
    echo '
' >> ${copy_script_loc}
    # create the regex for the ATAC files
    REGEX_ATAC=${full_rna_loc}'*'
    # check each file
    for object in ${REGEX_ATAC}
        do
        # get the basename
        object_basename=$(basename "${object}")
        # make the new full path
        object_new_loc=${full_atac_copied_loc}${ATAC_LANES_PREPEND}${object_basename}
        # paste the copy operation
        copy_operation_atac='cp '${object}' '${object_new_loc}
        # and echo this into the script
        echo ${copy_operation_atac} >> ${copy_script_loc}
    done
done