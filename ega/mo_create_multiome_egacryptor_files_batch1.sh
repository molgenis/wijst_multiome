#!/bin/bash

###################################################################
#Script Name	  : mo_create_multiome_egacryptor_files_batch1.sh
#Description	  : create jobs that encrypt the sequence files for EGA
#Author       	: Roy Oelen
###################################################################


# standard parameters
CORES='4'
MEMORY_GB='64'
TMP_SIZE='512mb'
RUNTIME='23:59:59'

# where the unencrypted files are
RNA_LANES_UNENCRYPTED_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/rna/sequence_data_unencrypted/'
ATAC_LANES_UNENCRYPTED_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/atac/sequence_data_unencrypted/'
# where the encrypted files will be placed
RNA_LANES_ENCRYPTED_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/rna/sequence_data_encrypted/'
ATAC_LANES_ENCRYPTED_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/atac/sequence_data_encrypted/'
# add we'll make a job per lane so we can do the encryption in batches
JOB_SCRIPTS_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/raw/ega/jobs/'
# we need EGA cryptor
EGACRYPTOR_LOC='/groups/umcg-franke-scrna/tmp02/software/egacryptor/EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar'

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
    # set the paths for the input and output
    full_rna_loc=${RNA_LANES_UNENCRYPTED_LOC}'/'${lane}'/'
    # and full ATAC location
    full_atac_loc=${ATAC_LANES_UNENCRYPTED_LOC}'/'${lane}'/'
    # get the full copied RNA location
    full_rna_encrypted_loc=${RNA_LANES_ENCRYPTED_LOC}'/'${lane}'/'
    # and full copied ATAC location
    full_atac_encrypted_loc=${ATAC_LANES_ENCRYPTED_LOC}'/'${lane}'/'
    
    # generate the commands
    rna_command='java -jar '${EGACRYPTOR_LOC}' -i '${full_rna_loc}' -o '${full_rna_encrypted_loc}
    atac_command='java -jar '${EGACRYPTOR_LOC}' -i '${full_atac_loc}' -o '${full_atac_encrypted_loc}
    
    # name the job
    JOB_NAME='mo_egacrypt_'${lane}
    # set paths to job and out/err
    JOB_LOC=${JOB_SCRIPTS_LOC}'/'${JOB_NAME}'_SBATCH.sh'
    JOB_OUT=${JOB_SCRIPTS_LOC}'/'${JOB_NAME}'.out'
    JOB_ERR=${JOB_SCRIPTS_LOC}'/'${JOB_NAME}'.err'

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
    # load Java in the job
    echo 'ml Java
' >> ${JOB_LOC}
    # send the encryption jobs
    echo ${rna_command} >> ${JOB_LOC}
    echo ${atac_command} >> ${JOB_LOC}
done