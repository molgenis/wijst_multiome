#!/bin/bash

###################################################################
#Script Name	  : mo_create_souping_samples_jobs.sh
#Description	  : create sbatch jobs for Souporcell
#Args           :
#Author       	: Roy Oelen
###################################################################


#directory and file listings
LANES_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/"
LANE_READGROUP_APPEND="outs/gex_possorted_bam.bam"
LANE_BARCODE_APPEND="outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
OUTPUT_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/souporcell_output/"
JOB_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/jobs/"
SOUPOR_IMAGE="/groups/umcg-franke-scrna/tmp02/software/sc-eqtlgen-consortium-pipeline/wg1/WG1-pipeline-QC_wgpipeline.simg"
GENOME_LOC="/groups/umcg-franke-scrna/tmp02/external_datasets/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa"
#SAMPLE_SHEET_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/lpmcv2_sample_sheet_full.tsv'
COMMON_VARIANTS_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/imputed_hg38_all_anc_mmaf005_chrprepend.vcf.gz'

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


for lane in ${LANES[*]} ; do

    lane_id=${lane}
    output_folder=${OUTPUT_DIR}/${lane_id}/
    mkdir -p ${output_folder}
    output_job=${JOB_DIR}/soup_${lane_id}_SBATCH.sh

    # grab the data for this lane
#    lane_data=$(grep ${lane} ${SAMPLE_SHEET_LOC})
#
#    # the samples are the second entry
#    samples=$(echo ${lane_data} | awk '{print $3}')
#    # split samples by comma
#    IFS=',' read -r -a array <<< "${samples}"
#    echo ${#array[@]}
#    # get the length of the array
#    nr_of_samples=${#array[@]}
    nr_of_samples=8

    echo -e "#!/usr/bin/env bash
#SBATCH --job-name=soup_${lane_id}
#SBATCH --output=soup_${lane_id}.out
#SBATCH --error=soup_${lane_id}.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
        set -e
        cd ${output_folder}
        gunzip -c ${LANES_DIR}/${lane_id}/${LANE_BARCODE_APPEND} > ${output_folder}/barcodes.tsv
        export SINGULARITY_BINDPATH=\"/groups/umcg-franke-scrna/tmp02/projects/multiome/,/groups/umcg-franke-scrna/tmp02/external_datasets/,/groups/umcg-franke-scrna/tmp02/software/ \"
        singularity exec ${SOUPOR_IMAGE} souporcell_pipeline.py \
-i ${LANES_DIR}/${lane_id}/${LANE_READGROUP_APPEND} \
-b ${output_folder}/barcodes.tsv \
-f ${GENOME_LOC} \
-t 8 \
-o ${output_folder}/ \
-k ${nr_of_samples} \
--skip_remap True \
--common_variants ${COMMON_VARIANTS_LOC}
        " > ${output_job}
done
