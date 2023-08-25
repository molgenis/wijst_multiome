#!/usr/bin/env bash

###################################################################
#Script Name	  : mo_create_demuxlet_jobs.sh
#Description	  : create SBATCH jobs scripts to do Demuxlet
#Author       	: Roy Oelen
###################################################################


#directory and file listings
LANES_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/"
LANE_BARCODE_APPEND="outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
BAM_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/popscle_tools_filtered/alignment/filtered_alignment/"
INDIVIDUAL_GENOTYPES="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/genotype_per_lane/"
OUTPUT_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/demuxlet/demuxlet_output/"
JOB_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/demuxlet/jobs/"
DEMUXLET_CMD="singularity exec --bind /groups/umcg-franke/tmp02/ /groups/umcg-franke-scrna/tmp02/software/sc-eqtlgen-consortium-pipeline/wg1/wg1-pipeline-20230308_2.simg popscle demuxlet"

#parameters used
TAG_GROUP="CB"
TAG_UMI="UB"
FIELD="GT"

for dir in "$LANES_DIR"/*lane*/ ; do

    lane_id=$(basename $dir)
    output_job=${JOB_DIR}"/demux_"${lane_id}"_SBATCH.sh"
    bam_loc=${BAM_DIR}/"mo_"${lane_id}"_genofiltered_gex_possorted_bam.bam"

    echo -e "#!/usr/bin/env bash
#SBATCH --job-name=demux_${lane_id}
#SBATCH --output=demux_${lane_id}.out
#SBATCH --error=demux_${lane_id}.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
	set -e
        ${DEMUXLET_CMD} \\
                --sam ${bam_loc} \\
                --tag-group ${TAG_GROUP} \\
                --tag-UMI ${TAG_UMI} \\
                --field ${FIELD} \\
                --vcf ${INDIVIDUAL_GENOTYPES}mo_bamsorted_${lane_id}_maf005.vcf.gz \\
                --out ${OUTPUT_DIR}${lane_id} \\
                --group-list ${LANES_DIR}${lane_id}/${LANE_BARCODE_APPEND}
        " > ${output_job}
done
