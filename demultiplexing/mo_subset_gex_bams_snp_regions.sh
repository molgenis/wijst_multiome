#!/usr/bin/env bash

###################################################################
#Script Name	  : mo_subset_gex_bams_snp_regions.sh
#Description	  : create SBATCH jobs scripts to filter the RNA alignment data for regions that overlap SNPs
#Author       	: Roy Oelen
###################################################################

#directory and file listings
LANES_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/"
LANE_READGROUP_APPEND="outs/gex_possorted_bam.bam"
BARCODE_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/souporcell/barcodes/"
GENOTYPES_LOC=""
OUTPUT_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/popscle_tools_filtered/alignment/filtered_alignment/"
JOB_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/popscle_tools_filtered/alignment/jobs/"
FILTER_CMD="/groups/umcg-franke-scrna/tmp02/software/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh"

# check each lane
for dir in "$LANES_DIR"/*lane*/ ; do

    # extrac the lane
    lane_id=$(basename $dir)
    # where we will place the job file
    output_job=${JOB_DIR}"/filter_bam_"${lane_id}"_SBATCH.sh"
    # where the unzipped barcodes are
    barcode_loc=${BARCODE_DIR}"/"${lane_id}"_barcodes.tsv"
    # where the original alignment file is
    original_bam_loc=${LANES_DIR}${lane_id}/${LANE_READGROUP_APPEND}
    # where we will place the filtered sequence file
    filtered_bam_loc=${OUTPUT_DIR}"/mo_"${lane_id}"_genofiltered_gex_possorted_bam.bam"

    echo -e "#!/usr/bin/env bash
#SBATCH --job-name=bam_filter_${lane_id}
#SBATCH --output=bam_filter_${lane_id}.out
#SBATCH --error=bam_filter_${lane_id}.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
	set -e
        
        ml SAMtools
        ml BEDTools
        
        ${FILTER_CMD} \\
                ${original_bam_loc} \\
                ${barcode_loc} \\
                ${GENOTYPES_LOC} \\
                ${filtered_bam_loc}
        " > ${output_job}
done
