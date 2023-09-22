#!/usr/bin/env bash

###################################################################
#Script Name	  : mo_create_popscle_sort_vcfs_jobs.sh
#Description	  : create SBATCH jobs scripts to sort the VCF files to be in the BAM order
#Author       	: Roy Oelen
###################################################################


#directory and file listings
LANES_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/"
FILTERED_BAM_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/popscle_tools_filtered/alignment/filtered_gex_alignment/"
INDIVIDUAL_GENOTYPES="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/genotype_per_lane/"
OUTPUT_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/genotype_per_lane/"
JOBS_DIR="/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/popscle_tools_filtered/genotype/jobs/"
SORT_CMD="/groups/umcg-franke-scrna/tmp02/software/popscle_helper_tools/sort_vcf_same_as_bam.sh"

# check each lane
for dir in "$LANES_DIR"/*lane*/ ; do

    # extrac the lane
    lane_id=$(basename $dir)
    # where we will place the job file
    output_job=${JOBS_DIR}"/sort_vcf_"${lane_id}"_SBATCH.sh"
    # where incorrectly sorted VCF is
    genotype_filter_loc=${INDIVIDUAL_GENOTYPES}"/mo_"${lane_id}"_maf005.vcf"
    # where the BAM is with our sort error
    filtered_bam_loc=${FILTERED_BAM_DIR}"/mo_"${lane_id}"_genofiltered_gex_possorted_bam.bam"
    # where we place the sorted VCF
    genotype_sorted_loc=${OUTPUT_DIR}"mo_gex_bamsorted_"${lane_id}"_maf005.vcf"
    # and which we will zip
    genotype_sorted_zipped_loc=${OUTPUT_DIR}"mo_gex_bamsorted_"${lane_id}"_maf005.vcf.gz"
    # create a temporary sorting directory
    genotype_sorted_tmp_directory=${OUTPUT_DIR}"mo_gex_bamsorted_"${lane_id}"_maf005_tmp/"

    echo -e "#!/usr/bin/env bash
#SBATCH --job-name=vcf_sort_gex_${lane_id}
#SBATCH --output=vcf_sort_gex_${lane_id}.out
#SBATCH --error=vcf_sort_gex_${lane_id}.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
	set -e

        ml SAMtools
        export PATH="/groups/umcg-franke-scrna/tmp02/software/bedtools/:'$PATH'"
        ml BCFtools
        ml HTSlib

        bgzip -c ${genotype_filter_loc} > ${genotype_filter_loc}.gz
        tabix -p vcf ${genotype_filter_loc}.gz

        mkdir -p ${genotype_sorted_tmp_directory}
        TMPDIR=${genotype_sorted_tmp_directory}

        ${SORT_CMD} \\
                ${filtered_bam_loc} \\
                ${genotype_filter_loc}.gz \\
                > ${genotype_sorted_zipped_loc} \\

        #bgzip -c ${genotype_sorted_loc} > ${genotype_sorted_zipped_loc}

        tabix -p vcf ${genotype_sorted_zipped_loc}

        #rm ${genotype_sorted_loc}

        rm -r ${genotype_sorted_tmp_directory}

        " > ${output_job}
done
