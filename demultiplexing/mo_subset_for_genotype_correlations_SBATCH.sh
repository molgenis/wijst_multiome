#!/bin/bash
#SBATCH --job-name=mo_subset_for_genotype_correlations
#SBATCH --output=mo_subset_for_genotype_correlations.out
#SBATCH --error=mo_subset_for_genotype_correlations.err
#SBATCH --time=167:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp=512mb

###################################################################
#Script Name	  : mo_subset_for_genotype_correlations_SBATCH.sh
#Description	  : subset the genotype file to have the genotypes per lane for demultiplexing
#Args           :
#Author       	: Roy Oelen
###################################################################

# the location of the full genotype file
FULL_GENO_FILE='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/imputed_hg38_all_anc_mmaf005_chrprepend.vcf.gz'
# where we want the subsetted genotypes
GENOTYPE_LANE_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/genotype/genotype_per_lane/'
# the prepend of the genotype files
GENOTYPE_LANE_PREPEND='mo_'
# the extentions of the genotype files per lane
GENOTYPE_LANE_UNFILTERED_APPEND='.vcf.gz'
GENPTYPE_LANE_FILTERED_APPEND='_maf005.vcf'
# where the file is that has the samples per lane
SAMPLES_LANE_LOC='/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/participant_per_lane/'
# the extention for the per-lane participant list
SAMPLES_LANE_APPEND='.txt'
# the lanes to consider
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

# load required libraries
ml BCFtools
ml HTSlib

# check each run
for lane in ${LANES[*]}
  do

    echo 'starting '${lane}
    # subset to just the samples of that lane
    bcftools view \
      -S ${SAMPLES_LANE_LOC}/${lane}${SAMPLES_LANE_APPEND} \
      -O z \
      -o ${GENOTYPE_LANE_LOC}/${GENOTYPE_LANE_PREPEND}${lane}${GENOTYPE_LANE_UNFILTERED_APPEND} \
      --force-samples \
      --threads 8 \
      ${FULL_GENO_FILE}
    # index the file
    tabix -p vcf ${GENOTYPE_LANE_LOC}/${GENOTYPE_LANE_PREPEND}${lane}${GENOTYPE_LANE_UNFILTERED_APPEND}
    # subset to SNPs with a minor allele frequency of 0.05 or more
    bcftools view \
      -q 0.05:minor \
      -O v \
      -o ${GENOTYPE_LANE_LOC}/${GENOTYPE_LANE_PREPEND}${lane}${GENPTYPE_LANE_FILTERED_APPEND} \
      --threads 8 \
      ${GENOTYPE_LANE_LOC}/${GENOTYPE_LANE_PREPEND}${lane}${GENOTYPE_LANE_UNFILTERED_APPEND}

    echo 'finished '${lane}

done
