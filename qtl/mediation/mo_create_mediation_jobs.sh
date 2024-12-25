#!/bin/bash

###################################################################
#Script Name	  : mo_create_mediation_jobs.sh
#Description	  : create SBATCH jobs scripts to do mediation per cell type object
#Args           :
#Author       	: Roy Oelen
#example        : ./mo_create_mediation_jobs.sh
###################################################################

# standard parameters
CORES='2'
MEMORY_GB='64'
TMP_SIZE='512mb'
RUNTIME='23:59:59'

# chromosomes
CHROMS=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22')
# cell types
CELL_TYPES=('B' 'CD4T' 'CD8T' 'DC' 'monocyte' 'NK')

# genotype file prepend
GENO_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr'
# expression file prepend before cell type
EXP_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/sc-eqtlgen/input/L1/UT/'
# expression file append after cell type
EXP_APPEND='.qtlInput.txt.gz'
# accessibility file prepend before cell type
ATAC_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/'
# accessibility file append after cell type
#ATAC_APPEND='.qtlInput.txt.gz'
ATAC_APPEND='.qtlInput.PcCorrectedResiduals.txt.gz'
# confinement file prepend before cell type
CONFINE_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/confinements/'
# confinement file append after cell type
CONFINE_APPEND='.confinement.tsv.gz'
# metadata file prepend before cell type
METADATA_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/input/atac_to_expression/metadata/'
# metadata file append after cell type
METADATA_APPEND='.metadata.tsv.gz'
# output file prepend before cell type
#OUT_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/output/atac_to_expression/UT/'
OUT_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/output/regressed_atac_to_gausnormed_expression/UT/'
# output file append after cell type
OUT_APPEND='.tsv.gz'
# output file prepend before cell type
#JOB_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/jobs/atac_to_expression/UT/'
JOB_PREPEND='/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/mediation/jobs/regressed_atac_to_gausnormed_expression/UT/'
# output file append after cell type
JOB_APPEND='_SBATCH.sh'

# location to the R script
SCRIPT_LOC='/groups/umcg-franke-scrna/tmp04/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_perform_qtl_mediation_analysis.R'

# R command
R_COMMAND='~/start_Rscript.sh'

# fixed effects string
FIXED_EFFECTS='RNA_UT_PC1,RNA_UT_PC2,RNA_UT_PC3,RNA_UT_PC4,RNA_UT_PC5,RNA_UT_PC6,RNA_UT_PC7,RNA_UT_PC8,RNA_UT_PC9,RNA_UT_PC10'
# random effects string
RANDOM_EFFECTS='donor'

# do each cell types
for cell_type in ${CELL_TYPES[*]}
do
      # paste together paths
      exp_loc=${EXP_PREPEND}${cell_type}${EXP_APPEND}
      acc_loc=${ATAC_PREPEND}${cell_type}${ATAC_APPEND}
      mtdt_loc=${METADATA_PREPEND}${cell_type}${METADATA_APPEND}
      confine_loc=${CONFINE_PREPEND}${cell_type}${CONFINE_APPEND}
      # and each chromosomes
      for chrom in ${CHROMS[*]}
      do
        # paste other paths
        geno_loc=${GENO_PREPEND}${chrom}
        out_loc=${OUT_PREPEND}${cell_type}'_chr'${chrom}${OUT_APPEND}
        job_loc=${JOB_PREPEND}${cell_type}'_chr'${chrom}${JOB_APPEND}
        job_name=${cell_type}'_chr'${chrom}
        job_out=${JOB_PREPEND}${job_name}'.out'
        job_err=${JOB_PREPEND}${job_name}'.err'
        # echo the header
        echo '#!/bin/bash
#SBATCH --job-name='${job_name}'
#SBATCH --output='${job_out}'
#SBATCH --error='${job_err}'
#SBATCH --time='${RUNTIME}'
#SBATCH --cpus-per-task='${CORES}'
#SBATCH --mem='${MEMORY_GB}'GB
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp='${TMP_SIZE}'

'> ${job_loc}
        # build the command
        run_command=${R_COMMAND}'
        '${SCRIPT_LOC}'
        --eqtl_file '${exp_loc}'
        --caqtl_file '${acc_loc}'
        --genotype_file '${geno_loc}'
        --metadata '${mtdt_loc}'
        --fixed_effects '${FIXED_EFFECTS}'
        --random_effects '${RANDOM_EFFECTS}'
        --confinement_list '${confine_loc}'
        --out '${out_loc}'
        --eqtl_gausnorm'

        # add to the job
        echo ${run_command} >> ${job_loc}
      done
done
