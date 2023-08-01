# Wijst Multi-ome study

This is the github repository of the Multiome study containing CA stimulated PBMCs

## Software

Here we will list the software used to generate the data

Cellranger 7.1.0\
Cellranger-arc 2.0.0\
R 4.3.1\
plink2-20230707

## Custom code

Below we will outline the steps taken to process and analyse the data

### alignment to reference genome

'*alignment/mo_create_cellranger_jobs_batch1.sh*' create SLURM jobs to align the RNA-seq data to the b38 human reference\
'*alignment/mo_create_multiome_csvs_batch1.sh*' create SLURM jobs to do a joint align of the RNA-seq and ATAC-seq data to the b38 human reference

### genotypeing

'*genotyping/mo_collect_previously_genotyped_individuals.R*' collect data on which individuals were previously genotyped\
'*genotyping/mo_get_previous_genotypes.R*' create list of previously genotyped individuals, create sample sheets, create age/sex metadata file\
'*genotyping/mo_preprocess_previous_genotypes.sh*' subset previously genotyped individuals genotype data, and convert to format for imputation\
'*genotyping/mo_create_ugli_psam.R*' create psam file required for imputation, for the previously genotyped individuals\
'*genotyping/PreImputation_ugli.yaml*' the imputation config used for imputing the previously genotyped UGLI individuals
