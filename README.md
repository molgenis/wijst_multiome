# Wijst Multi-ome study

This is the github repository of the Multiome study containing CA stimulated PBMCs

## Software

Here we will list the software used to generate the data

Cellranger 7.1.0\
Cellranger-arc 2.0.0\
R 4.3.1\
plink2-20230707\
CellBender v3\
Seurat 4.9.9.9059\
Souporcell v2\
Demuxlet v2

## Custom code

Below we will outline the steps taken to process and analyse the data

### alignment to reference genome

'*alignment/mo_create_cellranger_jobs_batch1.sh*' create SLURM jobs to align the RNA-seq data to the b38 human reference\
'*alignment/mo_create_multiome_csvs_batch1.sh*' create SLURM jobs to do a joint align of the RNA-seq and ATAC-seq data to the b38 human reference\
'*alignment/mo_create_cellranger_arc_jobs_batch1.sh*' create SLURM jobs to align the RNA+ATAC data to the b38 human reference\
'*alignment/mo_download_websummaries.sh*' download the CellRanger web summary files\
'*alignment/mo_get_cellranger_summaries.R*' created an Excel sheet describing the CellRanger results, from the CellRanger output directories

### genotypeing

'*genotyping/mo_collect_previously_genotyped_individuals.R*' collect data on which individuals were previously genotyped\
'*genotyping/mo_get_previous_genotypes.R*' create list of previously genotyped individuals, create sample sheets, create age/sex metadata file\
'*genotyping/mo_preprocess_previous_genotypes.sh*' subset previously genotyped individuals genotype data, and convert to format for imputation\
'*genotyping/mo_create_ugli_psam.R*' create psam file required for imputation, for the previously genotyped individuals\
'*genotyping/mo_create_batch2_psam.R*' create psam file required for imputation, for the second batch of individuals\
'*genotyping/PreImputation_ugli.yaml*' the imputation config used for imputing the previously genotyped UGLI individuals\
'*genotyping/PreImputation_mo_batch2.yaml*' the imputation config used for imputing the second batch of new individuals

### ambient RNA correction
'*ambient_rna_correction/mo_cellbend_rna.sh*' create sbatch jobs to do ambient RNA correction using CellBender on the RNA data

### demultiplexing and doublet detection
'*demultiplexing/mo_run_scrublet.py*' python script to run Scrublet on a CellBender corrected 10x lane\
'*demultiplexing/mo_test_scrublet.ipynb*' jupyter notebook to test Scrublet on individual 10x lanes\
'*demultiplexing/mo_subset_for_genotype_correlations.sh*' script to subset genotypes for each 10x experiment, using the individuals in present in that experiment\
'*demultiplexing/mo_subset_gex_bams_snp_regions.sh*' script to create jobs to subset the gene expression alignment files, to only have reads overlapping variants we genotyped\
'*demultiplexing/mo_create_popscle_sort_vcfs_jobs.sh*' script to create jobs to sort the per-lane VCF files by the order of chromosomes in the alignment files
'*demultiplexing/mo_create_demuxlet_jobs.sh*' script to create jobs to run Demuxlet on the gene expression alignments
'*demultiplexing/mo_create_souping_samples_jobs.sh*' script to create jobs that perform Souporcell on each 10x lane\
'*demultiplexing/mo_correlate_genotypes.R*' correlate the Souporcell cluster genotypes to the genotypes generated for the individuals, to do sample assignment

### scanpy preprocess
'*scanpy_preprocess/mo_h5_to_scanpy_objects.ipynb*' read the Seurat-incompatible h5 CellBender outputs into Scanpy H5AD objects

### seurat preprocess
'*seurat_preprocess/mo_scanpy_to_h5seurat.R*' covert the scanpy H5AD objects into Seurat compatible h5seurat objects
