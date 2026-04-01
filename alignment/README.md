### alignment to reference genome
'*alignment/mo_create_cellranger_jobs_batch1.sh*'   create SLURM jobs to align the RNA-seq data to the b38 human reference\
'*alignment/mo_create_multiome_csvs_batch1.sh*'     create SLURM jobs to do a joint align of the RNA-seq and ATAC-seq data to the b38 human reference\
'*alignment/mo_create_cellranger_arc_jobs_batch1.sh*'   create SLURM jobs to align the RNA+ATAC data to the b38 human reference\
'*alignment/mo_download_websummaries.sh*' download the  CellRanger web summary files\
'*alignment/mo_get_cellranger_summaries.R*'     created an Excel sheet describing the CellRanger results, from the CellRanger output directories\
'*alignment/mo_copy_trailing_cellranger_data.sh*'     copy cellranger folder back to tmp while keeping folder structure
