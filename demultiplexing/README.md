### demultiplexing and doublet detection
'*demultiplexing/mo_subset_for_genotype_correlations.sh*'   script to subset genotypes for each 10x experiment, using the individuals in present in that experiment\
'*demultiplexing/mo_subset_for_genotype_correlations.sh*'   script to subset genotypes for each 10x experiment, using the individuals in present in that experiment as job\
'*demultiplexing/mo_subset_gex_bams_snp_regions.sh*'    script to create jobs to subset the gene expression alignment files, to only have reads overlapping variants we genotyped\
'*demultiplexing/mo_create_popscle_sort_vcfs_jobs.sh*' s    cript to create jobs to sort the per-lane VCF files by the order of chromosomes in the alignment files\
'*demultiplexing/mo_create_demuxlet_jobs.sh*'   script to create jobs to run Demuxlet on the gene expression alignments\
'*demultiplexing/mo_create_souping_samples_jobs_externalbarcodes.sh*'   script to create jobs that perform Souporcell on each 10x lane while allowing external barcode files\
'*demultiplexing/mo_correlate_genotypes.R*'     correlate the Souporcell cluster genotypes to the genotypes generated for the individuals, to do sample assignment\
'*demultiplexing/mo_plot_demultiplexing_assignments.Rmd*'   plot the demultiplexing assignments\
'*demultiplexing/mo_sample_assignment.R*'   do sample assignment on the Seurat object, based on the correlated souporcell genotypes\
'*demultiplexing/mo_get_rna_qced_barcodes.R*'   extract the valid barcodes in the Seurat object

