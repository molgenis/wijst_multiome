#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_annotate_gwas_variants.R
# Function: use the immune GWAS data, and add the SNP format in the chrom-pos-alt-ref format
# 
############################################################################################################################

####################
# libraries        #
####################

# table format
library(data.table)


####################
# Functions        #
####################



####################
# Settings         #
####################

# luck seed
set.seed(7777)
# whether we are in debug mode
debug <- T


####################
# Main code        #
####################

# location of the GWAS file
gwas_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full.tsv'
# location of 1000g rsID to chrom:pos:alt:ref
rsid_to_chrompos_loc <- '/groups/umcg-franke-scrna/tmp04/external_datasets/sc-eqtlgen-imputation-ref-hg38/ref_panel_QC/30x-GRCh38-rs-to-alleles.tsv.gz'

# read these two tables
gwas <- fread(gwas_loc, header = T, sep = '\t')
rsid_to_chrompos <- fread(rsid_to_chrompos_loc, header = T, sep = '\t')

# add a new column that converts the rsid to chrom-pos-alt-ref format
gwas[['chromposaltref']] <- rsid_to_chrompos[match(gwas[['SNPS']], rsid_to_chrompos[['rsid']]), ][['chromposallele']]
# where it is NA, set back the original
gwas[is.na(gwas[['chromposaltref']]), ][['chromposaltref']] <- gwas[is.na(gwas[['chromposaltref']]), ][['SNPS']]

# write the result
gwas_chromposaltref_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/GWAS_enrichment/GWAS_vars/immune-gwas-catalog-download-associations-alt-full-chromposrefalt.tsv.gz'
write.table(gwas, gzfile(gwas_chromposaltref_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_sha256_for_file(gwas_chromposaltref_loc)
