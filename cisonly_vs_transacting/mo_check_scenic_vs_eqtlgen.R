#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_check_scenic_vs_eqtlgen.R
# Function: compare SCENIC+ to eQTLgen cis-trans pairs
#
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(mdfiver)


####################
# Functions        #
####################


####################
# Settings         #
####################

# luck seed
set.seed(7777)


####################
# Main code        #
####################

# location of the CREs identified by SCENIC
scenic_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/eRegulon_both.tsv.gz'
# location of trans-eQTLs
eqtlgen_trans_pairs <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/replication/eqtlgen/PositivesBasedOnColocalization.txt.gz'

# read scenic output
scenic_output <- fread(scenic_output_loc, header = T, sep = '\t')
# read the eqtlgen output
eqtlgen_output <- fread(eqtlgen_trans_pairs, sep = '\t', header = F)
# set correct header
colnames(eqtlgen_output) <- c('CisGene', 'TransGene', 'Direction', 'NrCredibleSetsAffectingBothCisAndTrans', 'SNP1', 'SNP2')

# read the location of the genes
gene_anno_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz'
gene_anno <- fread(gene_anno_loc, header = F, sep = '\t')
# add columns
colnames(gene_anno) <- c('ens', 'gs', 'modality', 'chrom', 'start', 'end')

# add gene symbols to the genes in scenic
scenic_output[['gene_ens']] <- gene_anno[match(scenic_output[['Gene']], gene_anno[['gs']]), ][['ens']]
scenic_output[['tf_ens']] <- gene_anno[match(scenic_output[['TF']], gene_anno[['gs']]), ][['ens']]
# and combine the tf to gene
scenic_output[['tf_to_gene_ens']] <- paste(scenic_output[['tf_ens']], scenic_output[['gene_ens']])

# add the cis to trans link in the the eQTLgen output
eqtlgen_output[['cis_to_trans']] <- paste(eqtlgen_output[['CisGene']], eqtlgen_output[['TransGene']])

# merge them
scenic_to_eqtlgen <- merge(x = unique(scenic_output[, c('tf_to_gene_ens', 'rho_TF2G', 'TF', 'Gene')]), y = eqtlgen_output[, c('cis_to_trans', 'Direction', 'CisGene', 'TransGene')], by.x = 'tf_to_gene_ens', by.y = 'cis_to_trans')
# check concordance
sum(scenic_to_eqtlgen[['Direction']] == sign(scenic_to_eqtlgen[['rho_TF2G']])) / nrow(scenic_to_eqtlgen)

# subset to columns we care about
scenic_to_eqtlgen <- scenic_to_eqtlgen[, c('CisGene', 'TransGene', 'TF', 'Gene', 'Direction', 'rho_TF2G')]
# write the results
scenic_to_eqtlgen_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/replication/eqtlgen/mo_eqtlgen_vs_scenic_tf_cistrans.tsv.gz'
write.table(scenic_to_eqtlgen, gzfile(scenic_to_eqtlgen_loc), row.names = F, col.names = T, sep = '\t', quote = F)
# and make checksum
mdfiver::create_md5_for_file(scenic_to_eqtlgen_loc)
