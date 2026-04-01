#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_scenic_to_sceqtlgen_overlap.R
# Function: compare SCENIC+ TF-gene relations to sc-eQTLgen cis-trans pairs
#
############################################################################################################################

####################
# libraries        #
####################


####################
# Functions        #
####################

extract_cis_trans_pairs <- function(cistrans_table, cis_column='phenotype', trans_column='sc_cis_ct_gene_list', trans_sep=', ', ct_gene_sep=';') {
  # check each row of the cistrans table
  cistrans_ct <- apply(cistrans_table, 1, function(x) {
    # extract the list
    sc_cistrans <- x[[trans_column]]
    # check if this is NA
    if (is.na(sc_cistrans)) {
      # then return nothing
      return(NULL)
    }
    else {
      # otherwise try to split
      sc_cistrans_split <- strsplit(sc_cistrans, trans_sep)[[1]]
      # make a list to keep all the cell types
      cistrans_per_ct <- list()
      # extract the cis gene
      cis_gene <- x[[cis_column]]
      # now check each of the entries
      for (cistrans_ct in unique(sc_cistrans_split)) {
        # now split this into cell type and gene
        ct_to_gene <-  strsplit(cistrans_ct, ct_gene_sep)[[1]]
        # get explicit columns
        ct <- ct_to_gene[[1]]
        gene <- ct_to_gene[[2]]
        # make into a table
        ct_to_gene_row <- data.frame('celltype' = c(ct), 'cisgene' = c(cis_gene), 'transgene' = c(gene))
        # put in the list
        cistrans_per_ct[[paste(cis_gene, ct, gene)]] <- ct_to_gene_row
      }
      # combined the rows for this row
      cistrans_all_cis <- do.call('rbind', cistrans_per_ct)
      return(cistrans_all_cis)
    }
  })
  # merge all of them
  cistrans_ct_table <- do.call('rbind', cistrans_ct)
  return(cistrans_ct_table)
}


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

# location of cis-trans outputs:
cistrans_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/overlaps/sceqtlgen_cistrans_pairs.tsv.gz'
# read the table
cistrans <- fread(cistrans_loc, header = T, sep = '\t')
# we're only interested where there is a single-cell overlap
cistrans <- cistrans[!is.na(sc_cis_ct_gene_list), ]
# now split the columns
cistrans_split <- extract_cis_trans_pairs(cistrans)
# write this to a new file
cistrans_split_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/overlaps/sceqtlgen_cistrans_pairs_split.tsv.gz'
write.table(cistrans_split, gzfile(cistrans_split_loc), row.names = F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file(cistrans_split_loc)

# add the cis to trans link in the the eQTLgen output
cistrans_split[['cis_to_trans']] <- paste(cistrans_split[['cisgene']], cistrans_split[['transgene']])
# merge them
scenic_to_sceqtlgen <- merge(x = unique(scenic_output[, c('tf_to_gene_ens', 'rho_TF2G', 'TF', 'Gene')]), y = cistrans_split[, c('cis_to_trans', 'celltype', 'cisgene', 'transgene')], by.x = 'tf_to_gene_ens', by.y = 'cis_to_trans')
# write the result
cistrans_split_to_sceqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/overlaps/mo_scenic_to_sceqtlgen_cistrans_pairs.tsv.gz'
write.table(scenic_to_sceqtlgen, gzfile(cistrans_split_to_sceqtl_loc), , row.names = F, col.names = T, sep = '\t', quote = F)
# read the sign
cistrans_sign_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/overlaps/sc_cis_trans_pairs_clean_V3.txt'
cistrans_sign <- read.table(cistrans_sign_loc, header = T, sep = '\t')
# column names
colnames(cistrans_sign) <- c('transgene', 'cisgene', 'direction')
# save in proper format
write.table(cistrans_sign, gzfile('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/overlaps/sc_cis_trans_pairs_clean_V3_properformat.tsv.gz'), row.names =F, col.names = T, sep = '\t', quote = F)
mdfiver::create_sha256_for_file('/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/overlaps/sc_cis_trans_pairs_clean_V3_properformat.tsv.gz')
# add cistrans again
cistrans_sign[['cis_to_trans']] <- paste(cistrans_sign[['cisgene']], cistrans_sign[['transgene']])
scenic_to_sceqtlgen <- merge(x = scenic_to_sceqtlgen, y = cistrans_sign[, c('cis_to_trans', 'direction')], by.x = 'tf_to_gene_ens', by.y = 'cis_to_trans')
# write the result
cistrans_split_to_sceqtl_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/overlaps/mo_scenic_to_sceqtlgen_cistrans_pairs.tsv.gz'
write.table(scenic_to_sceqtlgen, gzfile(cistrans_split_to_sceqtl_loc), , row.names = F, col.names = T, sep = '\t', quote = F)

