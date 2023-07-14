#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Maryna Korshevniuk
# Name: mo_collect_previously_genotyped_individuals.R
# Function: collect the IDs of the individuals that were genotyped before in other studies
############################################################################################################################

####################
# libraries        #
####################

####################
# Functions        #
####################

####################
# Main Code        #
####################

# full ID table of all NEXT participants
all_next_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LL_multiome/LL_ids_from_excel.tsv'
all_next <- read.table(all_next_loc, ,sep='\t',  fill=T, header=T)
# multiome participants with NEXT IDs, 71 people, a subset of all NEXT participants
multiome_next_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LL_multiome/Final-LL-NEXT-DAG2-GSA_info_multiome.tsv'
multiome_next <- read.table(multiome_next_loc, sep = '\t', header = T,  colClasses = rep('character', times = 9))
# make the colnames a bit better
colnames(multiome_next) <- c('achternaam', 'voorletters', 'gebdat', 'NEXT_ID', 'LL_ID', 'COVID_ID', 'ronde', 'DNA_avail', 'GSA_avail')
# add the genotypeing platform as it has been described in the all_next file
multiome_next[['next_genotype_platform']] <- all_next[match(multiome_next[['NEXT_ID']], all_next[['NEXT.ID']]), 'Genotyping.platform'] # we observe that the 'Genotyping.platform' and 'GSA_avail' are not always concordant
# this file has a mapping of the NEXT ID to the hash 10+ character ID of LifeLines
next_to_llhashid_loc <- '/groups/umcg-lifelines/tmp01/projects/ov20_0554/umcg-mkorshevniuk/ll_next_pairing/LL_NEXT_Pairing_File.csv'
next_to_llhashid <- read.table(next_to_llhashid_loc, sep = ',', header = T, colClasses = rep('character', times = 3))
# we need to pad zeroes to the beginning of the NEXT IDs
next_to_llhashid[['next_id_zeroed']] <- apply(next_to_llhashid, 1, function(x){
  # extract the value from the right column
  unpadded_id <- x['NEXT_ID']
  # we need six characters, so the current length minus that is how much to pad
  nr_to_pad <- 6 - nchar(unpadded_id)
  # pad some zeroes
  pad <- paste(rep('0', times = nr_to_pad), collapse = '')
  # now add that to the original
  padded_id <- paste(pad, unpadded_id, sep = '')
  return(padded_id)
})
# now we will add the hashid to the multiome table
multiome_next[['ll_pseudo_id']] <- next_to_llhashid[match(multiome_next[['NEXT_ID']], next_to_llhashid[['next_id_zeroed']]), 'PROJECT_PSEUDO_ID']
# read the table that maps the pseudo ll id to the 8-letter ll ID
pseudo_int_to_ext_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/phenotype_linkage_file_project_pseudo_id.txt'
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T)
# and add it to the multiome table
multiome_next[['ll_pseudo_ext']] <- pseudo_int_to_ext[match(multiome_next[['ll_pseudo_id']], pseudo_int_to_ext[['PROJECT_PSEUDO_ID']]), 'PSEUDOIDEXT']
# read the UGLI mapping
pseudo_int_to_ugli_mapping_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'
pseudo_int_to_ugli_mapping <- read.table(pseudo_int_to_ugli_mapping_loc, sep = '\t', header = T)
# then add the UGLI ID
multiome_next[['ugli_id']] <- pseudo_int_to_ugli_mapping[match(multiome_next[['ll_pseudo_ext']], pseudo_int_to_ugli_mapping[['PSEUDOIDEXT']]), 'UGLI_ID']
# read the cytoSNP mapping
pseudo_ext_to_cytoid_mapping_loc <- '/groups/umcg-lifelines/tmp01/releases/cytosnp_linkage_files/v4/cytosnp_linkage_file.dat'
pseudo_ext_to_cytoid_mapping <- read.table(pseudo_ext_to_cytoid_mapping_loc, sep = '\t', header = T)
# then add the cytoSNP ID
multiome_next[['cyto_id']] <- pseudo_ext_to_cytoid_map[match(multiome_next[['ll_pseudo_ext']], multiome_next[['PSEUDOIDEXT']]), 'cytosnp_ID']
# now check the number of entries for each genotyping
nr_next_genotyped_soesma <- nrow(multiome_next[!is.na(multiome_next[['GSA_avail']]) & multiome_next[['GSA_avail']] != '', ])
# 10
nr_next_genotyped_gsaplatform <- nrow(multiome_next[!is.na(multiome_next[['next_genotype_platform']]) & multiome_next[['next_genotype_platform']] != '', ])
# 27
# if you compare the Soesma list to the full NEXT file, it seems that what Soesma noted as present, is always Version 3.0
multiome_next[!is.na(multiome_next[['next_genotype_platform']]) & multiome_next[['next_genotype_platform']] != '', c('GSA_avail', 'next_genotype_platform')]
multiome_next[!is.na(multiome_next[['GSA_avail']]) & multiome_next[['GSA_avail']] == 'Yes', c('GSA_avail', 'next_genotype_platform')]
# on to checking numbers
nr_next_ugli_genotyped <- nrow(multiome_next[!is.na(multiome_next[['ugli_id']]) & multiome_next[['ugli_id']] != '', ])
# 17
# if you compare the full NEXT file to the UGLI, you see that where there is UGLI available, it is also always Version 1.0 in full NEXT file
nr_next_cytosnp_genotyped <- nrow(multiome_next[!is.na(multiome_next[['cyto_id']]) & multiome_next[['cyto_id']] != '', ])
# check how many are genotyped at all
nr_next_genotyped_atall <- nrow(multiome_next[
  (!is.na(multiome_next[['cyto_id']]) & multiome_next[['cyto_id']] != '') |
    (!is.na(multiome_next[['ugli_id']]) & multiome_next[['ugli_id']] != '') |
    (!is.na(multiome_next[['next_genotype_platform']]) & multiome_next[['next_genotype_platform']] != '')
  , ])
# 27