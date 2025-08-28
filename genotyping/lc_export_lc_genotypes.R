#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lc_export_lc_genotypes.R
# Function: export from the multiomics data, the genotypes that are specific for the LONG-Covid samples
############################################################################################################################


####################
# Main Code        #
####################

# annotation of the original MO annotations to the GSA annotations
realid_assignments_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'
# read the realids
realid_assignments <- read.table(realid_assignments_loc, header = T, sep = '\t')
# fix the sample IDs in the realid table so that if there is no reannotation, the original annotation is in the last column
realid_assignments[['sample_final']] <- realid_assignments$SampleID_CORRECT
realid_assignments[is.na(realid_assignments$SampleID_CORRECT) | realid_assignments$SampleID_CORRECT == '', 'sample_final']  <- realid_assignments[is.na(realid_assignments$SampleID_CORRECT) | realid_assignments$SampleID_CORRECT == '', 'sample']
realid_assignments[['realid_final']] <- realid_assignments$RealID
realid_assignments[is.na(realid_assignments$RealID) | realid_assignments$RealID == '', 'realid_final']  <- realid_assignments[is.na(realid_assignments$RealID) | realid_assignments$RealID == '', 'sample_final']

# annotations of the original MO annotations to the other IDs we have
other_id_assignments_p2_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_origsample_to_realid_p2.tsv'
other_id_assignments_p3_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_origsample_to_realid_p3.tsv'
# read these tables
other_id_assignments_p2 <- read.table(other_id_assignments_p2_loc, header = T, sep = '\t')
other_id_assignments_p3 <- read.table(other_id_assignments_p3_loc, header = T, sep = '\t')
# merge these
other_id_assignments <- rbind(other_id_assignments_p2, other_id_assignments_p3)

# read the annotation to the long-covid ID in the Lanting sheet
origid_to_lcid_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_origsample_to_lcid.tsv'
origid_to_lcid <- read.table(origid_to_lcid_loc, header = T, sep = '\t')
# make sure the sample ID annotation is the same
# origid_to_lcid[['Sample_ID']] <- paste0('MO', origid_to_lcid[['Sample_ID']])
# make sure the sample ID annotation is the same
origid_to_lcid[['Sample_ID']] <- paste0('LONG_', origid_to_lcid[['Sample_ID']])

# rename columns
colnames(other_id_assignments) <- c('realid1', 'sample_id_old')
# get the realid table for the old IDs we don't already have
realid_not_otherid <- realid_assignments[!(realid_assignments$sample %in% other_id_assignments$sample_id_old), c('realid_final', 'sample')]
# set same names here as well
colnames(realid_not_otherid) <- c('realid1', 'sample_id_old')
# add to the other assignments
other_id_assignments <- rbind(other_id_assignments, realid_not_otherid)
# add GSA ID
other_id_assignments[['GSA_ID']] <- realid_assignments[match(other_id_assignments[['sample_id_old']], realid_assignments[['sample']]), 'sample_final']
other_id_assignments[['realid2']] <- realid_assignments[match(other_id_assignments[['sample_id_old']], realid_assignments[['sample']]), 'realid_final']
other_id_assignments[['realid_final']] <- other_id_assignments[['realid1']]
other_id_assignments[is.na(other_id_assignments[['realid_final']]), 'realid_final'] <- other_id_assignments[is.na(other_id_assignments[['realid_final']]), 'realid2']

# # add the new, correct ID to the lcid table
# origid_to_lcid[['GSA_ID']] <- realid_assignments[match(origid_to_lcid[['Sample_ID']], realid_assignments[['sample']]), 'sample_final']
# # add the realid, based on the original ID
# origid_to_lcid[['realid1']] <- other_id_assignments[match(origid_to_lcid[['Sample_ID']], other_id_assignments[['Sample_ID']]), 'Sample_Code']
# # add the other realid, to check if it is the same
# origid_to_lcid[['realid2']] <- realid_assignments[match(origid_to_lcid[['Sample_ID']], realid_assignments[['sample']]), 'realid_final']
# 
# # add the final real id assignment
# origid_to_lcid[['realid_final']] <- origid_to_lcid[['realid1']]
# # taking the '2' one when we don't have the '1' version
# origid_to_lcid[is.na(origid_to_lcid[['realid_final']]), 'realid_final'] <- origid_to_lcid[is.na(origid_to_lcid[['realid_final']]), 'realid2']
# add 'LC' to the long covid ID
origid_to_lcid[['LongCovid_ID']] <- paste0('LC', origid_to_lcid[['LongCovid_ID']])

# add the GSA ID for origid
origid_to_lcid[['GSA_ID']] <- other_id_assignments[match(origid_to_lcid[['Sample_ID']], other_id_assignments[['realid_final']]), 'GSA_ID']

# write list of samples to subset our genotype data for
lc_gsa_inclusion_list_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/long_covid_export/genotype/mo_lc_gsa_ids.txt'
write.table(origid_to_lcid[!is.na(origid_to_lcid[['GSA_ID']]) & grep('LONG', origid_to_lcid[['Sample_ID']]), 'GSA_ID'], lc_gsa_inclusion_list_loc, row.names = F, col.names = F, quote = F)

# do the conversion
# /groups/umcg-franke-scrna/tmp04/software/plink2_amd/plink2 \
#  --bfile /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/genotype/imputed_hg38_all_anc \
#  --keep ./mo_lc_gsa_ids.txt \
#  --make-bed \
#  --out ./lc_imputed_hg38_all_anc

# read the fam file that was created
fam_plink_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/long_covid_export/genotype/lc_imputed_hg38_all_anc.fam'
fam_plink <- read.table(fam_plink_loc, header = F, sep = '\t')
# remap the second column to have the LC IDs that were required
fam_plink[['V2']] <- origid_to_lcid[match(fam_plink[['V2']], origid_to_lcid[['GSA_ID']]), 'LongCovid_ID']
# overwrite the original file
write.table(fam_plink, fam_plink_loc, row.names = F, col.names = F, sep = '\t', quote = F)

# also save this mapping table we painstakingly created
origid_to_lcid_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_lc_to_gsa_and_realid.tsv.gz'
write.table(origid_to_lcid, gzfile(origid_to_lcid_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_md5_for_file(origid_to_lcid_loc)

# and the other mapping table
gsa_id_to_realid_lcid_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_gsa_to_realid_lcid.tsv.gz'
write.table(other_id_assignments, gzfile(gsa_id_to_realid_lcid_loc), row.names = F, col.names = T, sep = '\t')
mdfiver::create_md5_for_file(gsa_id_to_realid_lcid_loc)
