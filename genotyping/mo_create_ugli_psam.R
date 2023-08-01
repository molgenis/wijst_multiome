#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_create_ugli_psam.R
# Function: 
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

# metadata file location
metadata_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_age_sex_ugli.tsv'
# read the metadata file
metadata <- read.table(metadata_loc, sep = '\t', header = T)

# create the required psam file
all_original_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_all.psam.original'
all_new_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_all.psam'
all_psam <- read.table(all_original_psam_loc, header = T, sep = '\t', check.names = F, comment.char = '')
# we remove the sex column, as it is empty right now
all_psam[['SEX']] <- 0
# needs to be numeric?
all_psam[['PAT']] <- 0
all_psam[['MAT']] <- 0
# grab the sex from the metadata
all_psam_sex <- metadata[match(all_psam[['IID']], metadata[['sample']]), 'sex']
# create the empty sex column in the psam
all_psam[['SEX']] <- NA
# we need to change the coding from M/F to 1/2 in the psam
all_psam[!is.na(all_psam_sex) & all_psam_sex == 'male', 'SEX'] <- 1
all_psam[!is.na(all_psam_sex) & all_psam_sex == 'female', 'SEX'] <- 2
all_psam[is.na(all_psam_sex), 'SEX' ] <- 0
# we don't know most of these
all_psam[['Provided_Ancestry']] <- 'EUR'
all_psam[['genotyping_platform']] <- 'GSA-MD-v3'
all_psam[['array_available']] <- 'N'
all_psam[['wgs_available']] <- 'N'
all_psam[['wes_available']] <- 'Y'
all_psam[['age']] <- metadata[match(all_psam[['IID']], metadata[['sample']]), 'age']
all_psam[['age_range']] <- apply(all_psam, 1, function(x){
  floor(as.numeric(x['age'])/10)*10
})
all_psam[['Study']] <- 'wijst_multiome'
all_psam[['smoking_status']] <- NA
all_psam[['hormonal_contraception_use_currently']] <- NA
all_psam[['menopause']] <- NA
all_psam[['pregnancy_status']] <- NA
write.table(all_psam, all_new_psam_loc, sep = '\t', row.names = F, col.names = T, quote = F)
