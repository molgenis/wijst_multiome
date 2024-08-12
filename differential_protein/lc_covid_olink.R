#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: lc_covid_olink.R
# Function: create L2 azimuth specific Seurat objects
############################################################################################################################

####################
# libraries        #
####################

# Load OlinkAnalyze
library(OlinkAnalyze)
# for plotting
library(ggplot2)

####################
# Functions        #
####################

plot_plate_layout <- function(olink_data, well_id_column='WellID', plate_id_column='PlateID', sample_column='SampleID') {
  # add the row
  olink_data[['row']] <- gsub('\\d+', '', olink_data[[well_id_column]])
  olink_data[['column']] <- gsub('[A-Z]+', '', olink_data[[well_id_column]])
  # order the rows in reverse order
  olink_data[['row']] <- factor(olink_data[['row']], levels = unique(olink_data[['row']])[rev(order(unique(olink_data[['row']])))])
  # also order the columns 
  olink_data[['column']] <- factor(olink_data[['column']], levels = as.character(1:max(as.numeric(olink_data[['column']]))))
  # subset to the data we can plot
  olink_data_samples <- olink_data[, c(plate_id_column, sample_column, 'row', 'column')]
  # make a tile plot
  p <- ggplot(data=olink_data_samples, aes(x=column, y=row, fill = .data[[plate_id_column]])) + 
    geom_tile() + geom_text(data = olink_data, mapping = aes(label=.data[[sample_column]])) +
    facet_grid(rows = vars(olink_data[[plate_id_column]]))
  # some options
  return(p)
}


####################
# Main Code        #
####################

# location of the protein data
protein_data_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/olink_protein/OLINK2023-037_INF_EXTENDED_NPX_2024-07-10.csv'
# location of the sample mapping
sample_mapping_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/olink_protein/OV21_00402_linkage_file_olink_20240719.csv'
# the full id table
mo_full_id_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_full_id_table.tsv.gz'
# get the other id table
mo_other_id_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_id_to_otherids.tsv.gz'
# and the 'realids'
mo_realid_assignments_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'
# the age/sex for the mo participant
mo_age_sex_loc <- '/groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_age_sex.tsv.gz'

# load the data
protein_data <- read_NPX(protein_data_loc)
# load the sample data
sample_mapping <- read.table(sample_mapping_loc, header = T, sep = ',')

# load other ID data
mo_full_id <- read.table(mo_full_id_loc, header = T, sep = '\t')
mo_other_id <- read.table(mo_other_id_loc, header = T, sep = '\t')
mo_realid_assignments <- read.table(mo_realid_assignments_loc, header = T, sep = '\t')
# fix the sample IDs in the realid table
mo_realid_assignments[['sample_final']] <- mo_realid_assignments$SampleID_CORRECT
mo_realid_assignments[is.na(mo_realid_assignments$SampleID_CORRECT) | mo_realid_assignments$SampleID_CORRECT == '', 'sample_final']  <- mo_realid_assignments[is.na(mo_realid_assignments$SampleID_CORRECT) | mo_realid_assignments$SampleID_CORRECT == '', 'sample']
mo_realid_assignments[['realid_final']] <- mo_realid_assignments$RealID
mo_realid_assignments[is.na(mo_realid_assignments$RealID) | mo_realid_assignments$RealID == '', 'realid_final']  <- mo_realid_assignments[is.na(mo_realid_assignments$RealID) | mo_realid_assignments$RealID == '', 'sample_final']
# replace spaces in other ID for the covid identifiers
mo_other_id[['Other.ID']] <- gsub(' ', '', mo_other_id[['Other.ID']])

# rename the plate
sample_mapping[['plate']] <- gsub('OV210402EDPL0', 'PlateLayout_Plate', sample_mapping$Release.Rack)
# also plate1 and plate2 are actually plate 4 and 5
sample_mapping[sample_mapping[['plate']] == 'PlateLayout_Plate1', 'plate'] <- 'PlateLayout_Plate4'
sample_mapping[sample_mapping[['plate']] == 'PlateLayout_Plate2', 'plate'] <- 'PlateLayout_Plate5'
# remove the 0 in the well name
sample_mapping[['Release.Position']] <- gsub('(0)(\\d+)$', '\\2', sample_mapping[['Release.Position']])
# add combination of well and plate
sample_mapping[['well_plate']] <- paste(sample_mapping[['plate']], sample_mapping[['Release.Position']], sep =  '_')

# the entries with just a number are MO ones, rename them
protein_data[['SampleID']] <- gsub('^(\\d+)$', 'MO\\1', protein_data[['SampleID']])
# get more descriptive ID for those
#protein_data_realid <- mo_realid_assignments[match(protein_data[['SampleID']], mo_realid_assignments[['sample_final']]), 'realid_final']
protein_data_realid <- mo_realid_assignments[match(protein_data[['SampleID']], mo_realid_assignments[['sample']]), 'realid_final']
# replace the sample name where possible
protein_data[!is.na(protein_data_realid), 'realid'] <- protein_data_realid[!is.na(protein_data_realid)]
# if the ID already starts with ll-next, we can just copy that
protein_data[!is.na(protein_data[['SampleID']]) & startsWith(protein_data[['SampleID']], 'LL-NEXT'), 'realid'] <- protein_data[!is.na(protein_data[['SampleID']]) & startsWith(protein_data[['SampleID']], 'LL-NEXT'), 'SampleID']

# now do the LLNEXT and COVID19 samples one as well
protein_data_next_id <- mo_other_id[match(protein_data[['SampleID']], mo_other_id[['Other.ID']]), 'PROJECT_PSEUDO_ID_NEW']
# replace the sample name where possible
protein_data[!is.na(protein_data_next_id) & protein_data_next_id != '' & protein_data_next_id != ' ', 'llnextcovid_pseudo'] <- protein_data_next_id[!is.na(protein_data_next_id) & protein_data_next_id != '' & protein_data_next_id != ' ']

# try from the other one we have as well
protein_data_next2_id <- mo_full_id[match(protein_data[['SampleID']], mo_full_id[['COVID_ID']]), 'll_pseudo_id']
# replace the sample name where possible
protein_data[!is.na(protein_data_next2_id) & protein_data_next2_id != '' & protein_data_next2_id != ' ', 'llnextcovid_pseudo'] <- protein_data_next2_id[!is.na(protein_data_next2_id) & protein_data_next2_id != '' & protein_data_next2_id != ' ']

# add the case/control status
protein_data[['case_control']] <- NA
# get combination of well and plate
protein_data_well_plate <- paste(protein_data[['PlateID']], protein_data[['WellID']], sep = '_')
# use that for case/control inference
protein_data[['case_control']] <- sample_mapping[match(protein_data_well_plate, sample_mapping[['well_plate']]), 'Case.Control']
# where we don't have data, the data was pre-pandemic
protein_data[['pandemic']] <- ifelse(is.na(protein_data[['case_control']]), 'prepandemic', 'postpandemic')
# furthermore, if it was pre-pandemic, the case/control status must be case
protein_data[is.na(protein_data[['case_control']]), 'case_control'] <- 'control'

# match the sample ID that we have for the data in the sample mapping
protein_data_sample_id <- sample_mapping[match(protein_data_well_plate, sample_mapping[['well_plate']]), 'project_pseudo_id']
# replace the sample name where possible
protein_data[!is.na(protein_data_sample_id), 'llnextcovid_pseudo'] <- protein_data_sample_id[!is.na(protein_data_sample_id)]

# read the age/sex file
mo_age_sex <- read.table(mo_age_sex_loc, header = T, sep = '\t')
# match to realid
protein_data_sex <- mo_age_sex[match(protein_data[['realid']], mo_age_sex[['sample']]), 'sex']
protein_data_age <- mo_age_sex[match(protein_data[['realid']], mo_age_sex[['sample']]), 'age']
# replace the sample name where possible
protein_data[!is.na(protein_data_sex), 'sex'] <- protein_data_sex[!is.na(protein_data_sex)]
protein_data[!is.na(protein_data_age), 'age'] <- protein_data_age[!is.na(protein_data_age)]
# make the sex lowercase in the sample mapping
sample_mapping[['Gender']] <- tolower(sample_mapping[['Gender']])
# get those as well 
protein_data_sex2 <- sample_mapping[match(protein_data[['SampleID']], sample_mapping[['SampleID']]), 'Gender']
protein_data_age2 <- sample_mapping[match(protein_data[['SampleID']], sample_mapping[['SampleID']]), 'Age']
# again add sex and age
protein_data[!is.na(protein_data_sex2), 'sex'] <- protein_data_sex2[!is.na(protein_data_sex2)]
protein_data[!is.na(protein_data_age2), 'age'] <- protein_data_age2[!is.na(protein_data_age2)]


# plot the plates to see how they line up
plot_plate_layout(protein_data) + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'case_control') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'pandemic') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'sex') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))
plot_plate_layout(protein_data, sample_column = 'age') + scale_fill_manual(values = roycols::get_color_list(protein_data[['PlateID']]))

# now do actual statistical analysis
olink_lmer(protein_data, variable = c('case_control', 'pandemic'), random = c('SampleID'))
