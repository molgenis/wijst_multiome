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
protein_data_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/olink_protein/OLINK2023-037_INF_EXTENDED_NPX_2024-07-10.csv'
# location of the sample mapping
sample_mapping_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/olink_protein/OV21_00402_linkage_file_olink_20240719.csv'
# the full id table
mo_full_id_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_full_id_table.tsv.gz'
# get the other id table
mo_other_id_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_id_to_otherids.tsv.gz'
# and the 'realids'
mo_realid_assignments_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/metadata/mo_sample_sheet_final.tsv'

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
protein_data[!is.na(protein_data_realid), 'SampleID'] <- protein_data_realid[!is.na(protein_data_realid)]

# get combination of well and plate
protein_data_well_plate <- paste(protein_data[['PlateID']], protein_data[['WellID']], sep = '_')
# match the sample ID that we have
protein_data_sample_id <- sample_mapping[match(protein_data_well_plate, sample_mapping[['well_plate']]), 'project_pseudo_id']
# replace the sample name where possible
protein_data[!is.na(protein_data_sample_id), 'SampleID'] <- protein_data_sample_id[!is.na(protein_data_sample_id)]

# now do the LLNEXT and COVID19 samples one as well
protein_data_next_id <- mo_other_id[match(protein_data[['SampleID']], mo_other_id[['Other.ID']]), 'PROJECT_PSEUDO_ID_NEW']
# replace the sample name where possible
protein_data[!is.na(protein_data_next_id) & protein_data_next_id != '' & protein_data_next_id != ' ', 'SampleID'] <- protein_data_next_id[!is.na(protein_data_next_id) & protein_data_next_id != '' & protein_data_next_id != ' ']

# try from the other one we have as well
protein_data_next2_id <- mo_full_id[match(protein_data[['SampleID']], mo_full_id[['COVID_ID']]), 'll_pseudo_id']
# replace the sample name where possible
protein_data[!is.na(protein_data_next2_id) & protein_data_next2_id != '' & protein_data_next2_id != ' ', 'SampleID'] <- protein_data_next2_id[!is.na(protein_data_next2_id) & protein_data_next2_id != '' & protein_data_next2_id != ' ']

# plot the plates to see how they line up
plot_plate_layout(protein_data)
