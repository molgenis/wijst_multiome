#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_region_topics_to_beds.R
# Function: take the topic membership table, and convert it to bed files to scenicplus
############################################################################################################################

####################
# libraries        #
####################

library(stringr)
library(mdfiver)


####################
# Functions        #
####################


region_memberships_to_bed_formats <- function(region_membership_table, region_column='region', topic_columns_pattern='^Topic\\d+', true_value='True') {
  # get the columns in the table
  topic_columns <- colnames(region_membership_table)
  # subset to the ones with the pattern
  topic_columns <- topic_columns[grep(topic_columns_pattern, topic_columns)]
  # and remove the region column
  topic_columns <- setdiff(topic_columns, region_column)
  # we'll store the table per topic
  bed_format_per_topic <- list()
  # check each the topics
  for (topic in topic_columns) {
    # extract the regions that are for this topic
    topic_regions <- region_membership_table[!is.na(region_membership_table[[topic]]) &region_membership_table[[topic]] == true_value, region_column]
    # make into chromosome, start and end
    topic_bed_format <- stringr::str_split_fixed(topic_regions, ':|-', 3)
    # turn into a dataframe
    topic_bed_format <- data.frame(topic_bed_format)
    # rename columns
    colnames(topic_bed_format) <- c('chromosome', 'start', 'end')
    # and add the name
    topic_bed_format[['name']] <- topic_regions
    # put in the list
    bed_format_per_topic[[topic]] <- topic_bed_format
  }
  return(bed_format_per_topic)
}

####################
# Main Code        #
####################

# the location of the region-to-topic
region_to_topic_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/models/all_nuclei_and_regions_major_minor/mo_topic_region_membership_120.tsv.gz'

# and where to place the bed files
region_to_topic_beds_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/topic_memberships/region_memberships/beds/'

# read the region to topic file
region_to_topic <- read.table(region_to_topic_loc, header = T, sep = '\t')

# turn into bed format
bed_to_topic <- region_memberships_to_bed_formats(region_to_topic)

# now write each topic
for (topic in names(bed_to_topic)) {
  # paste the output path together
  topic_bed_full_loc <- paste0(region_to_topic_beds_loc, '/', topic, '.bed')
  # write the table
  write.table(bed_to_topic[[topic]], topic_bed_full_loc, row.names = F, col.names = F, sep = '\t', quote = F)
  # and make a checksum
  mdfiver::create_md5_for_file(topic_bed_full_loc)
}
