#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_arc_metadata.R
# Function: get cellranger arc metadata
############################################################################################################################


####################
# libraries        #
####################

# none


####################
# Functions        #
####################

get_arc_metadata <- function(cellranger_loc, lanes, metadata_append='outs/per_barcode_metrics.csv') {
  # store per lane first
  metadata_per_lane <- list()
  # check each lane
  for (lane in lanes) {
    # paste the path together
    metadata_loc <- paste(cellranger_loc, '/', lane, '/', metadata_append, sep = '')
    # read the table
    metadata_table <- read.table(metadata_loc, sep = ',', header = T)
    # add lane as column
    metadata_table[['lane']] <- lane
    # create some metadata, for now, we'll first just store the lane here
    barcodes_short <- gsub('(-\\d+)', '', metadata_table[['barcode']])
    barcodes_lane <- paste(barcodes_short, rep(lane, times = length(barcodes_short)), sep = '_')
    # add short barcode and lane+barcode
    metadata_table[['barcode_1']] <- metadata_table[['barcode']]
    metadata_table[['barcode_lane']] <- barcodes_lane
    metadata_table[['barcode']] <- barcodes_short
    # put in list
    metadata_per_lane[[lane]] <- metadata_table
  }
  # merge all together
  metadata_all <- do.call('rbind', metadata_per_lane)
  rownames(metadata_all) <- metadata_all[['barcode_lane']]
  return(metadata_all)
}


####################
# Main Code        #
####################

# these are the lanes
lanes <- c('230105_lane1', '230105_lane2', '230105_lane3', '230105_lane4',
           '230105_lane5', '230105_lane6', '230105_lane7', '230105_lane8',
           '230112_lane1', '230112_lane2', '230112_lane3', '230112_lane4',
           '230112_lane5', '230112_lane6', '230112_lane7', '230112_lane8',
           '230120_lane1', '230120_lane2', '230120_lane3', '230120_lane4',
           '230120_lane5', '230120_lane6', '230120_lane7', '230120_lane8',
           '230127_lane1', '230127_lane2', '230127_lane3', '230127_lane4',
           '230127_lane5', '230127_lane6', '230127_lane7', '230127_lane8',
           '230202_lane1', '230202_lane2', '230202_lane3', '230202_lane4',
           '230202_lane5', '230202_lane6', '230202_lane7', '230202_lane8',
           '230209_lane1', '230209_lane2', '230209_lane3', '230209_lane4',
           '230209_lane5', '230209_lane6', '230209_lane7', '230209_lane8',
           '230216_lane1', '230216_lane2', '230216_lane3', '230216_lane4',
           '230216_lane5', '230216_lane6', '230216_lane7', '230216_lane8',
           '230223_lane1', '230223_lane2', '230223_lane3', '230223_lane4',
           '230223_lane5', '230223_lane6', '230223_lane7', '230223_lane8',
           '230302_lane1', '230302_lane2', '230302_lane3', '230302_lane4',
           '230302_lane5', '230302_lane6', '230302_lane7', '230302_lane8',
           '230316_lane1', '230316_lane2', '230316_lane3', '230316_lane4',
           '230316_lane5', '230316_lane6', '230316_lane7'
)

# get cellranger loc
cellranger_loc <- '/groups/umcg-franke-scrna/prm02/projects/multiome/processed/joint/alignment/b38/'
# get the metdatadata
arc_metadata <- get_arc_metadata(cellranger_loc, lanes)
# get where to store the metadata
arc_metadata_loc <- gz('/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/metadata/arc_metadata.tsv.gz')
# save the file
write.table(arc_metadata, arc_metadata_loc, sep = '\t', row.names = F, col.names = T)
