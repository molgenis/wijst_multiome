#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: mo_get_credible_atac_barcodes.R
# Function: various plot of the MO object
############################################################################################################################

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
           '230316_lane5', '230316_lane6', '230316_lane7', '230316_lane8'
)

# these are all credible barcodes
all_barcodes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/archr_preprocess_samples/rounded_fragments_project/barcodes_after_qc_doublets_removed.csv'
all_barcodes <- read.table(all_barcodes_loc, header = T, comment.char = '')

# here we will output them per lane
credible_barcodes_loc <- '/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/cpeaks_peak_calling/credible_barcodes/'

for (lane in lanes) {
  # progress
  print(lane)
  # get the ones for a lane
  barcodes_lane <- all_barcodes[grepl(lane, all_barcodes$x), , drop = F]
  if (nrow(barcodes_lane) > 0) {
    # remove the lane from the barcode name
    barcodes_lane[['x']] <- gsub(paste(lane, '#', sep = ''), '', barcodes_lane[['x']])
    # paste together the output path
    output_path <- paste(credible_barcodes_loc, '/', lane, '.txt', sep = '')
    # write the barcodes file
    write.table(barcodes_lane, output_path, row.names = F, col.names = F, quote = F)
  }
}
