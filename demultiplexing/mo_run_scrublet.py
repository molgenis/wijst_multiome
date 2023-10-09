"""
This script is for doing transcription-based doublet detection

authors: Roy Oelen
"""

# imports
import scrublet as scr
import scipy.io
import numpy as np
import pandas as pd
import argparse
import pathlib
import scanpy as sc

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-l', '--lane_location', type = str, help = 'location of the 10x lane (string)')
parser.add_argument('-o', '--output_location', type = str, help = 'location of the scrublet output file (string)')
parser.add_argument('-r', '--expected_doublet_rate', type = float, help = 'expected doublet rate (float)')
parser.add_argument('-m', '--min_gene_variability_pctl', type = float, help = 'minimal gene variability (float)')
args = parser.parse_args()

# location of the lane
lane_location = args.lane_location
# grab the lane
path = pathlib.PurePath(lane_location)
lane = path.name
# specific lane location
lane_append = '/outs/filtered_feature_bc_matrix/'
lane_append = ''
# we're expecting around 8% doublets, given the number of cells
#expected_doublet_rate=0.08
expected_doublet_rate=args.expected_doublet_rate

# state where to save to
scrub_save_loc = args.output_location

# get the full paths
matrix_loc = "".join([lane_location, lane_append, 'cellbender_remove_background_output_filtered.h5'])
barcodes_loc = "".join([lane_location, lane_append, 'cellbender_remove_background_output_cell_barcodes.csv'])
# load the files
counts_matrix = sc.read_10x_h5(matrix_loc)
#genes = np.array(scr.load_genes(features_loc, delimiter='\t', column=1))
barcodes = pd.read_csv(barcodes_loc, sep= '\t', header=None)
# perform scrublet
scrub = scr.Scrublet(counts_matrix.X, expected_doublet_rate=expected_doublet_rate)
# grabbing scores and assignments
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=args.min_gene_variability_pctl, n_prin_comps=30)
# add everything together in one frame
assignment = barcodes
assignment['doublet'] = predicted_doublets
assignment['doublet_score'] = doublet_scores
assignment['barcode'] = assignment[0].str.slice(0, 16)
assignment['lane_barcode'] = assignment['barcode']+'_'+lane
# remove original column
assignment = assignment.drop(0, axis=1)

# save to file
assignment.to_csv(scrub_save_loc, index=False, sep='\t')

# and make plot
scrub.plot_histogram();


# python mo_run_scrublet.py \
# -l /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/2023_09_12_cellbender-v0.3.0/default-run/230105_lane1/ \
# -o /groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/demultiplexing/scrublet/scrublet_output/mo_scrublet_230105_lane1.tsv \
# -r 0.1 \
# -m 0.85
