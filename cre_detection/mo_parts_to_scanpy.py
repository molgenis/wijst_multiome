"""
This script is for creating a scanpy object

authors: Roy Oelen

example usage:

python mo_parts_to_scanpy.py \
    --matrix_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/ASDC/RNA/counts/matrix.mtx \
    --barcodes_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/ASDC/RNA/counts/barcodes.tsv.gz \
    --features_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/ASDC/RNA/counts/features.tsv.gz \
    --output_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/reconstructed_objects/ASDC.h5ad \
    --metadata_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/ASDC/metadata.tsv.gz

"""

# imports
import numpy as np
import pandas as pd
import argparse
import pathlib
import scanpy as sc

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--matrix_location', type = str, help = 'location of the matrix.mtx file (string)')
parser.add_argument('-b', '--barcodes_location', type = str, help = 'location of the barcodes.tsv.gz file (string)')
parser.add_argument('-f', '--features_location', type = str, help = 'location of the features.tsv.gz file (string)')
parser.add_argument('-o', '--output_location', type = str, help = 'location of the resulting scanpy object (string)')
parser.add_argument('-d', '--metadata_location', type = str, help = 'location of metadata to add to the object, first column must be the index (string)', default = None)
args = parser.parse_args()

# read count data
object_counts = sc.read_mtx(args.matrix_location)
# read barcodes
object_bc = pd.read_csv(args.barcodes_location, header=None)
# read features
object_features = pd.read_csv(args.features_location, header=None)
# transpose to scanpy format
object_raw = sc.AnnData(object_counts.T)
# add barcodes and genes to obs and vars
object_raw.obs['cell_id'] = object_bc[0].tolist()
object_raw.var['gene_name'] = object_features[0].tolist()
# set indices for the obs and vars
object_raw.obs.index = object_raw.obs['cell_id']
object_raw.var.index = object_raw.var['gene_name']
# if we have metadata, read that as well
if args.metadata_location is not None:
    object_metadata = pd.read_csv(args.metadata_location, header=0, sep = '\t', index_col = 0, low_memory = False)
    # and join onto the scanpy object
    object_raw.obs = object_raw.obs.merge(object_metadata, how = 'left', left_index = True, right_index = True)
# backup the raw expression
object_raw.raw = object_raw
# QC steps
object_raw.var["mt"] = object_raw.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    object_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
object_raw = object_raw[object_raw.obs.n_genes_by_counts < 2500, :]
object_raw = object_raw[object_raw.obs.pct_counts_mt < 5, :].copy()
# normalization
sc.pp.normalize_total(object_raw, target_sum=1e4)
sc.pp.log1p(object_raw)
sc.pp.highly_variable_genes(object_raw, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(object_raw, max_value=10)
# dimensional reduction
sc.tl.pca(object_raw, svd_solver="arpack")
sc.pp.neighbors(object_raw, n_neighbors=10, n_pcs=40)
sc.tl.umap(object_raw)
# clustering
sc.tl.leiden(
    object_raw,
    resolution=0.9,
    random_state=0,
    n_iterations=2,
    directed=False,
)
# marker genes
sc.tl.rank_genes_groups(object_raw, "leiden", method="t-test")
# write the result
object_raw.write_h5ad(args.output_location)