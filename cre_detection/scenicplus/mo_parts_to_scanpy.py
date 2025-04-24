"""
python mo_parts_to_scanpy.py

This script is for creating a scanpy object

authors: Roy Oelen

example usage:

python mo_parts_to_scanpy.py \
    --matrix_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/all/RNA/counts/matrix.mtx \
    --barcodes_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/all/RNA/counts/barcodes.tsv.gz \
    --features_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/all/RNA/counts/features.tsv.gz \
    --output_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/reconstructed_objects/merged_major_and_minor_celltypes.h5ad \
    --metadata_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scanpy_objects/deconstructed_objects/all/metadata.tsv.gz \
    --barcode_include_location /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/merged_major_and_minor_celltypes/barcodes.tsv.gz \
    --skip_qc

"""

###########
# imports #
###########

import numpy as np
import pandas as pd
import argparse
import pathlib
import scanpy as sc
import hashlib


##########################################
# create md5 checksum for created object #
##########################################

def create_md5_file(input_file):
    """
    Creates an MD5 hash of the specified file and writes it to a new file with the same name but .md5 added to the extension.

    Args:
        input_file (str): The path to the input file for which the MD5 hash should be created.

    Returns:
        int: Returns 0 on success, 1 on failure.

    Raises:
        FileNotFoundError: Thrown if the input file does not exist.
        IOError: Thrown if there is an error reading the input file (like permission denied) or writing the output md5 file.
    """
    try:
        # get an md5 of the file
        digest = None
        with open(input_file, "rb") as f:
            # try Python 3.11+ method if it is available
            if callable(getattr(hashlib, 'file_digest', None)):
                # digest with one command
                digest = hashlib.file_digest(f, 'md5')
            # or the older 3.8+ method if we don't have the newer method
            else:
                # initialize digest
                digest = hashlib.md5()
                # read file in chunks
                while chunk := f.read(8192):
                    # update digestion
                    digest.update(chunk)       
        # get the output path of the md5
        output_md5_loc = ''.join([input_file, '.md5'])
        # and write that
        with open(output_md5_loc, "w") as m:
            m.write(digest.hexdigest())
        # upon success, return 0
        return 0
    except Exception as e:
        print(f"Exception occured upon md5 file creation: {e}")
        return 1


###################
# parse arguments #
###################

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--matrix_location', type = str, help = 'location of the matrix.mtx file (string)')
parser.add_argument('-b', '--barcodes_location', type = str, help = 'location of the barcodes.tsv.gz file (string)')
parser.add_argument('-f', '--features_location', type = str, help = 'location of the features.tsv.gz file (string)')
parser.add_argument('-o', '--output_location', type = str, help = 'location of the resulting scanpy object (string)')
parser.add_argument('-d', '--metadata_location', type = str, help = 'location of metadata to add to the object, first column must be the index (string)', default = None)
parser.add_argument('-i', '--barcode_include_location', type = str, help = 'text file containing barcodes of cells to include, will only keep these cells (string)', default = None)
parser.add_argument('-q', '--skip_qc', action='store_true', help = 'perform quality control on the the data (bool)', default = True)
args = parser.parse_args()


##############
# processing #
##############

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
# if we have a list of barcodes to include, read those
if args.barcode_include_location is not None:
    bc_inclusion_list = pd.read_csv(args.barcode_include_location, header=None)[0].tolist()
    # and filter on those barcodes
    object_raw = object_raw[bc_inclusion_list].copy()
# backup the raw expression
object_raw.layers['counts'] = object_raw.X.copy()


##############
# perform QC #
##############

# QC steps
object_raw.var["mt"] = object_raw.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    object_raw, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
# perform the QC if requested
if args.do_qc is not False:
    # removing cells with few reads
    object_raw = object_raw[object_raw.obs.n_genes_by_counts < 2500, :].copy()
    # and high MT content
    object_raw = object_raw[object_raw.obs.pct_counts_mt < 5, :].copy()


#################
# normalization #
#################

# normalization
sc.pp.normalize_total(object_raw, target_sum=1e4)
sc.pp.log1p(object_raw)
sc.pp.highly_variable_genes(object_raw, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(object_raw, max_value=10)


#####################################
# dimensional recution and clusters #
#####################################

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


#################
# create output #
#################

# write the result
object_raw.write_h5ad(args.output_location)
# and make an md5
create_md5_file(args.output_location)