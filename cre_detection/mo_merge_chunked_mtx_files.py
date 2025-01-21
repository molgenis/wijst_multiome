"""
This script is used to combine multiple mtx files into an npz file

authors: Roy Oelen

example usage:

python mo_merge_chunked_mtx_files.py \
    --mtx_folder /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/merged_major_celltypes/ \
    --mtx_regex 'matrix_chunk_\\d+.mtx' \
    --npz_out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/deconstruced_atac_objects/merged_major_celltypes/matrix_merged.npz

"""

#############
# libraries #
#############

import os
from scipy.io import mmread, mmwrite
from scipy.sparse import vstack, save_npz, load_npz
import re
import argparse
import scanpy as sc

#############
# functions #
#############


def combine_mtx_files(mtx_file_list, npz_output_file):
    """
    Combine multiple .mtx files into a single sparse matrix and save as a .npz file.

    Args:
        mtx_file_list (list of str): List of paths to the .mtx files to be combined.
        npz_output_file (str): Path to the output .npz file.

    Returns:
        0 upon success
    """
    
    # initialize the matrix
    combined_matrix = None
    
    # read each file
    for file in mtx_file_list:
        # read and convert to csc
        #matrix = mmread(file).X.tocsc()
        matrix = sc.read_mtx(file).X.tocsc()
        # set to be the matrix if no matrix was set yet
        if combined_matrix is None:
            combined_matrix = matrix
        # otherwise, vstack onto the existing matrix
        else:
            combined_matrix = vstack([combined_matrix, matrix])
    # write the result
    save_npz(npz_output_file, combined_matrix)
    return 0



####################
# argument parsing #
####################

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--mtx_folder', type = str, help = 'location of the mtx files (string)')
parser.add_argument('-r', '--mtx_regex', type = str, help = 'regex pattern for the mtx files (string)')
parser.add_argument('-o', '--npz_out', type = str, help = 'full output location of the npz output file (string)')
args = parser.parse_args()


#############
# Main code #
#############

# List all files in the directory
files = os.listdir(args.mtx_folder)
# compile a pattern
pattern = re.compile(f'{args.mtx_regex}')
# filter by pattern
filtered_files = [f for f in files if os.path.isfile(os.path.join(args.mtx_folder, f)) and pattern.match(f)]
# do the conversion
combine_mtx_files(files, args.npz_out)