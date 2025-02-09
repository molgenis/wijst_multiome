"""
This script is used to combine multiple mtx files into an npz file

authors: Roy Oelen

"""

#############
# libraries #
#############

import os
from scipy.io import mmread, mmwrite
from scipy.sparse import vstack, save_npz, load_npz
import re
import argparse
import pickle
import gzip


#############
# functions #
#############


def combine_mtx_files(mtx_file_list):
    """
    Combine multiple .mtx files into a single sparse matrix and save as a .npz file.

    Args:
        mtx_file_list (list of str): List of paths to the .mtx files to be combined.

    Returns:
        npz matrix
    """
    
    # initialize the matrix
    combined_matrix = None
    
    # read each file
    for file in mtx_file_list:
        print(' '.join(['reading file', str(file)]))
        # read and convert to csc
        matrix = mmread(file).tocsc()
        #matrix = sc.read_mtx(file).X.tocsc()
        # set to be the matrix if no matrix was set yet
        if combined_matrix is None:
            combined_matrix = matrix
        # otherwise, vstack onto the existing matrix
        else:
            combined_matrix = vstack([combined_matrix, matrix])
    # return the result
    return combined_matrix



##################
# argument setup #
##################

args_mtx_folder = '/'
args_mtx_regex = 'matrix_chunk_\d+\.mtx.gz'
args_npz_out = 'matrix_merged.mtx'

#############
# Main code #
#############

# List all files in the directory
files = os.listdir(args_mtx_folder)
# compile a pattern
pattern = re.compile(f'{args_mtx_regex}')
# filter by pattern
filtered_files = [f for f in files if os.path.isfile(os.path.join(args_mtx_folder, f)) and pattern.match(f)]
# order the files
filtered_files.sort()
# add the path
filtered_files = [os.path.join(args_mtx_folder, f) for f in filtered_files]
# do the conversion
combined_matrix = combine_mtx_files(filtered_files)

# this is where we will output:
npz_output_file = args_npz_out

# write this to a file
if npz_output_file.endswith('.npz'):
    # to npz
    save_npz(npz_output_file, combined_matrix)
elif npz_output_file.endswith('.pickle'):
    # or as pickle
    with open(npz_output_file, 'wb') as handle:
        pickle.dump(combined_matrix, handle, protocol=pickle.HIGHEST_PROTOCOL)
elif npz_output_file.endswith('.mtx.gz'):
    # or mtx.gz
    with gzip.open(npz_output_file, 'wb') as f:
        mmwrite(f, combined_matrix)
else:
    # default to npz
    print(''.join(['not recognizing output format, saving as npz']))
    save_npz(npz_output_file, combined_matrix)