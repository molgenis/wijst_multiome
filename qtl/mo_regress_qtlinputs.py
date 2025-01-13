#!/usr/bin/env python3

"""
This script is used to regress out principal components from gene expression of open chromatin data
authors: Roy Oelen

example usage:

python mo_regress_qtlinputs.py \
    --qtl_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/NK.qtlInput.txt.gz \
    --pc_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/NK.qtlInput.Pcs.txt.gz \
    --output_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/caqtl/sc-eqtlgen/input/L1/UT/NK.qtlInput.PcCorrectedResiduals.txt.gz

"""

###############
# libraries   #
###############

# for checking file locations
import os
import glob
import os.path
# use warnings
import warnings
# use regex
import re
# for pandas
import pandas as pd
# do regression
from sklearn.linear_model import LinearRegression
# for multicore processing
import concurrent.futures
# argumenting
import argparse
# md5 checksum creation
import hashlib


###############
# functions   #
###############

# method to do a single gene
def regress_out_pc(gene):
    y = input_data.loc[gene].values
    X = rna_pcs_matched.values
    model.fit(X, y)
    residuals = y - model.predict(X)
    return gene, residuals

# method to create md5
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

###############
# optparse    #
###############

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-q', '--qtl_loc', type = str, help = 'location of the QTL input file (string)')
parser.add_argument('-p', '--pc_loc', type = str, help = 'location of the PCs file to regress out (string)')
parser.add_argument('-o', '--output_loc', type = str, help = 'location of the regressed QTL output (string)')
args = parser.parse_args()


###############
# main        #
###############

# now for some atac data
input_data_loc = args.qtl_loc
# read the atac data
input_data = pd.read_csv(input_data_loc, header = 0, sep = '\t', index_col = 0)
# location of the PCs
input_pcs_loc = args.pc_loc
# read the PC data
input_pcs = pd.read_csv(input_pcs_loc, header = 0, sep = '\t', index_col = 0)
# get the samples present in both
samples_both = list(set(input_data.columns).intersection(input_pcs.index))
# subset to the samples present in both
input_data_matched = input_data.loc[:, samples_both]
input_pcs_matched = input_pcs.loc[samples_both, ]
# create the model
model = LinearRegression()
# Create a new dataframe to store the residuals
input_data_matched_residuals = pd.DataFrame(index=input_data_matched.index, columns=input_data_matched.columns)
# method to do a single gene
def regress_out_pc(gene):
    y = input_data_matched.loc[gene].values
    X = input_pcs_matched.values
    model.fit(X, y)
    residuals = y - model.predict(X)
    return gene, residuals

# do this concurrently
with concurrent.futures.ThreadPoolExecutor() as executor:
    futures = {executor.submit(regress_out_pc, gene): gene for gene in input_data_matched.index}
    # submit each job
    for future in concurrent.futures.as_completed(futures):
        gene, residuals = future.result()
        input_data_matched_residuals.loc[gene] = residuals

# write the result
residuals_loc = args.output_loc
input_data_matched_residuals.to_csv(residuals_loc, sep = '\t', header = True, index = True, compression = 'gzip')

# write this to a file
if residuals_loc.endswith('.gz'):
    input_data_matched_residuals.to_csv(residuals_loc, sep = '\t', header = True, index = True, compression = 'gzip')
else:
    input_data_matched_residuals.to_csv(residuals_loc, sep = '\t', header = True, index = True)

# and make an md5
create_md5_file(residuals_loc)
