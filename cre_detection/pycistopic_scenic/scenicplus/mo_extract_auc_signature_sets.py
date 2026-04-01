"""
mo_extract_auc_signature_sets.py

This script is used to extract the gene and region signatures of the eRegulons from the output of SCENIC+ and save them in a tab-separated file because SCENIC+ for some reason doesn't save these anymore
authors: Roy Oelen

"""

###########
# imports #
###########

# for path operations
import os
import glob
from pathlib import Path
# for reading the output of SCENIC+
import scanpy as sc
import anndata
import mudata
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
# for saving results
import joblib
# for visualizing results
import pandas as pd
# for making a checksum
import hashlib
# convert to matrix
from scipy.sparse import csr_array
# write zipped mtx file
import gzip
# the libraries from the original code excerpt
from pycisTopic.diff_features import *
from pycisTopic.signature_enrichment import *
from pyscenic.binarization import binarize
# to make lists of signatures
from itertools import repeat

#############
# functions #
#############

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


def get_eRegulons_as_signatures(scplus_obj,
                                eRegulon_metadata_key: str = 'eRegulon_metadata',
                                key_added: str = 'eRegulon_signatures'):
    """
    Format eRegulons for scoring
    taken from https://scenicplus.readthedocs.io/en/stable/_modules/scenicplus/eregulon_enrichment.html

    Parameters
    ----------
    scplus_obj: `class::SCENICPLUS`
        A SCENICPLUS object with eRegulons metadata computed.
    eRegulon_metadata_key: str, optional
        Key where the eRegulon metadata is stored (in `scplus_obj.uns`)
    key_added: str, optional
        Key where formated signatures will be stored (in `scplus_obj.uns`)
    """
    region_signatures = {x: list(set(scplus_obj.uns[eRegulon_metadata_key][scplus_obj.uns[eRegulon_metadata_key].Region_signature_name == x]['Region'])) for x in list(
        set(scplus_obj.uns[eRegulon_metadata_key].Region_signature_name))}
    gene_signatures = {x: list(set(scplus_obj.uns[eRegulon_metadata_key][scplus_obj.uns[eRegulon_metadata_key].Gene_signature_name == x]['Gene'])) for x in list(
        set(scplus_obj.uns[eRegulon_metadata_key].Gene_signature_name))}
    if not key_added in scplus_obj.uns.keys():
        scplus_obj.uns[key_added] = {}
    scplus_obj.uns[key_added]['Gene_based'] = gene_signatures
    scplus_obj.uns[key_added]['Region_based'] = region_signatures


###################
# paths of inputs #
###################

# location of the mu that has the eregulon data
eregulon_mu_loc = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/scplus_pipeline_merged_major_and_minor_celltypes/output/scplusmdata.h5mu'

###################
# load the inputs #
###################

# read the scenic output metadata
scplus_mdata = mudata.read(eregulon_mu_loc)

# add direct eregulon info
get_eRegulons_as_signatures(scplus_mdata, eRegulon_metadata_key = 'direct_e_regulon_metadata', key_added = 'direct_e_regulon_signatures')
# and indirect
get_eRegulons_as_signatures(scplus_mdata, eRegulon_metadata_key = 'extended_e_regulon_metadata', key_added = 'extended_e_regulon_signatures')

# make into lists of signatures
eregs_to_source_to_modality_to_signatures = []
# read the sources
for source in ['direct_e_regulon_signatures', 'extended_e_regulon_signatures']:
    # read the modalities
    for modality in ['Gene_based', 'Region_based']:
        # extract the signatures
        signatures = scplus_mdata.uns[source][modality]
        # check signatures
        for (signature_name, signature_genes_or_regions) in signatures.items():
            if len(signature_genes_or_regions) > 0:
                # make into a pandas dataframe
                signature_df = pd.DataFrame(
                    {
                        'source': list(repeat(source, len(signature_genes_or_regions))),
                        'modality': list(repeat(modality, len(signature_genes_or_regions))),
                        'signature_name': list(repeat(signature_name, len(signature_genes_or_regions))),
                        'gene_or_region': list(signature_genes_or_regions)
                    }
                )
                eregs_to_source_to_modality_to_signatures.append(signature_df)

# concatenate all signatures into one dataframe
eregs_to_source_to_modality_to_signatures_df = pd.concat(eregs_to_source_to_modality_to_signatures, ignore_index=True)
# write the result
output_signature_df_loc = os.path.join(os.path.dirname(eregulon_mu_loc), 'eRegulon_signatures.tsv.gz')
eregs_to_source_to_modality_to_signatures_df.to_csv(output_signature_df_loc, sep='\t', index=False)
# and make an md5 of the output file
create_md5_file(output_signature_df_loc)
