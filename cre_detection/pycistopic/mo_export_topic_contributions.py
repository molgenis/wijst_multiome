"""

mo_export_topic_contributions.py

This script exports the topic contributions of each region or cell to a topic

authors: Roy Oelen

"""

#############
# libraries #
#############

import pycisTopic
import pandas as pd
import pickle
import hashlib


##################
# reading object #
##################

# location the object was stored
pycistopic_object_wtopics_loc = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/objects/merged_major_and_minor_celltypes_160topics.pkl'

# read the object
with open(pycistopic_object_wtopics_loc, 'rb') as f:
   cistopic_obj = pickle.load(f)


###########################
# exporting contributions #
###########################

# export the topic contribution matrix for the cells
topic_contribution_cells_loc = ''.join(['/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/models/all_nuclei_and_regions_major_minor/', 'mo_cell_topic_contributions_160_topics.tsv.gz'])
cistopic_obj.selected_model.cell_topic.to_csv(topic_contribution_cells_loc, sep = '\t', header = True, index = True, compression = 'gzip')
# and for the regions
region_contribution_cells_loc = ''.join(['/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/scenicplus_workdir/pycistopic/models/all_nuclei_and_regions_major_minor/', 'mo_region_topic_contributions_160_topics.tsv.gz'])
cistopic_obj.selected_model.region_topic.to_csv(region_contribution_cells_loc, sep = '\t', header = True, index = True, compression = 'gzip')


#########################################
# create md5 checksum for created files #
#########################################

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


# create md5 checksum files for both input files
create_md5_file(topic_contribution_cells_loc)
create_md5_file(region_contribution_cells_loc)
