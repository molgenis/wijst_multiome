"""
mo_pycistopic_find_variable_features.py

This script is used to normalize the accessibility

authors: Roy Oelen, Martijn van der Werf

"""

###########
# imports #
###########

# object
from pycisTopic.cistopic_class import CistopicObject
import pickle
# imputation
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
# for storing the result
import joblib
# for making the md5
import hashlib


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


##############################
# read the pycistopic object #
##############################

# load the object from disk
pycistopic_object_wimpnorm_loc = '/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/objects/merged_major_and_minor_celltypes_120topics_impnorm.joblib'
normalized_imputed_acc_obj = joblib.load(pycistopic_object_wimpnorm_loc)


########################
# get variable regions #
########################

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)

##########################
# save variable regions #
#########################

# location to store the object
var_features_loc = '/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/dar_detection/merged_major_and_minor_celltypes_120topics/wilcoxon/merged_major_and_minor_celltypes_120topics_varfeatures.joblib'

# save the object
joblib.dump(variable_regions, var_features_loc)

# make a checksum
create_md5_file(var_features_loc)
