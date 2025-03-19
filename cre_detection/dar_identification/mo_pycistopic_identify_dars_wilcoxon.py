"""
mo_pycistopic_identify_dars_wilcoxon.py

This script is used to identify DARs using the Wilcoxon method

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


##########################################
# read the non-imputed pycistopic object #
##########################################

# location to store the object
pycistopic_object_loc = '/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/objects/merged_major_and_minor_celltypes_120topics.pkl'
# use symlinks due to path size limitations
pycistopic_object_loc = './120'

# save the object
with open(pycistopic_object_loc, 'rb') as f:
   cistopic_obj = pickle.load(f)


######################################
# read the imputed pycistopic object #
######################################

# load the object from disk
pycistopic_object_wimputations_jl_loc = '//scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/objects/merged_major_and_minor_celltypes_120topics_imputed.joblib'
imputed_acc_obj = joblib.load(pycistopic_object_wimputations_jl_loc)


##################################################
# run DAR identification using wilcoxon-rank-sum #
##################################################

# Run DAR analysis
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=4,
    _temp_dir=os.environ["TMPDIR"],
    split_pattern = '_'
)

# location to store the object
markers_loc = '/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/dar_detection/merged_major_and_minor_celltypes_120topics/wilcoxon/merged_major_and_minor_celltypes_120topics_dars.joblib'

# save the object
joblib.dump(markers_dict, markers_loc)

# make a checksum
create_md5_file(markers_loc)
