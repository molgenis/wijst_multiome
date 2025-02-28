#!/bin/bash

#################################################################################
#Script Name	  : mo_split_sample_and_celltype_each.sh                        #
#Description	  : get expression and accessibility for each donor separately  #
#Author       	: Roy Oelen                                                     #
#################################################################################

# these are the cell type
CELL_TYPES=('b' 'cd4t' 'cd8t' 'dc' 'monocyte' 'nk')

# check each celltype
for celltype in ${CELL_TYPES[*]}; do
    # using the script
    ~/start_Rscript.sh /groups/umcg-franke-scrna/tmp04/users/umcg-roelen/singularity/rstudio-server/simulated_home/mo_split_sample_and_celltype.R \
        --cell_type ${celltype} \
        --cell_type_column celltype_imputed_lowerres \
        --seurat_object_path /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_multimodal_${celltype}_1_80_20240521.rds \
        --cre_pairs_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/confinements/region_to_peak_variant_overlaps.tsv.gz \
        --seurat_assignment_column sample_final,lane,condition_final \
        --output_folder /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/cre_eqtl/eqtl_caqtl_overlap/matrices/

done