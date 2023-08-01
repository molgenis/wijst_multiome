#!/bin/bash

############################################################################################################################
# Authors: Roy Oelen
# Name: mo_preprocess_previous_genotypes.sh
# Function: create md5 checksums for each file
############################################################################################################################

# the lanes we want to do
CHROMS=('1' '2' '3' '4' '5' \
'6' '7' '8' '9' '10' \
'11' '12' '13' '14' '15' \
'16' '17' '18' '19' '20' \
'21' '22' \
'X' 'XY' \
)

# check chromosome
for chr in ${CHROMS[*]}
  do
    # subset to the individuals we care about
    /groups/umcg-franke-scrna/tmp01/software/plink2_20230707/plink2 --bfile /groups/umcg-lifelines/prm03/releases/gsa_genotypes/v1/Data/UGLI_QCed_genotypes/chr_${chr} --keep /groups/umcg-franke-scrna/tmp01/projects/multiome/ongoing/metadata/mo_ugli_plinkfilter.tsv --make-bed --out /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_${chr}
    # convert to pgen
    # /groups/umcg-franke-scrna/tmp01/software/plink2_20230707/plink2 --bfile /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_${chr} --make-pgen --out /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_${chr}
done

# also create an annotation of all the plink files
rm -f /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chrom_list.txt
for chr in ${CHROMS[*]}
  do
    echo chr_${chr} >> /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chrom_list.txt
done

# merge all the bgens together
plink -bfile /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed//chr_1 --make-bed --merge-list chrom_list.txt --out /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_all
# convert to pgen
/groups/umcg-franke-scrna/tmp01/software/plink2_20230707/plink2 --bfile /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_all --make-pgen --out /groups/umcg-franke-scrna/tmp01/projects/multiome/processed/genotype/ugli/unimputed/chr_all
