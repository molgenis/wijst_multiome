#!/bin/bash

###################################################################
#Script Name	  : mo_prepare_onek1k_tfieqtl.sh
#Description	  : take the OneK1K data and convert to format required for TFa-i-eQTL analysis
#Args           : 
#Author       	: Roy Oelen
###################################################################


# convert object into chunks
Rscript mo_create_hybrid_cre_inputs_rna.R \
    --in /groups/umcg-franke-scrna/tmp04/external_datasets/onek1k/seurat_objects/onek1k_major_cts_snumber.rds \
    --out /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/cre_detection/limix_sc/input/tf_interaction/onek1k/L1/all/ \
    --confinement /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eQTA/featureVariantFile.w150k.filtered0.0001_cts.txt \
    --gene_annotation_file /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/eqtl/annotations/cellranger_arc_gene_annotations.tsv.gz \
    --donor_annotation_column sample_final

# convert genotype data from bgen to bed
/groups/umcg-franke-scrna/tmp04/software/plink2_amd/plink2 \
  --bgen /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered.bgen ref-last \
  --make-bed \
  --out /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered

# backup original fam file
mv /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered.fam /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered.fam.original

# add 'S' prefix to sample IDs
R --slave -e '
fam <- read.table("/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered.fam.original", header = T, sep = "\t");
fam[["X1"]] <- paste0("S", fam[["X1"]]);
write.table(fam, "/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered.fam", row.names = F, col.names = F, sep = "\t", quote = F);
mdfiver::create_sha256_for_file("/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered.fam");
'

# split genotype data into chromosome-specific files
for chr in {1..22} X Y; do
  /groups/umcg-franke-scrna/tmp04/software/plink2_amd/plink2 --bfile /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered --chr $chr --make-bed --out /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/genotype_input/EUR_imputed_hg38_varFiltered$chr
done
