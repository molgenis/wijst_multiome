"""
This Match GWAS variants to eQTL if within 1MB

authors: Dan Kaptijn, Roy Oelen

example usage:
python 1_gwas_eqtl_snps_match.py \
    --gwas_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/GRN-reconstruction/GWAS/gwas_filtered_updated_08032024.tsv.gz \
    --variant_feature_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_14/Finemapping/finemapped_monocyte_enrichment_input.tsv.gz \
    --out_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/gwas_enrichment/eqtl/sceqtlgen/finemapped/finemapped_monocyte_enrichment_snp_gwas_qtl_for_ld_output.tsv.gz


"""

#############
# libraries #
#############

import pandas as pd
import numpy as np
import argparse


#############
# main code #
#############

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gwas_loc', type = str, help = 'location of gwas file (string)')
parser.add_argument('-v', '--variant_feature_loc', type = str, help = 'location of variant-gene list file (string)')
parser.add_argument('-o', '--out_loc', type = str, help = 'location to write the output')
args = parser.parse_args()


#Location to save file of variants to test for LD
output = '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/gwas_enrichment/eqtl/sceqtlgen/finemapped/finemapped_monocyte_enrichment_output.tsv.gz'
output = args.out_loc
# location of the variant-feature file to check
var_feature_loc = '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_14/Finemapping/finemapped_monocyte_enrichment_input.tsv.gz'
var_feature_loc = args.variant_feature_loc
# location of the GWAS file to check
gwas_loc = '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/GRN-reconstruction/GWAS/gwas_filtered_updated_08032024.tsv.gz'
gwas_loc = args.gwas_loc

#Input file with list of eQTL SNPs to test
### Following loads only the first column of the file which contains: "snp;gene_pair"
df = pd.read_csv(var_feature_loc, sep = '\t',usecols = [0])
#Input file with list of GWAS SNPs to test
gwas = pd.read_csv(gwas_loc, sep = '\t')

### Select the SNP from the column which is in format: chr:pos:alt:ref
snps_to_test = [i.split(";")[0] for i in df.coQTL]
del df

snps_to_test = list(set(snps_to_test))
gwas_snps = [i for i in gwas.chr_pos]
gwas_snps = np.unique(gwas_snps)

gwas_df=pd.DataFrame()
gwas_df['chr']=[i.split('_')[0] for i in gwas_snps]
gwas_df['pos']=[i.split('_')[1] for i in gwas_snps]

snps_df=pd.DataFrame()
snps_df["chr"]=[i.split(":")[0] for i in snps_to_test]
snps_df["pos"]=[int(i.split(":")[1]) for i in snps_to_test]
snps_df["snp_id"]=[i for i in snps_to_test]

results_dict = {}
for chr in range(1,23):
  oneMB = 1000000
  snps_chr = snps_df[snps_df.chr == str(chr)].sort_values(by='pos')
  gwas_pos = sorted([int(i) for i in gwas_df[gwas_df.chr == str(chr)].pos])
  for pos,id in zip(snps_chr.pos, snps_chr.snp_id):
    if id not in results_dict:
      results_dict[id] = ';'.join([str(i) for i in gwas_pos if i>= pos-oneMB and i<= pos+oneMB])


results_df = pd.DataFrame()
for chr_pos in results_dict:
  results_df = results_df._append({'snp': chr_pos,'gwas_loc':results_dict[chr_pos]},ignore_index=True)

# check if the output file ends in gz
if args.out_loc.endswith('.gz'):
    results_df.to_csv(output,sep='\t',index=None, compression = 'gzip')
else:
    results_df.to_csv(output,sep='\t',index=None)


