"""
Match eQTL variants to eQTL if within 2MB

authors: Dan Kaptijn, Roy Oelen

example usage:
python 1_eqtl_eqtl_snps_match.py \
    --variant_feature_loc /groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/Meta_14/Finemapping/finemapped_monocyte_enrichment_input.tsv.gz \
    --out_loc /groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/gwas_enrichment/eqtl/sceqtlgen/finemapped/finemapped_monocyte_snp_snp_for_ld_output.tsv.gz


"""

#############
# libraries #
#############

import pandas as pd
import numpy as np
import argparse
from tqdm import tqdm

#############
# main code #
#############

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-v', '--variant_feature_loc', type = str, help = 'location of variant-gene list file (string)')
parser.add_argument('-o', '--out_loc', type = str, help = 'location to write the output')
args = parser.parse_args()

#Location to save file of variants to test for LD
output = args.out_loc

#Input file with list of eQTL SNPs to test
variant_feature_loc = args.variant_feature_loc
### Following loads only the first column of the file which contains: "snp;gene_pair"
df = pd.read_csv(variant_feature_loc,sep='\t',usecols=[0])

### Select the SNP from the column which is in format: chr:pos:alt:ref
snps_to_test = list(set([i.split(";")[0] for i in df.coQTL]))
del df

snps_df=pd.DataFrame()
snps_df["chr"]=[i.split(":")[0] for i in snps_to_test]
snps_df["pos"]=[int(i.split(":")[1]) for i in snps_to_test]
snps_df["snp_id"]=[i for i in snps_to_test]

results_dict = {}
for chr in tqdm(range(1,23)):
  twoMB = 2000000
  snps_chr = snps_df[snps_df.chr == str(chr)].sort_values(by='pos')
  snps_pos = sorted([int(i) for i in snps_chr.pos])
  for pos,id in zip(snps_chr.pos, snps_chr.snp_id):
    if id not in results_dict:
      results_dict[id] = ';'.join([str(i) for i in snps_pos if i>= pos-twoMB and i<= pos+twoMB])


results_df = pd.DataFrame()
for chr_pos in results_dict:
  results_df = results_df._append({'snp': chr_pos,'gwas_loc':results_dict[chr_pos]},ignore_index=True)


# check if the output file ends in gz
if args.out_loc.endswith('.gz'):
    results_df.to_csv(output,sep='\t',index=None, compression = 'gzip')
else:
    results_df.to_csv(output,sep='\t',index=None)