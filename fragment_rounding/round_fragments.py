"""
This script is use to round atac fragments from CellRanger ARC, to represent fragment instead of read counts
authors: Martijn van der Werf

example usage:

python round_fragments.py 20230105
   
"""

import gzip
import sys

command_line_lane = sys.argv[1].strip()

lane = f'/groups/umcg-franke-scrna/tmp02/projects/multiome/processed/joint/alignment/b38/{command_line_lane}/outs/atac_fragments.tsv.gz' 
outfile = f'/groups/umcg-franke-scrna/tmp02/projects/multiome/ongoing/rounded_fragments/{command_line_lane}_rounded_fragments.tsv.gz'
fh = gzip.open(outfile,'wt')
with gzip.open(lane, 'rt') as lane:
    while True:  
        line = lane.readline()
        if line.startswith('#'):
            fh.write(line)
        else:
            line = line.strip().split('\t')
            if int(line[4]) %2:
                line[4] = (int(line[4]) + 1) / 2
            fh.write(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+str(int(line[4]))+"\t"+"\n")
        if not line:
            break
fh.close()
