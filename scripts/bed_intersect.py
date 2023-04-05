import pandas as pd
import sys
import numpy as np
from pybedtools import BedTool
pd.set_option('display.max_columns', None)

def clean(x):
    split = x.str.split('.', 4, expand=True)
    split = split[0]
    return (split)

tmap = pd.read_csv(sys.argv[1], sep='\t')
gmap_align_prime = pd.read_csv(sys.argv[2], sep='\t', header=None)

# take only lncRNA classes from tmap file
tmap['group'] = np.nan
tmap['group'][tmap['class_code'].isin(['='])] = 'true lncrna'
tmap.dropna(subset = ["group"], inplace=True)
tmap = tmap[tmap['qry_gene_id'].str.contains("path1")]
print(tmap['group'].value_counts())
# filter duplicates
query = clean(tmap['qry_gene_id'])
gmap_align_prime['name'] = clean(gmap_align_prime[3])
query = query.drop_duplicates()
# take lncRNA from alignment file
if str(sys.argv[3]) == 'on':
   print('TE finder ON')
   gmap_align_prime = gmap_align_prime[gmap_align_prime['name'].isin(query)]
   gmap_align = gmap_align_prime[gmap_align_prime[7].str.contains("gene")]
   gmap_align.to_csv(sys.argv[4], sep='\t', index=False, header=None)
# BedTool intersect
   genes = BedTool(sys.argv[4])
   loci = genes.intersect(wa=True, wb=True, b=sys.argv[5])
   loci_df = pd.read_table(loci.fn, sep='\t', names=['chr','lncRNA_start','lncRNA_end','lncRNA',4,'lncRNA_strand',6,7,8,9,10,'chr_TE','TE_start','TE_end'])
   loci_df = loci_df[['chr','lncRNA_start','lncRNA_end','lncRNA','lncRNA_strand','chr_TE','TE_start','TE_end']]
   loci_df = loci_df.drop_duplicates(subset=['lncRNA'])
   loci_name = loci_df['lncRNA'].str.split('.',4, expand=True)
   true_lncrna = gmap_align_prime[~gmap_align_prime['name'].isin(loci_name[0])]
   loci_df.to_csv(sys.argv[6], sep='\t', index=False, header=None)
   true_lncrna.to_csv(sys.argv[7], sep='\t', index=False, header=None)
else:
   gmap_align_prime = gmap_align_prime[gmap_align_prime['name'].isin(query)]
   gmap_align_prime.to_csv(sys.argv[4], sep='\t', index=False, header=None)
