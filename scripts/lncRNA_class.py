import pandas as pd
import sys
import numpy as np
from pybedtools import BedTool

def clean(x):
    split = x.str.rsplit('.', 2, expand=True)
    split = split[0]
    return (split)

tmap = pd.read_csv(sys.argv[1], sep='\t')
gmap_align = pd.read_csv(sys.argv[2], sep='\t', header=None)
# take only lncRNA classes from tmap file
tmap['group'] = np.nan
tmap['group'][tmap['class_code'].isin(['x'])] = "exon antisense"
tmap['group'][tmap['class_code'].isin(['i'])] = "intron"
tmap['group'][tmap['class_code'].isin(['u'])] = "intergenic"
tmap.dropna(subset = ["group"], inplace=True)
tmap = tmap[tmap['qry_gene_id'].str.contains("path1")]
print(tmap['group'].value_counts())
# filter duplicates
query = clean(tmap['qry_gene_id'])
gmap_align['name'] = clean(gmap_align[3])
query = query.drop_duplicates()
# take lncRNA from alignment file
gmap_align = gmap_align[gmap_align['name'].isin(query)]
gmap_align = gmap_align[gmap_align[7].str.contains("gene")]
gmap_align.to_csv(sys.argv[3], sep='\t', index=False, header=None)
# lncRNA into loci
genes = BedTool(sys.argv[3])
loci = genes.merge(s=True, c=[4,5,6], o=['collapse','mean','distinct'])
loci_df = pd.read_table(loci.fn, sep='\t', names=[0,1,2,3,4,5])
loci_df = loci_df[(loci_df[2] - loci_df[1])>200]
loci_df['numb'] = list(range(0, len(loci_df[3])))
loci_df['loc'] = 'LOC'
loci_df['loci'] = loci_df['loc'] + '_' + loci_df['numb'].astype(str)
loci_df = loci_df[[0,1,2,'loci',4,5,3]]
loci_df.to_csv(sys.argv[4], sep='\t', index=False, header=None)

