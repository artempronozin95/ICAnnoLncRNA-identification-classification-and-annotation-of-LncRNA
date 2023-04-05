import pandas as pd
import sys

def coord(x):
    uni = x.drop_duplicates(subset=[1])
    Row_list_query = []
    for index, rows in uni.iterrows():
        my_list = [rows[1], rows[2]]
        Row_list_query.append(my_list)
    return(Row_list_query)

gmap  = sys.argv[1]
known_seq = sys.argv[2]
out = sys.argv[3]

gmap_align = pd.read_csv(gmap, sep='\t', header=None)
target = pd.read_csv(known_seq, sep='\t', header=None)

target_coords = coord(target[target[7].str.contains("transcript")])
name = gmap_align[9].str.split(';',4, expand=True)
name = name[1].str.split('=',2, expand=True)
gmap_align['name'] = name[1]
gene = gmap_align[gmap_align[3].str.contains("path1")]
exon = gmap_align[gmap_align[3].str.contains("exon")]
gr = exon[[0,1,2,'name']].groupby(by='name')
intron_length = []
intron_data = []
for key, item in gr:
    n = 0
    if 0 < len(item['name']) - 2:
        intron = abs(item[2].iloc[n] - item[1].iloc[n + 1]) - 1
        intron_data.append([item[0].iloc[n], item[2].iloc[n], item[1].iloc[n + 1], key + '_intron' + str(n+1),intron])
        n = n + 1
        if intron > 200000:
            large = gene[gene['name'].str.contains(key)]
            target_large = coord(large)
            for start_stop1 in target_coords:
                start1, stop1 = map(int, start_stop1)
                for start_stop in target_large:
                    start, stop = map(int, start_stop)
                    if start1 in range(start, stop):
                        intron_length.append(key)
        else:
           pass
    else:
       pass
df = pd.DataFrame(intron_data)
df.to_csv('data/output/lncRNA_structure/itron_coordin.tsv', sep='\t', index=False, header=None)
gmap_align = gmap_align[~gmap_align['name'].isin(intron_length)]
gmap_align['length'] = gmap_align[2] - gmap_align[1]
gmap_align = gmap_align[gmap_align['length'] < 5000]
gmap_align_big = gmap_align[gmap_align['length'] > 5000]
gmap_align_big.to_csv('data/output/lncRNA_structure/long_transcripts.bed', index=False, sep='\t', header=None)
gmap_align = gmap_align[~gmap_align['name'].isin(gmap_align_big['name'])]
gmap_align = gmap_align[[0,1,2,3,4,5,6,7,8,9]]
gmap_align.to_csv(out, index=False, sep='\t', header=None)

