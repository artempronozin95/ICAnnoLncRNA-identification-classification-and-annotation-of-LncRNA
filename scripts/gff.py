import pandas as pd
import sys
from collections import Counter
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib.ticker import FuncFormatter
from functools import partial
import math
import statistics
from scipy.stats import lognorm
import numpy as np

pd.options.mode.chained_assignment = None  # default='warn'

def to_percent(y, position, n):
    s = str(round(100 * y / n, 3))
    return s + '%'

def coord(x):
    Row_list_query = []
    for index, rows in x.iterrows():
        my_list = [rows[1], rows[2]]
        Row_list_query.append(my_list)
    return(Row_list_query)

def find_same_coords(x):
    x['exon_num'] = x['name'] +"_"+ exon[7].astype(str).str.cat(x['number'].astype(str), sep='_')
    uni3 = x.drop_duplicates(subset=[1])
    target_coord = coord(uni3)
    for start_stop1 in query_coord:
        start1, stop1 = map(int, start_stop1)
        for start_stop in target_coord:
            start, stop = map(int, start_stop)
            if start1 in range(start-100, start + 100) and stop1 in range(stop - 100, stop+100):
                same_query = uni[1] == start1
                same_target = uni3[1] == start
                index = uni3[same_target]
                lines = uni[same_query]
                chrom.append([' '.join('{0}'.format(i) for i in lines[3]),
                                  ' '.join('{0}'.format(i) for i in index['exon_num']), str(round(start1)), str(round(stop1)), str(start),
                                  str(stop), ' '.join('{0}'.format(i) for i in index['number']),
                                  ' '.join('{0}'.format(i) for i in lines['length']), ' '.join('{0}'.format(i) for i in index[5]), ' '.join('{0}'.format(i) for i in lines[5])])

query_bed = sys.argv[1]
reference_bed = sys.argv[2]
gff_tmap = sys.argv[3]

query = pd.read_csv(query_bed, sep='\t', header=None)
target = pd.read_csv(reference_bed, sep='\t', header=None)
tmap = pd.read_csv(gff_tmap, sep='\t', header=None, skiprows=1)
query = query[pd.to_numeric(query[0], errors='coerce').notnull()]
target = target[pd.to_numeric(target[0], errors='coerce').notnull()]

name = query[3].str.split('.',4, expand=True)
query['name'] = name[0]
name_tmap = tmap[3].str.split('.',4, expand=True)
tmap['name'] = name_tmap[0]

tmap = tmap[tmap[2].isin(['x', 'i', 'u'])]
query = query[query['name'].isin(tmap['name'])]


# exon structure of lncRNA
chart = query[query[3].str.contains("exon")]
spl = chart[3].str.split('.',2, expand=True)
size = spl.pivot_table(index = [0], aggfunc ='size').to_dict()
res = Counter(size.values())
dict={}
for k, v in res.items():
    dict[k] = v/len(size)*100
f, ax = plt.subplots(figsize=(15, 5))
ax.bar(list(dict.keys()), dict.values(), color='b')
#ax.set_xlim([0, 20])
plt.xticks(np.sort(range(max(list(dict.keys())) + 1)))
plt.ylabel('Proportion of lncRNAs (%)', size=15)
plt.xlabel('Number of exon', size=15)
plt.savefig("data/output/number_of_exon.png", dpi=500)

#exon size chart
chart['exon_length'] = (chart[2] - chart[1])
step = int(math.ceil(chart['exon_length'].mean() / 100.0)) * 100
f, ax = plt.subplots(figsize=(10, 5))
density = stats.gaussian_kde(chart['exon_length'])
n, x, _ = ax.hist(chart['exon_length'],  density=True, bins = np.logspace(np.log10(min(chart['exon_length'])),np.log10(max(chart['exon_length'])), 50),  color='b')
ax.plot(x, density(x), color='r')
ax.set_xscale('log')
ax.set_ylabel('Probability', size=15)
ax.set_xlabel('Exon size (bp)', size=15)
plt.savefig("data/output/exon_size.png", dpi=500)

# intron size chart
chart['name'] = spl[0]
gr = chart[[1,2,'name']].groupby(by='name')
intron_length = []
name = []
for key, item in gr:
    n=0
    if n <= len(item['name']) - 2:
        intron = abs(item[2].iloc[n] - item[1].iloc[n + 1]) - 1
        n = n + 1
        intron_length.append(intron)
    else:
        intron_length.append(0)

intron_length = [n for n in intron_length if n < 100000]
intron_length = list(filter(lambda x: x != 0, intron_length))
step = int(math.ceil(statistics.mean(intron_length) / 100.0)) * 100
f, ax = plt.subplots(figsize=(10, 5))
density = stats.gaussian_kde(intron_length)
n ,x , _ = ax.hist(intron_length, density=True, bins=np.logspace(np.log10(min(intron_length)),np.log10(max(intron_length)), 50), color='b')
ax.plot(x, density(x), color="r")
ax.set_xscale('log')
ax.set_ylabel('Probability', size=15)
ax.set_xlabel('Intron size (bp)', size=15)
plt.savefig("data/output/intron_size.png", dpi=500)

# choose only query transcripts
query = query[query[3].str.contains("path1")]
query['length'] = abs(query[1] - query[2])

# chart lcnRNA distribution across chromosome
f, ax = plt.subplots(figsize=(10, 5))
labels, counts = np.unique(query[0].astype(int), return_counts=True)
ax.bar(labels, counts, align='center', color='b')
plt.gca().set_xticks(labels)
plt.ylabel('Number of lncRNAs', size=15)
plt.xlabel('Chromosome of organism', size=15)
plt.savefig("data/output/number_of_lncRNA.png", dpi=500)

chromasom_all = pd.unique(target[0].astype(int))

# choose lncRNA type
name_target = target[9].str.split(';',4, expand=True)
target['name'] = name_target[0]
new_tmap = tmap[tmap[2].isin(['x'])]
# filter seq by lncRNA type
tmap_name = new_tmap[3].str.split('.',4, expand=True)
query_anti = query[query['name'].isin(tmap_name[0].to_list())]
# choose only target exons
exon = target[target[7].str.contains("exon")]
exon = exon[exon['name'].groupby(exon['name']).transform('size')>1]
uni = query_anti.drop_duplicates(subset=[1])

chrom = []
for w in chromasom_all:
    print(w)
    query_new = uni[uni[0].astype(int).isin([w])]
    target_new = exon[exon[0].astype(int).isin([w])]
    query_coord = coord(query_new)

    strand = target_new.groupby([5])
    for key, item in strand:
        if key == "-":
            item['number'] = item.groupby(['name']).cumcount(ascending=False) + 1
            find_same_coords(item)
        elif key == "+":
            item['number'] = item.groupby(['name']).cumcount() + 1
            find_same_coords(item)

df = pd.DataFrame(chrom)
df['number'] = df.groupby([6]).cumcount() + 1
df[6] = df[6].astype(int)
df = df.sort_values(by=[6])
df.rename(columns={0: 'transcript_id', 1: 'protein_coding_id', 2:'protein_coding_id_start', 3:'protein_coding_id_end', 4:'transcript_id_start',
                   5:'transcript_id_end', 6:'protein_coding_exon', 7:'transcript_length', 8:'protein_coding_strand', 9:'transcript_strand', "number":'number_of_transcripts'}, inplace=True)
df.to_csv('data/output/statistic_bed.tsv', sep='\t', index=False)

f, ax = plt.subplots(figsize=(15, 5))
ax.bar(df['protein_coding_exon'], (df['number_of_transcripts']/(len(df['protein_coding_exon']))*100), color='b')
plt.ylabel('Number of aligned lnRNAs', size=15)
plt.xlabel('Exon of the gene', size=15)
plt.xticks(range(max(df['protein_coding_exon'])+1))
plt.savefig('data/output/anti.png')

