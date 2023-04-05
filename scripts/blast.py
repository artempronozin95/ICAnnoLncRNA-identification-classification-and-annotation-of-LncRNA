

import os
import subprocess
import sys
from Bio import SeqIO
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def data_base(x):
    x_n = x.rsplit('/', 1)[1]
    x_n = x_n.rsplit('.', 1)[0]
    dir_check = os.path.isfile('data/reference/data_index/' + x_n + '.nin')
    if dir_check is True:
        print('index exist')
        index_path = os.path.join('data/reference/data_index/', x_n)
        return index_path
    else:
        print('build index')
        index_path = os.path.join('data/reference/data_index/', x_n)
        command = 'makeblastdb -in {db} -dbtype nucl -parse_seqids -out {index}'. format(db=x, index=index_path)
        exit_code = subprocess.call(command, shell=True)
        return index_path

def fasta(x):
    query = x['name']
    for w in query:
        try:
            old_id = old_new[old_new['new_id'].isin([w])]
            print(record_dict[str(old_id['old_id'].to_list()[0])].format('fasta'), end='', file=lncrna)
        except KeyError:
            continue

def alingment(x,y):

    path = x.rsplit('/', 1)[0]
    outfmt = os.path.join(path, 'blast' + '.outfmt6')
    aling = 'blastn -query {q} -db {dbw} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"  -evalue {evalue}  -max_target_seqs {max_target} -perc_identity {identity} -num_threads {threads} -out {outfmt}'. format(q=x, dbw=y, evalue=evalue, max_target=max_target, identity=identity, threads =threads, outfmt=outfmt)
    #| sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > {outfmt}'. format(q=x, dbw=y, outfmt=outfmt)
    exit_code = subprocess.call(aling, shell=True)

# input data
query = pd.read_csv(sys.argv[1], sep='\t', header=None)
gff_tmap = sys.argv[2]
data = sys.argv[3]
record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[4], "fasta"))
lncrna = open(sys.argv[5], 'w', encoding='utf-8')
evalue = sys.argv[6]
max_target = sys.argv[7]
identity = sys.argv[8]
threads = sys.argv[9]
old_new = pd.read_csv(sys.argv[10], sep='\t')
tmap = pd.read_csv(gff_tmap, sep='\t')
# lncRNA classification graf
query = query[query[7].str.contains("gene")]
tmap = tmap[tmap['qry_gene_id'].str.contains("path1")]
name_tmap = tmap['qry_gene_id'].str.rsplit('.',1, expand=True)
tmap[['name','seq']] = name_tmap[[0,1]]
#tmap = tmap[tmap['seq'].str.contains("path1")]
tmap = tmap[tmap['name'].isin(query[10])]
tmap= tmap.drop_duplicates(subset=['name'])
tmap['group'] = np.nan
tmap['group'][tmap['class_code'].isin(['x'])] = "exon antisense"
tmap['group'][tmap['class_code'].isin(['i'])] = "intron"
tmap['group'][tmap['class_code'].isin(['u'])] = "intergenic"
tmap.dropna(subset = ["group"], inplace=True)
print(tmap['group'].value_counts())
f, ax = plt.subplots(figsize=(13, 13))
tmap['group'].value_counts().plot(kind='bar')
plt.xticks(size=15, rotation=30)
plt.yticks(size=15)
plt.ylabel('LncRNA transcripts number', size=15)
plt.xlabel('LncRNA classes', size=15)
plt.savefig("data/output/new_lncRNA/classes.png", dpi=500)
# BLAST
fasta(tmap)
indx = data_base(data)
alingment(sys.argv[5], indx)
