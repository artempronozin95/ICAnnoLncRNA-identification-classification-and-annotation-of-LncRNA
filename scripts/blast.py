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
    query = x['qry_gene_id'].str.rsplit('.', 1)
    for w in query:
        try:
            print(record_dict[w[0]].format('fasta'), end='', file=lncrna)
        except KeyError:
            continue

def alingment(x,y):
    path = x.rsplit('/', 1)[0]
    outfmt = os.path.join(path, 'blast' + '.outfmt6')
    aling = 'blastn -query {q} -db {dbw} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"  -evalue 1e-5  -max_target_seqs 1000 -perc_identity 30 -num_threads 1 -out {outfmt}'. format(q=x, dbw=y, outfmt=outfmt)
    #| sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > {outfmt}'. format(q=x, dbw=y, outfmt=outfmt)
    exit_code = subprocess.call(aling, shell=True)

gff_tmap = sys.argv[1]
data = sys.argv[2]
record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[3], "fasta"))
lncrna = open(sys.argv[4], 'w', encoding='utf-8')
tmap = pd.read_csv(gff_tmap, sep='\t')
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
plt.savefig("data/output/classes.png", dpi=500)
fasta(tmap)
indx = data_base(data)
alingment(sys.argv[4], indx)
