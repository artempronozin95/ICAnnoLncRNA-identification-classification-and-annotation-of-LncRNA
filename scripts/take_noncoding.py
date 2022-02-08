import pandas as pd
import numpy as np
from Bio import SeqIO
import sys

lncfinder = pd.read_csv(sys.argv[1], sep=',', header=None)
record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))
lncrna = open(sys.argv[3], 'w', encoding='utf-8')
coding = open('data/output/Coding.fasta', 'w', encoding='utf-8')

lncfinder = lncfinder[[0,1]]
finder_cod = lncfinder[lncfinder[1].str.contains("NonCoding")==False]
finder_non = lncfinder[lncfinder[1].str.contains("NonCoding")]

for w in finder_non[0]:
    try:
      print(record_dict[w].format('fasta'), end='', file=lncrna)
    except KeyError:
      continue

for w in finder_cod[0]:
    try:
      print(record_dict[w].format('fasta'), end='', file=coding)
    except KeyError:
      continue
