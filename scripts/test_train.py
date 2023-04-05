import numpy as np
from sklearn.model_selection import train_test_split
from Bio import SeqIO
import sys
import os

lnc = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))
cds = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))
train_lnc = open('./data/input/test_train/train_lnc.fasta', 'w', encoding='utf-8')
print(len(lnc.keys()))
print(len(cds.keys()))
test_lnc = len(lnc.keys())*0.25
test_mrna = test_lnc*2
n=0
keys = list(cds.keys())
pers_lnc=test_lnc/len((list(lnc.keys())))
Y_train, Y_test = train_test_split(list(lnc.keys()), test_size=pers_lnc, shuffle=True)
print(len(Y_test))
if len(Y_train)*2 <= len(cds.keys())/5:
    train_mrna = len(Y_train)*2
else:
    train_mrna = len(cds.keys())/5.2
print(train_mrna)
for w in Y_train:
    print(lnc[w].format('fasta'), end='', file=train_lnc)

while n <= 4:
    os.mkdir('./data/input/test_train/' + str(n))
    train_cds = open('./data/input/test_train/' + str(n) + '/train_mrna.fasta', 'w', encoding='utf-8')
    test_file = open('./data/input/test_train/' + str(n) + '/test.fasta', 'w', encoding='utf-8')
    compare = open('./data/input/test_train/' + str(n) + '/compare.csv', 'w', encoding='utf-8')
    pers_mrna=train_mrna/(len(keys))
    print(pers_mrna)
    X_train, X_test = train_test_split(keys, test_size=pers_mrna, shuffle=True)
    print('train_set',len(X_train))
    print('test_set',len(X_test))
    if len(X_test) > test_mrna:
    	pers_mrna_test = test_mrna/(len(X_test))
    else:
        true_per = len(X_test)*0.25
        pers_mrna_test = true_per/len(X_test)
    pers_mrna_test = test_mrna/(len(X_test))
    train, test = train_test_split(X_test, test_size=pers_mrna_test, shuffle=True)
    for w in test:
        print(cds[w].format('fasta'), end='', file=test_file)
        print(cds[w].id , 'mrna', sep='\t', file=compare)
    for w in train:
        print(cds[w].format('fasta'), end='', file=train_cds)
    for w in Y_test:
        print(lnc[w].format('fasta'), end='', file=test_file)
        print(lnc[w].id, 'lnc', sep='\t', file=compare)
    keys = [w for w in keys if w not in X_test]
    print('new_set',len(keys))
    n=n+1
