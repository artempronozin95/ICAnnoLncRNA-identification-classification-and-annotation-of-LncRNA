import os
import subprocess
import sys
from Bio import SeqIO
import pandas as pd



def data_base(x):
    x_n = x.rsplit('/', 1)[1]
    x_n = x_n.rsplit('.', 1)[0]
    dir_check = os.path.isfile('data/reference/data_base/' + x_n + '.dmnd')
    if dir_check is True:
        print('index exist')
        index_path = os.path.join('data/reference/data_base/', x_n + '.dmnd')
        return index_path
    else:
        print('build index')
        index_path = os.path.join('data/reference/data_base/', x_n + '.dmnd')
        command = 'diamond makedb --in {db} -d {index}'. format(db=x, index=index_path)
        exit_code = subprocess.call(command, shell=True)
        return index_path

def alingment(x,y):
    path = x.rsplit('/', 1)[0]
    out_path = os.path.join(path, 'diamond' + '.daa')
    outfmt = os.path.join(path, 'diamond' + '.outfmt6')
    aling = 'diamond blastx -q {q} -d {dbw} -e 1e-3 --threads 10 -a {daa}'. format(q=x, dbw=y, daa= out_path)
    exit_code = subprocess.call(aling, shell=True)
    transform = "diamond view -a {daa} -o {outfmt6}". format(daa= out_path, outfmt6=outfmt)
    exit_code = subprocess.call(transform, shell=True)

data = sys.argv[1]
query = sys.argv[2]
indx = data_base(data)
alingment(query, indx)
