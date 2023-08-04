import pandas as pd
import os
import sys


strc = str(sys.argv[3])[1:-1]

if strc == 'own':
   print('Own table')
   with open(sys.argv[9], 'w') as f:
        f.write('Own table')
else:
  print('Build table')
  tissue = pd.read_csv(sys.argv[1], sep='\t', header=None)
  path = str(sys.argv[2])[1:-1]
  new_old =  pd.read_csv(sys.argv[4], sep='\t')
  cons =  pd.read_csv(sys.argv[5], sep=' ', header=None)
  noncons =  pd.read_csv(sys.argv[6], sep=' ', header=None)
  cod = pd.read_csv(sys.argv[7], sep=' ', header=None)
  cod = pd.merge(cod, new_old, how='inner', left_on=[0], right_on=['new_id']).drop(columns = [0])
  cod['type'] = 'mRNA'
  noncons['type'] = 'noncons'
  cons['type'] = 'consv'

  if strc == 'Kallisto':
     list_dir = os.listdir(path)
     lib = {}
     for w in list_dir:
         abudance = pd.read_csv(path + '/' + w +'/' + "abundance.tsv", sep='\t')
         lib[w] = abudance['tpm']
  if strc == 'Htseq':
     list_dir = os.listdir(path)
     lib = {}
     for w in list_dir:
         id = w.split('.')
         abudance = pd.read_csv(path + '/' + w, sep='\t', header=None)
         abudance = abudance.rename(columns={0: 'target_id'})
         lib[id[0]] = abudance[1]

  lib_data = pd.DataFrame.from_dict(lib)
  colnames = lib.keys()
  tissue = tissue[tissue[1].isin(colnames)]
  tissue_gr = tissue.groupby([2])
  median_tis = {}
  for k,v in tissue_gr:
     lib_data_gr = lib_data[v[1]]
     median_tis[k] = lib_data_gr.median(axis=1)

  median_tis_df = pd.DataFrame.from_dict(median_tis)
  median_tis_df['target_id'] = abudance['target_id']
  cod = pd.merge(cod, median_tis_df, how='inner', left_on=['new_id'], right_on=['target_id']).drop(columns = ['old_id', 'new_id'])
  noncons = pd.merge(noncons, median_tis_df, how='inner', left_on=[0], right_on=['target_id']).drop(columns = [0])
  cons = pd.merge(cons, median_tis_df, how='inner', left_on=[0], right_on=['target_id']).drop(columns = [0])
  full_table = pd.concat([cod, cons, noncons])
  full_table.to_csv(sys.argv[8], sep='\t', index=False)
  with open(sys.argv[9], 'w') as f:
        f.write(strc)
