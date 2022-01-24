import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

blast = pd.read_csv('data/output/tissue/cons.txt', header=None)
blast_non = pd.read_csv('data/output/tissue/non.txt', header=None)
blast_cod = pd.read_csv('data/output/tissue/cod.txt', header=None)

name = blast[0].str.rsplit('_', 1,  expand=True)
blast['srx'] = name[1]

name_non = blast_non[0].str.rsplit('_', 1,  expand=True)
blast_non['srx'] = name_non[1]

name_cod = blast_cod[0].str.rsplit('_', 1,  expand=True)
blast_cod['srx'] = name_cod[1]

df = pd.read_csv('tissue/SRX_all_org.tsv', sep='\t', header=None)
df = df[df[1].isin([str(sys.argv[1])])]

same = pd.merge(df, blast, how='inner', left_on=[2], right_on=['srx'])
same_non = pd.merge(df, blast_non, how='inner', left_on=[2], right_on=['srx'])
same_cod = pd.merge(df, blast_cod, how='inner', left_on=[2], right_on=['srx'])

transc = same[2].value_counts().reset_index()
transc_non = same_non[2].value_counts().reset_index()
transc_cod = same_cod[2].value_counts().reset_index()

transc = pd.merge(df, transc, how='inner', left_on=[2], right_on=['index'])
transc_non = pd.merge(df, transc_non, how='inner', left_on=[2], right_on=['index'])
transc_cod = pd.merge(df, transc_cod, how='inner', left_on=[2], right_on=['index'])

transc.to_csv('data/output/tissue/transc.csv')
transc_non.to_csv('data/output/tissue/transc_non.csv')
transc_cod.to_csv('data/output/tissue/transc_cod.csv')


pr = same[3].value_counts().reset_index()
pr_non = same_non[3].value_counts().reset_index()
pr_cod = same_cod[3].value_counts().reset_index()
df_tis = df[3].value_counts().reset_index()

pr = pr.rename(columns = {3: 'Conserved'}, inplace = False)
pr_non = pr_non.rename(columns = {3: 'Nonconserved'}, inplace = False)
pr_cod = pr_cod.rename(columns = {3: 'Coding'}, inplace = False)

full_table = pd.merge(pr, pr_cod, how='outer', left_on=['index'], right_on=['index'])
full_table = pd.merge(full_table, pr_non, how='outer', left_on=['index'], right_on=['index'])
full_table = pd.merge(full_table, df_tis, how='inner', left_on=['index'], right_on=['index'])

full_table = full_table.fillna(0)
full_table['con'] = full_table['Conserved'] / full_table[3]
full_table['noncon'] = full_table['Nonconserved'] / full_table[3]
full_table['cod'] = full_table['Coding'] / full_table[3]

full_table['lncRNA conserved'] = full_table['con'] / len(same[3])
full_table['lncRNA nonconseved'] = full_table['noncon'] / len(same_non[3])
full_table['mRNA'] = full_table['cod'] / len(same_cod[3])

full_table = full_table[['index', 'lncRNA conserved', 'lncRNA nonconseved', 'mRNA']]

full_table.to_csv('data/output/tissue/tissue_org.csv')
full_table = full_table.set_index('index')
#print(full_table1)

fig,ax = plt.subplots(1,1,figsize=(8,10))
sns.heatmap(full_table, annot=True, annot_kws={'size': 10}, xticklabels=True, cmap='Blues', ax=ax, cbar_kws={"shrink": .82}, yticklabels=True)
plt.savefig('data/output/tissue/tissue_org.png', bbox_inches='tight')
