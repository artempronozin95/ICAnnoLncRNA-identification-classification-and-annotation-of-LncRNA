import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def value_prepare(x, y, z):
    name = x[0].str.rsplit('_', 1, expand=True)
    x['srx'] = name[1]
    same = pd.merge(df, x, how='inner', left_on=[2], right_on=['srx'])
    transcript = same[2].value_counts().reset_index()
    transcript = pd.merge(df, transcript, how='inner', left_on=[2], right_on=['index'])
    transcript.to_csv('data/output/tissue/' + y + '.csv')
    pr = same[3].value_counts().reset_index()
    pr = pr.rename(columns={3: z}, inplace=False)
    return (pr, same)


#LncAPDB_vs_blast table
blast = pd.read_csv(sys.argv[2], sep='\t', header=None)
LncAPDB = pd.read_csv(sys.argv[3], sep=' ', header=None)
old_new = pd.read_csv(sys.argv[4], sep='\t', header=None, skiprows=1)
LncAPDB = LncAPDB.drop_duplicates(subset=[1])
blast = blast[blast[2].isin(['100.000'])]
LncAPDB_vs_blast = pd.merge(blast, LncAPDB, how='inner', left_on=[1], right_on=[1])
LncAPDB_vs_blast = LncAPDB_vs_blast[['0_x',1,'2_x',10, '0_y','2_y']]
LncAPDB_vs_blast.rename(columns={'0_x': 'transcript_id', 1: 'library_id', '2_x':'percent_identity', 10:'e_value','0_y':'id_of_database', '2_y':'database'}, inplace=True)
LncAPDB_vs_blast.to_csv('data/output/new_lncRNA/LncAPDB_vs_blast.csv', sep='\t', index=False)

#tissue analysis
blast = pd.read_csv('data/output/tissue/cons.txt', header=None)
blast_non = pd.read_csv('data/output/tissue/non.txt', header=None)
blast_cod = pd.read_csv('data/output/tissue/cod.txt', header=None)
df = pd.read_csv('tissue/SRX_all_org.tsv', sep='\t', header=None)
df = df[df[1].isin([str(sys.argv[1])])]
df_tis = df[3].value_counts().reset_index()
blast_cod = old_new[old_new[1].isin(blast_cod[0])]

try:

    transc, same = value_prepare(blast, 'transc' , 'Conserved')
    transc_non, same_non = value_prepare(blast_non, 'transc_non', 'Nonconserved')
    transc_cod, same_cod = value_prepare(blast_cod, 'transc_cod', 'Coding')

    full_table = pd.merge(transc, transc_cod, how='outer', left_on=['index'], right_on=['index'])
    full_table = pd.merge(full_table, transc_non, how='outer', left_on=['index'], right_on=['index'])
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

    fig,ax = plt.subplots(1,1,figsize=(8,10))
    sns.heatmap(full_table, annot=True, annot_kws={'size': 10}, xticklabels=True, cmap='Blues', ax=ax, cbar_kws={"shrink": .82}, yticklabels=True)
    plt.savefig('data/output/tissue/tissue_org.png', bbox_inches='tight')

except KeyError:
    print('Please check format of your sequence ID. It should be look like "MSTRG.1031.1_SRX123456" or "TRINITY_DN195_c0_g1_i1_SRX339783". The defining parameter is "_SRX#######" add it after your main ID.')

