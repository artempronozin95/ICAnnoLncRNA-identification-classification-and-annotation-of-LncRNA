import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t', header=None)

trans = df[4].str.split('=', 2, expand=True)
df[6] = trans[1]
dict =pd.Series(df[6].values,index=df[0])

not_trans = []
trans = []
for k,v in dict.items():
    if int(v) > 0:
        k = k.rsplit('_', 1)[0]
        trans.append(k)
    else:
        k = k.rsplit('_', 1)[0]
        not_trans.append(k)

uni_trans = list(set(trans))
uni_nottrans = list(set(not_trans))
same = set(uni_trans).intersection(uni_nottrans)
nottransmem = pd.DataFrame(uni_nottrans)
transmem = pd.DataFrame(uni_trans)
nottransmem = nottransmem[~nottransmem[0].isin(same)]
nottransmem.to_csv(sys.argv[3], index=None, header=None)
transmem.to_csv(sys.argv[2], index=None, header=None)


