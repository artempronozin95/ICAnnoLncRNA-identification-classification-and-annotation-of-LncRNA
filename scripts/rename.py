from Bio import SeqIO
import sys
import pandas as pd

fasta = SeqIO.parse(sys.argv[1], "fasta")
corrected_file = sys.argv[2]
new_old_id = []
with open(corrected_file, 'w') as corrected:

 n = 1
 for rec in fasta:
    if len(rec.seq) >= 200:
        new_name =  "Transcript_" + str(n)
        new_old_id.append([rec.id, new_name])
        rec.id = new_name
        rec.description = new_name  # <- Add this line
        n = n + 1
        SeqIO.write(rec, corrected, 'fasta')
    else:
        pass

new_old_id_df = pd.DataFrame(new_old_id, columns=['old_id', 'new_id'])
new_old_id_df.to_csv(sys.argv[3], sep='\t', index=False)