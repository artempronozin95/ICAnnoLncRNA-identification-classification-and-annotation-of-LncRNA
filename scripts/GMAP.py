import os
import sys
import subprocess
import pandas as pd

data_base = sys.argv[1]
target = sys.argv[2]
parameters = sys.argv[3]
out_file = sys.argv[5]
opt = sys.argv[7:]


def aling(reference, transcripts, min, max, thr):
    gmap_run = 'gmap'

    path, base = os.path.split(reference)
    ref_label = os.path.splitext(base)[0]
    out = os.path.join(out_file)
    gmap_run_logger_err_path = os.path.join('./data/output/', gmap_run + '.out.log')
    command = '{gmap} -D {tmp_dir} -d {ref_index_name} {transcripts} --min-intronlength={min_intron_length} --intronlength={intron_length} --cross-species ' \
                  '--format=gff3_gene --split-large-introns --npaths=1 -t {threads} {option} > {alignment_out} 2>> {log_out_2}'.\
            format(gmap=gmap_run, tmp_dir=path, ref_index_name=ref_label, transcripts=transcripts,
                   threads=thr, alignment_out=out, log_out_2=gmap_run_logger_err_path, min_intron_length=min[0], intron_length=max[0], option=' '.join('{0}'.format(w) for w in opt))
    exit_code = subprocess.call(command, shell=True)

dir_check = os.path.isdir('data/reference/intron_data/' + sys.argv[4])
if dir_check is True:
   print('use of existing statistics')
   gff_df = pd.read_csv('data/reference/intron_data/' + sys.argv[4] + '/' + sys.argv[4] + '.tsv', sep='\t', header=None)
   minn = gff_df[gff_df[0].str.contains("Min_intron_length")][1].to_list()
   maxx = gff_df[gff_df[0].str.contains("Max_intron_length")][1].to_list()
   thr = 40
   aling(data_base, target, minn, maxx, thr)
   print(target)
else:
   print('use of reference statistics')
   gff_df = pd.read_csv(parameters, sep='\t', header=None)
   minn = gff_df[gff_df[0].str.contains("Min_intron_length")][1].to_list()
   maxx = gff_df[gff_df[0].str.contains("Max_intron_length")][1].to_list()
   thr = 40
   aling(data_base, target, minn, maxx, thr)
