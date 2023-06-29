import os
import sys
import subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np
from gtfparse import read_gtf

data_base = sys.argv[1]
gff = sys.argv[2]


def build(reference, kmer):
    args_alignment = []
    gmap_build = 'gmap_build'
    path, base = os.path.split(reference)
    ref_label = os.path.splitext(base)[0]
    gmap_build_logger_out_path = os.path.join('./data/output/GMAP', gmap_build + '.out.log')
    gmap_build_logger_err_path = os.path.join('./data/output/GMAP', gmap_build + '.err.log')


    command = '{gmap_build} -D {tmp_dir} -d {ref_index_name} -k {kmer_value} {reference} >> {log_out_1} 2>> {log_out_2}'.\
            format(gmap_build=gmap_build, tmp_dir=path, ref_index_name=ref_label, reference=reference, log_out_1=gmap_build_logger_out_path,
                   log_out_2=gmap_build_logger_err_path, kmer_value=kmer)
    exit_code = subprocess.call(command, shell=True)


gff_df = read_gtf(gff)
df = SeqIO.to_dict(SeqIO.parse(data_base, "fasta"))
size = round(os.path.getsize(data_base) / (1024*1024))
statistic = open(sys.argv[3], 'w', encoding='utf-8')

if size <= 64:
    k = 12
elif 64 < size <= 256:
    k = 13
elif 256 < size <= 1000:
    k = 14
elif 1000 < size <= 10000:
    k = 15
print("k" ,'\t', str(k), file=statistic)

total_len = []
for name,seq in df.items():
    total_len.append(len(seq))
print('Golden_Path_Length' ,'\t', str(np.sum(total_len)), file=statistic)

exon_only = gff_df[gff_df['feature'].str.contains('exon')]
exon_length = (exon_only['end'] - exon_only['start']).mean()
print('Mean_exon_length','\t', str(round(exon_length)), file=statistic)
transcript_only = gff_df[gff_df['feature'].str.contains('transcript')]
transcript_length = (transcript_only['end'] - transcript_only['start']).mean()
print('Mean_transcript_length','\t', str(round(transcript_length)), file=statistic)

gr = exon_only[['start','end','gene_id']].groupby(by='gene_id')
intron_length = []
for key, item in gr:
    n=0
    while n <= len(item['gene_id']) - 2:
        intron = abs(item['end'].iloc[n] - item['start'].iloc[n + 1]) - 1
        n = n + 1
        intron_length.append(intron)

mean_intron_length = np.mean(intron_length)
max_intron_length = np.max(intron_length)
min_intron_length = np.min(intron_length)
print('Mean_intron_length' ,'\t', str(round(mean_intron_length)), file=statistic)
print('Max_intron_length' ,'\t', str(round(max_intron_length)), file=statistic)
print('Min_intron_length' ,'\t', str(round(min_intron_length)), file=statistic)

build(data_base, k)
