from Bio import SeqIO
import glob
import os
import pandas as pd

sp = set()
og = set()
for f in glob.glob('clusters_*/*'):
    og.add(os.path.basename(f.split('.')[0]))
    for rec in SeqIO.parse(f, 'fasta'):
        sp.add(rec.id.split('..')[0])
df = pd.DataFrame(columns=sp, index=og).fillna(0)
for f in glob.glob('clusters_*/*'):
    og = os.path.basename(f.split('.')[0])
    for rec in SeqIO.parse(f, 'fasta'):
        sp = rec.id.split('..')[0]
        df.loc[og, sp] += 1

df = df[df.sum(axis=1) > 3]
df['sum'] = df.sum(axis=1)
df = df.sort_values(by='sum', ascending=False)
df = df.drop(columns='sum')
df.to_csv('larger_clusters_korarchaeota.csv', index=True, header=True, sep='\t')
