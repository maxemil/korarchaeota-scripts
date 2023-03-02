from Bio import SeqIO
import glob
import os

with open('protein_to_cluster_id.tsv', 'w') as out:
    for f in glob.glob('clusters_*/*'):
        og = os.path.basename(f.split('.')[0])
        for rec in SeqIO.parse(f, 'fasta'):
            sp = rec.id.split('..')[0]
            print(rec.id.split('..')[1], sp, og, sep='\t', file=out)
