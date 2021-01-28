from Bio import SeqIO
import glob
import os

for f in glob.glob("Subset/*.faa"):
    id = os.path.basename(f).replace('.faa', '').replace('.', _)
    with open('{}.faa'.format(id), 'w') as out:
        for rec in SeqIO.parse(f, 'fasta'):
            rec.id = "{}..{}".format(id, rec.id)
            rec.description = ""
            SeqIO.write(rec, out, 'fasta')
