import glob
from Bio import SeqIO

with open('species_gene_list.txt', 'w') as out:
    for f in glob.glob('../01_proteomes/proteomes/*.faa'):
        rec = list(SeqIO.parse(f, 'fasta'))[0]
        print(rec.id.split('..')[0], file=out)


import ete3

names = {line.split()[1]:line.split()[0] for line in open('../gene_map.txt')}

species_names = []
tree = ete3.PhyloTree('Kor-noI.fasta.treefile',format=1)
for n in tree.get_leaves():
    species_names.append(n.name)


with open("../species_map.txt", 'w') as out:
    for c, sp in names.items():
        for n in species_names:
            if n == sp:
                print(sp, c, sep='\t', file=out)
            elif n.replace('.', '') == sp:
                print(n, c, sep='\t', file=out)
