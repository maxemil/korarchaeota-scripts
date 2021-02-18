import glob
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO

counts = []
for f in glob.glob('clusters_*/*'):
    count = len(list(SeqIO.parse(f, 'fasta')))
    counts.append(count)

c4 = [c for c in counts if c > 3]

print("number of clusters:", len(counts))
print("number of clusters over 3 members: ", len(c4))
print("largest cluster: ", max(c4))


sns.histplot(c4)
plt.show()
