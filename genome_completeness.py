import pandas as pd

names = {line.split()[0]:line.split()[1] for line in open('../05_ALE_reconstruction/species_map.txt')}

df = pd.read_csv('Completeness_KorProject.csv', sep=',')
df = df[df['Phyla'] != 'Idunnarchaeota']


kor = pd.read_csv('Korarchaeota_checkM.csv', sep='\t')
kor = kor[kor['Bin Id (WGS)'].isin(names)]
kor = kor.rename(columns={"Bin Id (WGS)": "Bin Id", "GC%": "GC(%)", "No. of contigs":"Contigs", "N50 contig length":"N50", "Longest contig":"Largest"})
kor['Phyla'] = 'Korarchaeota'
kor['Superphyla'] = 'TACK'

all = pd.concat([df, kor], join="inner")

all.to_csv('genome_stats_all.csv', sep='\t', header=True, index=False)
all['name codes'] = all['Bin Id'].apply(lambda x: names[x])
all['fraction missing'] = all['Completeness'].apply(lambda x: (100-x)/100)

all[['name codes', 'fraction missing']].to_csv('fraction_missing.txt', 
                        sep=':', header=False, index=False, float_format='%.2f')
