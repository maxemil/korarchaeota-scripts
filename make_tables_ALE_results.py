import pandas as pd
from math import modf
import ete3

def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor

# tree_suffix = 'TACK_base'
# tree_suffix = 'TACK_within_SR4'
tree_suffix = 'NM'
# tree_suffix = 'RP56'

# names = {line.split()[1]:line.split()[0] for line in open('species_map.txt')}
names = {line.split()[1]:line.split()[0] for line in open('../05_ALE_reconstruction/species_map.txt')}
groupings = pd.read_csv('../06_genome_completeness/genome_stats_all.csv', 
                        sep='\t', header=0)

# tree = ete3.PhyloTree('species_trees/Kor_noI_TACK_base_clean.tree', format=1)
# tree = ete3.PhyloTree('species_trees/Kor_noI_{}_clean.tree'.format(tree_suffix), format=1)
tree = ete3.PhyloTree('species_trees/{}-ChiSquare-PhyloBayes-EuryRooted_forALE_clean.tree'.format(tree_suffix), format=1)
for l in tree.get_leaves():
    l.name = names[l.name]
for p in set(groupings['Phyla']):
    pspecs = []
    for bin in groupings.loc[groupings['Phyla'] == p, 'Bin Id']:
        pspecs.append(tree.get_leaves_by_name(bin)[0])
    ancestor = get_ancestor(pspecs)
    if not ancestor.is_leaf():
        names[ancestor.name] = "L{}CA".format(p)
    
for p in set(groupings['Superphyla']):
    pspecs = []
    for bin in groupings.loc[groupings['Superphyla'] == p, 'Bin Id']:
        pspecs.append(tree.get_leaves_by_name(bin)[0])
    ancestor = get_ancestor(pspecs)
    if not ancestor.is_leaf():
        names[ancestor.name] = "L{}CA".format(p)

# make mapping tables for all names
# for n in tree.traverse():
#     if not n.name in names and not n.is_leaf():
#         names[n.name] = n.name
# pd.DataFrame(names.values(),index=names.keys()).to_csv(f'{tree_suffix}_mapping_all_names_ancestors.tsv', sep='\t', index=True, header=False)


for n in tree.traverse():
    if n.name in names:
        n.name = names[n.name]
tree.ladderize(direction=1)
def my_layout(n):
    n.img_style['size'] = 0
    name_face = ete3.TextFace(n.name)
    n.add_face(name_face, column=0, position='branch-right')
ts = ete3.TreeStyle()
ts.show_leaf_name = False
ts.show_branch_support = False
ts.layout_fn = my_layout
tree.render('Kor_noI_{}_longnames.pdf'.format(tree_suffix),tree_style=ts)
tree.write(outfile='Kor_noI_{}_longnames.tree'.format(tree_suffix), format=8)

header = ['species_tree', 'cluster', 'node', 'duplications',
          'transfers', 'losses', 'originations', 'copies']
events = pd.read_csv('korarchaeota-arCOG-{}/events.txt'.format(tree_suffix), 
                            sep='\t', names=header)
del events['species_tree']
events['cluster'] = events['cluster'].str.replace('.clean.ufboot.ale', '')

events['node'] = events['node'].apply(lambda x: names[x] if x in names else x)

for e in ['duplications', 'transfers', 'losses', 'originations', 'copies']:
    event = pd.pivot_table(events, values=e, index=['cluster'], columns=['node'])
    event = event.loc[event.sum(axis=1).sort_values(ascending=False).index]
    event.fillna(0, inplace=True)
    event.to_excel('{}_{}_raw.xlsx'.format(e, tree_suffix), engine='openpyxl', 
                header=True, index=True)
    event = event.applymap(lambda x: int(modf(x)[1] + int(modf(x)[0] >= 0.3)))
    event.to_excel('{}_{}.xlsx'.format(e, tree_suffix), engine='openpyxl', 
                header=True, index=True)
