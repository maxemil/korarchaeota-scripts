import pandas as pd
import ete3

def get_ancestor(nodes):
    ancestor = nodes[0]
    for n in nodes:
        ancestor = ancestor.get_common_ancestor(n)
    return ancestor

species_map = {line.split()[0]: line.split()[1] for line in open("species_map_NM.txt")}
treeO = ete3.Tree("species_trees/NM-ChiSquare-PhyloBayes-EuryRooted_forALE_root.tree")
treeC = ete3.Tree("species_trees/NM-ChiSquare-PhyloBayes-EuryRooted_forALE_clean.tree")

for n in treeO.traverse():
        if n.is_leaf():
            pass
        else:
            anc = get_ancestor(
                [treeC.get_leaves_by_name(species_map[l.name])[0] for l in n.get_leaves()])
            #n.name = name
            n.name = str(int(anc.support))

# needs to be added manually for now:
#treeO.name = str(int(treeC.support))

with open("NM-ChiSquare-PhyloBayes-EuryRooted_forALE_bl_names.tree", 'w') as outhandle:
    print(treeO.write(format=3), file=outhandle)
