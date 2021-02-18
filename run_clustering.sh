python ../korarchaeota-scripts/make_EggNOG_clusters.py --f "../01_proteomes/proteomes/*" -a "../02_eggnog5_annotations/*.annotations" -e clusters_arCOG -o faa_orthofinder

./orthofinder.sif orthofinder -f faa_orthofinder -X -os -M msa -t 20 -a 10 -o OrthoFinder

mv OrthoFinder/Results_Feb04/Orthogroup_Sequences clusters_orthofinder
