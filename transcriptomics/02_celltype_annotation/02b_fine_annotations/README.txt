
The coarse cell types were annotated in 02a_cluster_annotation.ipynb (additionally in final_annotations.ipynb there are
small adjustments to the coarse annotations).

For each coarse cell type (level1) we now annotate finer cell types (level2) --> 1 ipynb for each coarse cell type to 
cluster the sub states and investigate marker genes + trajectories.

Finally in the notebook final_annotations.ipynb all level2 annotations are summarised and saved.

In 02c_final_annotations_umaps_colors.ipynb an adata object is generated with the fine annotations and secondary umaps
for (groups) of coarse annotations.