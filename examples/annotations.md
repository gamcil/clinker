# Example annotation table

`annotations.tsv` matches the five GenBank example clusters by their input
filename. `organism`, `biosynthetic_class`, and `reported_bua_homologue` are
labels for demonstrating categorical heatmaps. The cluster-length and gene-count
values are illustrative only; they are not measurements from the example files.

Open the five example GenBank files in clinker, select **Add data → Import annotations**,
and choose `annotations.tsv`. The importer should identify `source_file` as the
matching column automatically. Choose a heatmap or bar chart for the numeric
columns and heatmaps for the categorical columns.
