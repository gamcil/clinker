# Example annotation table

`annotations.tsv` matches the five GenBank example clusters by their input
filename. `organism`, `biosynthetic_class`, and `reported_bua_homologue` are
labels for demonstrating categorical heatmaps. The cluster-length and gene-count
values are illustrative only; they are not measurements from the example files.

Open the five example GenBank files in clinker, select **Add data → Import annotations**,
and choose `annotations.tsv`. The importer should identify `source_file` as the
matching column automatically. Keep **Separate tracks** selected: the numeric
columns start as separate bar charts, and categorical columns as separate
heatmaps. Cluster length and gene count retain independent scales because they
have different units. Review the scale beside each column before importing.

For a table of comparable measurements, such as expression values for several
genes, choose **One heatmap** instead. Selected columns form one panel with a
shared scale. Numeric columns start selected; text columns start skipped.
Automatic scaling recognises numeric 0/1 values as discrete presence/absence;
other numeric values use a continuous scale. Categories can be selected
explicitly to share a category legend. Blank cells stay missing, and values are
never normalised automatically. This mixed metadata example is intended for
tracks, rather than a single shared-scale matrix.
