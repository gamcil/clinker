# GCF 2 annotation table

`annotations.tsv` covers all 169 GenBank fixtures in this directory, keyed by
their local filenames in `source_file`. BGC identifiers and radius values come
from `source.tsv`; collection names and sequence counts come from
`downloads.json`.

After plotting the fixtures in clinker, choose **Add data → Import annotations**
and select this table. Keep **Annotation tracks** for independent scales: use
heatmaps or bar charts for the numeric columns and a categorical heatmap for
`collection`. Identifier and URL columns start skipped. Review the panel and
scale preview before import. These metadata columns have different meanings
and units, so they should not share a heatmap-matrix scale.
