# Example inputs

The five GenBank files are burnettramic-acid (*bua*) cluster homologues used in
the manuscript figure and main README. The *A. burnettii* cluster was reported by
Li et al. (2019), [doi:10.1021/acs.orglett.8b04042](https://doi.org/10.1021/acs.orglett.8b04042).

```sh
clinker examples/*.gbk --plot
```

`gff3/` contains equivalent annotations and same-stem FASTA files. Select matching
`.gff3` / `.fna` pairs in the web app, or run:

```sh
clinker examples/gff3 --no-align --plot gff3-examples.html
```

See [annotations.md](annotations.md) for the illustrative annotation table.
