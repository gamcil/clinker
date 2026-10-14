# Example clusters used in manuscript figure
The GenBank files in this folder are the ones used to generate the figure shown
in the manuscript and in the main repository README.

These files are all homologues of the burnettramic acids gene cluster (*bua*)
previously reported by our group (Li et al. 2019, doi: [10.1021/acs.orglett.8b04042](https://doi.org/10.1021/acs.orglett.8b04042)).

They are:
	- A. burnettii MST-FP2249.gbk (*bua*)
	- A. alliaceus CBS 536.65.gbk
	- A. mulundensis DSM 5745.gbk
	- A. versicolor CBS 583.65.gbk
	- P. vexata CBS 129021.gbk

To generate the visualisation shown in the manuscript, the following command was used:

```
clinker examples/*.gbk -p
```

Equivalent GFF3 annotations and same-stem FASTA files are in `examples/gff3/`.
They preserve the records, gene labels, and protein translations from these
GenBank files. To test the Rust parser against all five converted inputs:

```
cd rust
cargo run -p clinker-cli -- ../examples/gff3 --no-align --plot gff3-examples.html
```

In the web application, select one or more matching `.gff3` and `.fna` pairs
from that directory.
