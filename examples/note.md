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
clinker examples/* -p
```
