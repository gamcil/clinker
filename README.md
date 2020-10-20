# clinker
Gene cluster comparison figure generator

## What is it?
clinker is a pipeline for easily generating publication-quality gene cluster
comparison figures.

<img src="images/figure.png" alt="bua cluster and homologues" width=700>

Given a set of GenBank files, clinker will automatically extract protein translations,
perform global alignments between sequences in each cluster, determine the
optimal display order based on cluster similarity, and generate an interactive
visualisation (using ![clustermap.js](https://github.com/gamcil/clustermap.js))
that can be extensively tweaked before being exported as an SVG file.

## Installation
clinker can be installed directly through pip:

`pip install clinker`

Or by cloning the source code from GitHub:

```
git clone https://github.com/gamcil/clinker.git
cd clinker
pip install .
```

## Usage
Running clinker can be as simple as:

`clinker clusters/*.gbk`

This will run the clinker pipeline on all GenBank files within the
folder and open a new tab in your web browser with the visualisation
application.

See `-h/--help` for more information:

```
usage: clinker [-h] [-i IDENTITY] [-f] [-o OUTPUT] [-dl DELIMITER]
               [-dc DECIMALS] [-hl] [-ha]
               files [files ...]

clinker: Automatic creation of publication-ready gene cluster comparison figures.

clinker generates gene cluster comparison figures from GenBank files.
It performs pairwise local or global alignments between every sequence
in every unique pair of clusters and generates interactive, to-scale
comparison figures using the clustermap.js library.

positional arguments:
  files                 Gene cluster GenBank files

optional arguments:
  -h, --help            show this help message and exit

Alignment options:
  -i IDENTITY, --identity IDENTITY
                        Minimum alignment sequence identity

Output options:
  -f, --force           Overwrite previous output file
  -o OUTPUT, --output OUTPUT
                        Save alignments to file
  -dl DELIMITER, --delimiter DELIMITER
                        Character to delimit output by
  -dc DECIMALS, --decimals DECIMALS
                        Number of decimal places in output
  -hl, --hide_link_headers
                        Hide alignment column headers
  -ha, --hide_aln_headers
                        Hide alignment cluster name headers

Example usage
-------------
Align clusters, plot results and print scores to screen:
  $ clinker files/*.gbk

Cameron Gilchrist, 2020
```
