# clinker native command line

Extract the archive for your operating system, then run `clinker --help`
(`.\clinker.exe --help` in Windows PowerShell). Rust and Python are not
required to run the executable.

Available archives target Linux x86-64, Windows x86-64, and Apple Silicon macOS.
The Linux executable dynamically links glibc and is built on the GitHub
`ubuntu-latest` runner; older Linux distributions and Alpine are not guaranteed
to support it. Intel macOS and Linux ARM binaries are not currently provided.

```sh
./clinker input-a.gbk input-b.gbk --plot
./clinker input-a.gbk input-b.gbk --plot figure.html
./clinker cluster.gff3 --plot figure.html
./clinker input-a.gbk input-b.gbk --session analysis.json
./clinker --session analysis.json --plot figure.html
./clinker input-a.gbk input-b.gbk --project comparison.json
./clinker --open-project comparison.json
./clinker --open-project comparison.json --no-browser
```

GFF3 files need embedded FASTA or a same-stem companion FASTA in the same
directory. If several companion FASTA files match, supply one explicitly alongside
the GFF file. Each input keeps its own directory for pairing, even when different
directories contain identical filenames. Directory discovery is sorted and skips
symlinked directories; repeated paths to the same input are analysed once.

Place input paths before `--plot`: without a filename it starts the
shared local viewer and opens your default browser. The CLI prints a loopback
URL and remains running until Ctrl+C. `--no-browser` prints the URL without
launching a browser, which is useful on remote/headless machines. The server
binds to `127.0.0.1` and serves only embedded assets and this run's JSON data.

The local viewer includes figure editing, tree and annotation import, theme
settings, project saving, SVG and similarity CSV exports. It renders native
results directly; it does not run analysis in WASM. File → Open project can
reopen another project while the CLI is running. Similarity values are derived
natively from the retained analysis; protein alignments are not rerun.

`--open-project` restores edited figure data, trees, annotations, appearance and
layout from a saved project. It does not reconstruct the figure from analysis.
Saving in the viewer downloads a project; it does not overwrite a file on disk.

With a filename, `--plot figure.html` writes the separate self-contained,
renderer-only HTML export, which opens directly in a browser. It does not include
the full editor. `--plot data.json` continues to write raw plot data.

`--project comparison.json` writes a complete version 2 project that opens in the
browser app through **File → Open project**. It retains the analysis data needed
for future exports. `--session` remains the native analysis continuation format;
it does not store interactive figure edits. See [project formats](project-format.md).

When adding files with `--session`, omitted alignment flags inherit the saved
identity cutoff, prefilter parameters and alignment mode. A session created with
`--no-align` continues without alignment. `--jobs` can change freely.

Changing analysis settings requires `--reanalyse`. This rebuilds all comparisons
from saved protein translations, including any newly supplied inputs:

```sh
clinker --session analysis.json --reanalyse --identity 0.5
clinker --session analysis.json --reanalyse --no-prefilter
```

Use `--prefilter` to enable filtering; omitted k-mer parameters retain saved
values, or default to 3 for a new filter. Reanalysis also enables alignment for a
session created with `--no-align`. Conversely, `--reanalyse --no-align` clears
saved links and records a session with alignment disabled. Older sessions without
alignment metadata can still be opened and exported, but require `--reanalyse` before adding files or
changing settings. Their previous filtering history cannot be recovered;
reanalysis defaults to no prefilter unless one is requested.

Output filenames and directories are checked before analysis. Destinations must
be distinct from each other and from input files. Existing exports require
`--force`; output symlinks and read-only files are rejected. Updating an existing
`--session` is intentional and does not require `--force`.

Each output is staged in its destination directory before being published.
Existing files are replaced atomically; new files never overwrite a file that
appeared during the run. The session is saved after the exports. Writes are
atomic per file, rather than a transaction across every output in a run.

The archive includes the clinker, clustermap.js and D3 licence notices.
