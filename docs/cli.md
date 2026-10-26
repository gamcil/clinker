# clinker native command line

Extract the archive for your operating system, then run `clinker-rs --help`
(`.\clinker-rs.exe --help` in Windows PowerShell). Rust and Python are not
required to run the executable.

Available archives target Linux x86-64, Windows x86-64, and Apple Silicon macOS.
The Linux executable dynamically links glibc and is built on the GitHub
`ubuntu-latest` runner; older Linux distributions and Alpine are not guaranteed
to support it. Intel macOS and Linux ARM binaries are not currently provided.

```sh
./clinker-rs input-a.gbk input-b.gbk --plot figure.html
./clinker-rs cluster.gff3 --plot figure.html
./clinker-rs input-a.gbk input-b.gbk --session analysis.json
./clinker-rs --session analysis.json --plot figure.html
./clinker-rs input-a.gbk input-b.gbk --project comparison.json
```

GFF3 files need embedded FASTA or a same-stem companion FASTA in the same
directory. `--plot` requires an output filename; HTML figures bundle their
rendering assets and can be opened directly in a browser.

`--project comparison.json` writes a complete version 2 project that opens in the
browser app through **File → Open project**. It retains the analysis data needed
for future exports. `--session` remains the native analysis continuation format;
it does not store interactive figure edits. See [project formats](project-format.md).

The archive includes the clinker, clustermap.js and D3 licence notices.
