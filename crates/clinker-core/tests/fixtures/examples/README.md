# Parser regression examples

Frozen copies of the GenBank and GFF3/FASTA examples used by parser and workflow
regressions. They live inside the core crate so published source packages can
compile and run their integration tests without a checkout of the repository.
Changes to the public demonstration inputs do not automatically change these
regression expectations. Update a snapshot deliberately when changing a test.
