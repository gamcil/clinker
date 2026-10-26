use clinker_core::{ParseOptions, parse_genbank, parse_gff3};

macro_rules! assert_example_parity {
    ($stem:literal) => {{
        let genbank = parse_genbank(
            concat!($stem, ".gbk"),
            include_bytes!(concat!("../../../examples/", $stem, ".gbk")),
        )
        .expect("bundled GenBank example parses");
        let gff = parse_gff3(
            concat!($stem, ".gff3"),
            include_bytes!(concat!("../../../examples/gff3/", $stem, ".gff3")),
            Some(include_bytes!(concat!(
                "../../../examples/gff3/",
                $stem,
                ".fna"
            ))),
            &ParseOptions::default(),
        )
        .expect("converted GFF3 example parses");

        assert_eq!(gff.name, genbank.name);
        assert_eq!(gff.loci.len(), genbank.loci.len());
        let genbank_genes = genbank
            .loci
            .iter()
            .flat_map(|locus| &locus.genes)
            .collect::<Vec<_>>();
        let gff_genes = gff
            .loci
            .iter()
            .flat_map(|locus| &locus.genes)
            .collect::<Vec<_>>();
        assert_eq!(gff_genes.len(), genbank_genes.len());
        assert_eq!(
            gff_genes
                .iter()
                .map(|gene| &gene.translation)
                .collect::<Vec<_>>(),
            genbank_genes
                .iter()
                .map(|gene| &gene.translation)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            gff_genes
                .iter()
                .map(|gene| (&gene.label, gene.start, gene.end, gene.strand))
                .collect::<Vec<_>>(),
            genbank_genes
                .iter()
                .map(|gene| (&gene.label, gene.start, gene.end, gene.strand))
                .collect::<Vec<_>>()
        );
    }};
}

#[test]
fn converted_gff3_examples_preserve_genes_proteins_and_cds_coordinates() {
    assert_example_parity!("A. alliaceus CBS 536.65");
    assert_example_parity!("A. burnettii MST-FP2249");
    assert_example_parity!("A. mulundensis DSM 5745");
    assert_example_parity!("A. versicolor CBS 583.65");
    assert_example_parity!("P. vexata CBS 129021");
}
