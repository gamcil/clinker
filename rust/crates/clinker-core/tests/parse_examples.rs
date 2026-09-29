use clinker_core::parse_genbank;

#[test]
fn parses_the_existing_versicolor_example_into_one_locus() {
    let bytes = include_bytes!("../../../../examples/A. versicolor CBS 583.65.gbk");

    let cluster = parse_genbank("A. versicolor CBS 583.65.gbk", bytes)
        .expect("the bundled example is valid GenBank");

    assert_eq!(cluster.name, "A. versicolor CBS 583.65");
    assert_eq!(cluster.loci.len(), 1);
    assert_eq!(cluster.loci[0].name, "KV878126");
    assert!(cluster.loci[0].end > 10_000);
    assert_eq!(cluster.loci[0].genes.len(), 8);

    let first = &cluster.loci[0].genes[0];
    // clinker gives protein_id precedence over locus_tag when choosing labels.
    assert_eq!(first.label, "OJI99004.1");
    assert_eq!(first.start, 0);
    assert_eq!(first.end, 1_355);
    assert_eq!(first.strand, -1);
    assert!(first.translation.starts_with("MHKPTILF"));
}

#[test]
fn derives_translations_for_forward_and_reverse_cds_features() {
    let bytes = br#"LOCUS       TEST                      18 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="forward"
     CDS             complement(10..18)
                     /locus_tag="reverse"
ORIGIN
        1 atggcttaattaagccat
//
"#;

    let cluster = parse_genbank("derived.gbk", bytes).expect("valid GenBank");
    let genes = &cluster.loci[0].genes;

    assert_eq!(genes[0].translation, "MA*");
    assert_eq!(genes[0].strand, 1);
    assert_eq!(genes[1].translation, "MA*");
    assert_eq!(genes[1].strand, -1);
}

#[test]
fn handles_reverse_spliced_cds_with_complements_inside_a_join() {
    let bytes = include_bytes!("../../../../examples/A. burnettii MST-FP2249.gbk");
    let cluster = parse_genbank("A. burnettii MST-FP2249.gbk", bytes)
        .expect("the bundled example is valid GenBank");
    let genes = &cluster.loci[0].genes;
    let atg22 = genes
        .iter()
        .find(|gene| gene.label == "ncbi:ETB97_008321-T1")
        .expect("ATG22 CDS is present");

    // The location is join(complement(5441..5600), complement(5095..5373),
    // complement(3634..5024)); its genomic envelope must not be inverted.
    assert_eq!((atg22.start, atg22.end, atg22.strand), (3_633, 5_600, -1));
    assert!(genes.iter().all(|gene| gene.start <= gene.end));
}
