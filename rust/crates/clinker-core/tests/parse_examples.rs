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
