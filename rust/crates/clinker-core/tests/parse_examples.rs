use std::collections::BTreeMap;

use clinker_core::{GenBankOptions, SequenceRange, parse_genbank, parse_genbank_with_options};

#[test]
fn parses_the_existing_versicolor_example_into_one_locus() {
    let bytes = include_bytes!("../../../../examples/A. versicolor CBS 583.65.gbk");

    let cluster = parse_genbank("A. versicolor CBS 583.65.gbk", bytes)
        .expect("the bundled example is valid GenBank");

    assert_eq!(cluster.name, "A. versicolor CBS 583.65");
    assert_eq!(
        cluster.source_file.as_deref(),
        Some("A. versicolor CBS 583.65.gbk")
    );
    assert_eq!(
        cluster.organism.as_deref(),
        Some("Aspergillus versicolor CBS 583.65")
    );
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

    let adjacent = &cluster.loci[0].genes[5..=6];
    assert_eq!((adjacent[0].start, adjacent[0].end), (9_734, 11_362));
    assert_eq!((adjacent[1].start, adjacent[1].end), (11_385, 11_979));
    assert!(adjacent[0].end <= adjacent[1].start);
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

    assert_eq!(genes[0].translation, "MA");
    assert_eq!(genes[0].strand, 1);
    assert_eq!(genes[1].translation, "MA");
    assert_eq!(genes[1].strand, -1);
}

#[test]
fn uses_cds_coordinates_instead_of_the_enclosing_gene() {
    let bytes = br#"LOCUS       PARENT                    15 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     gene            2..12
                     /locus_tag="parent"
     CDS             4..9
                     /locus_tag="parent"
                     /translation="MA"
ORIGIN
        1 aaaatggcttaaaaa
//
"#;

    let cluster = parse_genbank("parent.gbk", bytes).expect("valid GenBank");
    let gene = &cluster.loci[0].genes[0];

    assert_eq!((gene.start, gene.end), (3, 9));
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

#[test]
fn can_emit_records_as_separate_clusters() {
    let bytes = br#"LOCUS       FIRST                      9 bp    DNA     linear   UNA 01-JAN-2000
SOURCE      First organism
  ORGANISM  First organism
            Test lineage.
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="first"
                     /translation="MA"
ORIGIN
        1 atggcttaa
//
LOCUS       SECOND                     9 bp    DNA     linear   UNA 01-JAN-2000
SOURCE      Second organism
  ORGANISM  Second organism
            Test lineage.
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="second"
                     /translation="MA"
ORIGIN
        1 atggcttaa
//
"#;
    let options = GenBankOptions {
        as_separate_clusters: true,
        ..GenBankOptions::default()
    };

    let clusters = parse_genbank_with_options("records.gbk", bytes, &options).unwrap();

    assert_eq!(clusters.len(), 2);
    assert_eq!(clusters[0].name, "First organism");
    assert_eq!(clusters[1].name, "Second organism");
    assert_eq!(clusters[0].organism.as_deref(), Some("First organism"));
    assert_eq!(clusters[0].source_file.as_deref(), Some("records.gbk"));
    assert_eq!(clusters[0].loci[0].name, "FIRST");
    assert_eq!(clusters[1].loci[0].name, "SECOND");
    assert_eq!(clusters[0].loci[0].genes[0].label, "first");
    assert_eq!(clusters[1].loci[0].genes[0].label, "second");
}

#[test]
fn range_selection_keeps_only_fully_contained_features() {
    let bytes = br#"LOCUS       RANGED                    18 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..6
                     /locus_tag="outside"
                     /translation="M"
     CDS             7..12
                     /locus_tag="inside"
                     /translation="M"
     CDS             12..18
                     /locus_tag="overlap"
                     /translation="M"
ORIGIN
        1 atggctatggctatggct
//
"#;
    let options = GenBankOptions {
        ranges: BTreeMap::from([("RANGED".to_owned(), SequenceRange { start: 6, end: 12 })]),
        ..GenBankOptions::default()
    };

    let clusters = parse_genbank_with_options("range.gbk", bytes, &options).unwrap();
    let locus = &clusters[0].loci[0];

    assert_eq!((locus.start, locus.end), (6, 12));
    assert_eq!(locus.genes.len(), 1);
    assert_eq!(locus.genes[0].label, "inside");
    assert_eq!((locus.genes[0].start, locus.genes[0].end), (6, 12));
}

#[test]
fn rotates_a_fixture_backed_circular_origin_crossing_feature() {
    let bytes = br#"LOCUS       CIRCULAR                  30 bp    DNA     circular UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             join(25..30,1..9)
                     /locus_tag="wraps"
                     /translation="MMMMM"
     CDS             10..18
                     /locus_tag="middle"
                     /translation="MMM"
ORIGIN
        1 atgatgatgatgatgatgatgatgatgatg
//
"#;

    let shifted =
        parse_genbank_with_options("circular.gbk", bytes, &GenBankOptions::default()).unwrap();
    let unchanged = parse_genbank_with_options(
        "circular.gbk",
        bytes,
        &GenBankOptions {
            set_origin: false,
            ..GenBankOptions::default()
        },
    )
    .unwrap();

    let shifted_wrap = &shifted[0].loci[0].genes[0];
    let unchanged_wrap = &unchanged[0].loci[0].genes[0];
    assert_eq!((shifted_wrap.start, shifted_wrap.end), (0, 15));
    assert_eq!((unchanged_wrap.start, unchanged_wrap.end), (0, 30));
    assert_eq!(shifted_wrap.translation, unchanged_wrap.translation);
    assert!(
        shifted[0].loci[0]
            .genes
            .iter()
            .all(|gene| gene.start <= gene.end)
    );
}
