use std::collections::BTreeMap;

use clinker_core::{ParseOptions, SequenceRange, parse_genbank, parse_genbank_with_options};

#[test]
fn parses_the_existing_versicolor_example_into_one_locus() {
    let bytes = include_bytes!("fixtures/examples/A. versicolor CBS 583.65.gbk");

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
fn standalone_fuzzy_positions_support_bounds_strands_and_translation() {
    let bytes = br#"LOCUS       FUZZY                      9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             <1
                     /protein_id="RDA95433.1"
                     /translation="M"
     CDS             >9
                     /translation="A"
     CDS             complement(<1)
                     /translation="R"
     CDS             join(<1,2..9)
     CDS             complement(join(<1,2..9))
ORIGIN
        1 atggcttaa
//
"#;
    let cluster = parse_genbank("fuzzy.gbk", bytes).unwrap();
    let genes = &cluster.loci[0].genes;
    assert_eq!(genes.len(), 5);
    assert_eq!(genes[0].label, "RDA95433.1");
    assert_eq!(
        genes
            .iter()
            .map(|gene| (gene.start, gene.end, gene.strand))
            .collect::<Vec<_>>(),
        [(0, 1, 1), (8, 9, 1), (0, 1, -1), (0, 9, 1), (0, 9, -1)]
    );
    assert_eq!(genes[0].translation, "M");
    assert_eq!(genes[3].translation, "MA");
    assert_eq!(genes[4].translation, "LSH");

    let text = std::str::from_utf8(bytes).unwrap();
    for invalid in ["<0", "<10", "<ACCESSION", "OTHER"] {
        assert!(parse_genbank("invalid.gbk", text.replace("<1", invalid).as_bytes()).is_err());
    }

    let circular = include_str!("fixtures/circular.gbk").replace(
        "ORIGIN",
        "     CDS             complement(<1)\n                     /translation=\"M\"\nORIGIN",
    );
    let rotated = parse_genbank("circular.gbk", circular.as_bytes()).unwrap();
    let last = rotated.loci[0].genes.last().unwrap();
    assert_eq!((last.start, last.end, last.strand), (6, 7, -1));
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
fn translation_fallback_honours_frame_and_genetic_code_after_splicing() {
    let bytes = include_bytes!("fixtures/translation.gbk");
    let cluster = parse_genbank("framed.gbk", bytes).unwrap();
    assert_eq!(cluster.loci[0].genes[0].translation, "MW");
    assert_eq!(cluster.loci[0].genes[1].translation, "MW");
    assert_eq!(
        (cluster.loci[0].genes[1].start, cluster.loci[0].genes[1].end),
        (10, 20)
    );
}

#[test]
fn invalid_translation_qualifiers_fail_only_when_translation_is_derived() {
    for (frame, table) in [("0", "1"), ("4", "1"), ("bad", "1"), ("1", "999")] {
        let bytes = format!(
            "LOCUS       TEST                       9 bp    DNA     linear   UNA 01-JAN-2000\nFEATURES             Location/Qualifiers\n     CDS             1..9\n                     /locus_tag=\"example\"\n                     /codon_start={frame}\n                     /transl_table={table}\nORIGIN\n        1 atgtgataa\n//\n"
        );
        let error = parse_genbank("invalid.gbk", bytes.as_bytes()).unwrap_err();
        assert!(error.to_string().contains("translate CDS example"));
        let supplied = bytes.replace(
            "ORIGIN",
            "                     /translation=\"M W *\"\nORIGIN",
        );
        let cluster = parse_genbank("supplied.gbk", supplied.as_bytes()).unwrap();
        assert_eq!(cluster.loci[0].genes[0].translation, "MW");
    }
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
    let bytes = include_bytes!("fixtures/examples/A. burnettii MST-FP2249.gbk");
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
    let options = ParseOptions {
        as_separate_clusters: true,
        ..ParseOptions::default()
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
    let options = ParseOptions {
        ranges: BTreeMap::from([("RANGED".to_owned(), SequenceRange { start: 6, end: 12 })]),
        ..ParseOptions::default()
    };

    let clusters = parse_genbank_with_options("range.gbk", bytes, &options).unwrap();
    let locus = &clusters[0].loci[0];

    assert_eq!((locus.start, locus.end), (6, 12));
    assert_eq!(locus.genes.len(), 1);
    assert_eq!(locus.genes[0].label, "inside");
    assert_eq!((locus.genes[0].start, locus.genes[0].end), (6, 12));
    let invalid = String::from_utf8(bytes.to_vec()).unwrap().replace("ORIGIN", "     CDS             20..30\n                     /locus_tag=\"invalid\"\n                     /translation=\"M\"\nORIGIN");
    assert!(parse_genbank_with_options("range.gbk", invalid.as_bytes(), &options).is_err());
}

#[test]
fn circular_rotation_rejects_missing_sequence_and_invalid_locations() {
    use clinker_core::ParseError;

    let circular = std::str::from_utf8(include_bytes!("fixtures/circular.gbk")).unwrap();
    let contig = format!(
        "{}CONTIG      join(OTHER:1..30)\n//\n",
        circular.split_once("ORIGIN").unwrap().0
    );
    assert!(matches!(
        parse_genbank("contig.gbk", contig.as_bytes()),
        Err(ParseError::MissingRotationSequence { .. })
    ));
    let original = parse_genbank_with_options(
        "contig.gbk",
        contig.as_bytes(),
        &ParseOptions {
            set_origin: false,
            ..ParseOptions::default()
        },
    )
    .unwrap();
    assert_eq!(original[0].loci[0].genes[0].translation, "MMMMM");
    let whole_record = contig.replace("join(25..30,1..9)", "1..30");
    assert!(parse_genbank("whole-record.gbk", whole_record.as_bytes()).is_ok());
    let between = circular.replace("ORIGIN", "     misc_feature    30^1\nORIGIN");
    assert!(parse_genbank("between.gbk", between.as_bytes()).is_ok());

    let invalid = circular.replace(
        "ORIGIN",
        "     misc_feature    9223372036854775799..9223372036854775807\nORIGIN",
    );
    assert!(matches!(
        parse_genbank("invalid.gbk", invalid.as_bytes()),
        Err(ParseError::RotationLocation { .. })
    ));
}

#[test]
fn rotates_a_fixture_backed_circular_origin_crossing_feature() {
    let bytes = include_bytes!("fixtures/circular.gbk");

    let shifted =
        parse_genbank_with_options("circular.gbk", bytes, &ParseOptions::default()).unwrap();
    let unchanged = parse_genbank_with_options(
        "circular.gbk",
        bytes,
        &ParseOptions {
            set_origin: false,
            ..ParseOptions::default()
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

#[test]
fn explicit_circular_trim_uses_original_record_coordinates() {
    let options = ParseOptions {
        ranges: [("CIRCULAR".to_owned(), SequenceRange { start: 9, end: 18 })].into(),
        ..ParseOptions::default()
    };
    let clusters = parse_genbank_with_options(
        "circular.gbk",
        include_bytes!("fixtures/circular.gbk"),
        &options,
    )
    .unwrap();
    let locus = &clusters[0].loci[0];
    assert_eq!((locus.start, locus.end), (9, 18));
    assert_eq!(locus.genes.len(), 1);
    assert_eq!(locus.genes[0].label, "middle");
    assert_eq!((locus.genes[0].start, locus.genes[0].end), (9, 18));
}
