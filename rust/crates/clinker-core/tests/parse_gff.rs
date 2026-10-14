use std::collections::BTreeMap;

use clinker_core::{
    GenBankOptions, GffError, InputFile, SequenceRange, parse_gff3, parse_input_files_with_options,
};

const GFF: &[u8] = include_bytes!("fixtures/features.gff3");
const FASTA: &[u8] = include_bytes!("fixtures/features.fna");

#[test]
fn parses_multirecord_compound_cds_and_parent_gene_coordinates() {
    let cluster = parse_gff3(
        "features.gff3",
        GFF,
        Some(FASTA),
        &GenBankOptions::default(),
    )
    .expect("valid GFF3 and FASTA");

    assert_eq!(cluster.name, "features");
    assert_eq!(cluster.source_file.as_deref(), Some("features.gff3"));
    assert_eq!(cluster.loci.len(), 2);
    assert_eq!((cluster.loci[0].start, cluster.loci[0].end), (100, 130));

    let forward = &cluster.loci[0].genes[0];
    assert_eq!(forward.label, "PROT_A");
    assert_eq!((forward.start, forward.end, forward.strand), (100, 121, 1));
    assert_eq!(forward.translation, "MAKEF");

    let reverse = &cluster.loci[1].genes[0];
    assert_eq!(reverse.label, "PROT_B");
    assert_eq!((reverse.start, reverse.end, reverse.strand), (0, 18, -1));
    assert_eq!(reverse.translation, "MAKE");

    let orphan = &cluster.loci[1].genes[1];
    assert_eq!(orphan.label, "ORPHAN");
    assert_eq!((orphan.start, orphan.end), (18, 24));
    assert_eq!(orphan.translation, "MK");
}

#[test]
fn parses_an_embedded_fasta_section() {
    let input =
        b"##gff-version 3\nseq\t.\tCDS\t1\t9\t.\t+\t0\tID=embedded\n##FASTA\n>seq\nATGGCTTAA\n";
    let cluster = parse_gff3("embedded.gff3", input, None, &GenBankOptions::default())
        .expect("embedded FASTA is accepted");

    assert_eq!(cluster.loci[0].genes[0].translation, "MA");
}

#[test]
fn keeps_split_codons_at_compound_cds_boundaries() {
    let input = b"##gff-version 3\nseq\t.\tCDS\t1\t4\t.\t+\t0\tID=split\nseq\t.\tCDS\t5\t9\t.\t+\t2\tID=split\n##FASTA\n>seq\nATGGCTTAA\n";
    let cluster = parse_gff3("split.gff3", input, None, &GenBankOptions::default()).unwrap();

    assert_eq!(cluster.loci[0].genes[0].translation, "MA");
}

#[test]
fn decodes_gff3_attributes_with_rust_bio() {
    let input = b"##gff-version 3\nseq\t.\tCDS\t1\t9\t.\t+\t0\tID=encoded;protein_id=alpha%20beta\n##FASTA\n>seq\nATGGCTTAA\n";
    let cluster = parse_gff3("encoded.gff3", input, None, &GenBankOptions::default()).unwrap();

    assert_eq!(cluster.loci[0].genes[0].label, "alpha beta");
}

#[test]
fn accepts_fasta_records_without_cds_annotations() {
    let cluster = parse_gff3(
        "empty.gff3",
        b"##gff-version 3\n",
        Some(b">seq\nATGGCTTAA\n"),
        &GenBankOptions::default(),
    )
    .unwrap();

    assert!(cluster.loci[0].genes.is_empty());
}

#[test]
fn applies_absolute_ranges_to_offset_sequences() {
    let options = GenBankOptions {
        ranges: BTreeMap::from([(
            "chr1".to_owned(),
            SequenceRange {
                start: 110,
                end: 130,
            },
        )]),
        ..GenBankOptions::default()
    };
    let cluster = parse_gff3("features.gff3", GFF, Some(FASTA), &options).unwrap();

    assert_eq!((cluster.loci[0].start, cluster.loci[0].end), (110, 130));
    assert!(cluster.loci[0].genes.is_empty());
}

#[test]
fn uses_sequence_region_as_the_default_range_for_a_full_fasta_record() {
    let input = b"##gff-version 3\n##sequence-region seq 10 18\nseq\t.\tCDS\t10\t18\t.\t+\t0\tID=selected\n";
    let fasta = b">seq\nCCCCCCCCCATGGCTTAACCCCCCCCCCCC\n";
    let cluster = parse_gff3(
        "region.gff3",
        input,
        Some(fasta),
        &GenBankOptions::default(),
    )
    .unwrap();

    assert_eq!((cluster.loci[0].start, cluster.loci[0].end), (9, 18));
    assert_eq!(cluster.loci[0].genes[0].translation, "MA");
}

#[test]
fn rejects_zero_based_gff_coordinates_without_panicking() {
    let input = b"##gff-version 3\nseq\t.\tCDS\t0\t9\t.\t+\t0\tID=invalid\n";
    let result = parse_gff3(
        "invalid.gff3",
        input,
        Some(b">seq\nATGGCTTAA\n"),
        &GenBankOptions::default(),
    );

    assert!(matches!(result, Err(GffError::InvalidCoordinates { .. })));
}

#[test]
fn reports_missing_and_mismatched_sequence_data() {
    assert!(matches!(
        parse_gff3("features.gff3", GFF, None, &GenBankOptions::default()),
        Err(GffError::MissingFasta)
    ));
    assert!(matches!(
        parse_gff3(
            "features.gff3",
            GFF,
            Some(b">chr1\nATGGCTAAACCCGAATTTTAAGGGCCCAAA\n"),
            &GenBankOptions::default()
        ),
        Err(GffError::MissingSequence { record }) if record == "chr2"
    ));
}

#[test]
fn pairs_gff_and_fasta_inputs_by_stem() {
    let clusters = parse_input_files_with_options(
        &[
            InputFile {
                name: "features.gff3",
                bytes: GFF,
            },
            InputFile {
                name: "features.fna",
                bytes: FASTA,
            },
        ],
        &GenBankOptions::default(),
    )
    .unwrap();

    assert_eq!(clusters.len(), 1);
    assert_eq!(clusters[0].loci.len(), 2);
}
