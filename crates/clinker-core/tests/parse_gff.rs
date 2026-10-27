use std::collections::BTreeMap;

use clinker_core::{
    GenBankOptions, GffError, InputFile, SequenceRange, parse_gff3, parse_input_files_with_options,
};

const GFF: &[u8] = include_bytes!("fixtures/features.gff3");
const FASTA: &[u8] = include_bytes!("fixtures/features.fna");

#[test]
fn parses_multirecord_compound_cds_coordinates() {
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
fn overlapping_parent_genes_do_not_expand_cds_arrows_or_exclude_trimmed_cds() {
    let input = b"##gff-version 3\n\
seq\t.\tgene\t1\t15\t.\t+\t.\tID=gene_a\n\
seq\t.\tmRNA\t1\t15\t.\t+\t.\tID=tx_a;Parent=gene_a\n\
seq\t.\tCDS\t4\t9\t.\t+\t0\tID=cds_a;Parent=tx_a\n\
seq\t.\tgene\t7\t24\t.\t-\t.\tID=gene_b\n\
seq\t.\tmRNA\t7\t24\t.\t-\t.\tID=tx_b;Parent=gene_b\n\
seq\t.\tCDS\t13\t18\t.\t-\t0\tID=cds_b;Parent=tx_b\n\
##FASTA\n>seq\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    let cluster = parse_gff3("utrs.gff3", input, None, &GenBankOptions::default()).unwrap();
    let genes = &cluster.loci[0].genes;
    assert_eq!(genes.len(), 2);
    assert_eq!((genes[0].start, genes[0].end, genes[0].strand), (3, 9, 1));
    assert_eq!(
        (genes[1].start, genes[1].end, genes[1].strand),
        (12, 18, -1)
    );
    assert!(genes[0].end <= genes[1].start);

    let options = GenBankOptions {
        ranges: BTreeMap::from([("seq".into(), SequenceRange { start: 3, end: 9 })]),
        ..GenBankOptions::default()
    };
    let trimmed = parse_gff3("utrs.gff3", input, None, &options).unwrap();
    assert_eq!(trimmed.loci[0].genes.len(), 1);
    // GFF attributes are unordered; independent parses may iterate them differently.
    let mut trimmed_gene = trimmed.loci[0].genes[0].clone();
    let mut original_gene = genes[0].clone();
    trimmed_gene.names.sort();
    original_gene.names.sort();
    assert_eq!(trimmed_gene, original_gene);
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
fn initial_phase_can_span_short_cds_segments_on_either_strand() {
    for (strand, sequence, first, second) in [
        ("+", "CCATGGCTTAA", (1, 1, 2), (2, 11, 1)),
        ("-", "TTAAGCCATGG", (11, 11, 2), (1, 10, 1)),
    ] {
        let input = format!(
            "##gff-version 3\nseq\t.\tCDS\t{}\t{}\t.\t{strand}\t{}\tID=split\nseq\t.\tCDS\t{}\t{}\t.\t{strand}\t{}\tID=split\n##FASTA\n>seq\n{sequence}\n",
            first.0, first.1, first.2, second.0, second.1, second.2
        );
        let cluster = parse_gff3(
            "partial.gff3",
            input.as_bytes(),
            None,
            &GenBankOptions::default(),
        )
        .unwrap();
        assert_eq!(
            cluster.loci[0].genes[0].translation, "MA",
            "strand {strand}"
        );
    }
}

#[test]
fn genetic_code_applies_only_to_derived_translation() {
    let input = "##gff-version 3\nseq\t.\tCDS\t1\t9\t.\t+\t0\tID=example;transl_table=4\n##FASTA\n>seq\nATGTGATAA\n";
    let cluster = parse_gff3(
        "code.gff3",
        input.as_bytes(),
        None,
        &GenBankOptions::default(),
    )
    .unwrap();
    assert_eq!(cluster.loci[0].genes[0].translation, "MW");
    let invalid = input.replace("transl_table=4", "transl_table=999");
    assert!(matches!(
        parse_gff3(
            "invalid.gff3",
            invalid.as_bytes(),
            None,
            &GenBankOptions::default()
        ),
        Err(GffError::Translation { .. })
    ));
    let supplied = invalid.replace("transl_table=999", "transl_table=999;translation=MW*");
    let cluster = parse_gff3(
        "supplied.gff3",
        supplied.as_bytes(),
        None,
        &GenBankOptions::default(),
    )
    .unwrap();
    let gene = &cluster.loci[0].genes[0];
    assert_eq!(gene.translation, "MW");
    assert!(gene.names.iter().all(|(key, _)| key != "translation"));
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
fn large_absolute_coordinates_survive_parsing_arrangement_and_project_round_trips() {
    use clinker_core::{
        Analysis, AnalysisSession, ClinkerProject, GeneRef, Link, MAX_SEQUENCE_COORDINATE,
    };
    for origin in [
        3_000_000_000_u64,
        5_000_000_000,
        MAX_SEQUENCE_COORDINATE - 18,
    ] {
        let annotation = format!(
            "##gff-version 3\n##sequence-region seq {} {}\nseq\t.\tCDS\t{}\t{}\t.\t+\t0\tID=first\nseq\t.\tCDS\t{}\t{}\t.\t+\t0\tID=second\n",
            origin + 1,
            origin + 18,
            origin + 1,
            origin + 9,
            origin + 10,
            origin + 18
        );
        let clusters = [
            b">seq\nATGGCTTAAATGAAATAA\n".as_slice(),
            b">seq\nATGAAATAAATGGCTTAA\n",
        ]
        .into_iter()
        .map(|fasta| {
            parse_gff3(
                "large.gff3",
                annotation.as_bytes(),
                Some(fasta),
                &GenBankOptions::default(),
            )
            .unwrap()
        })
        .collect();
        let reference = |cluster, gene| GeneRef {
            cluster,
            locus: 0,
            gene,
        };
        let analysis = Analysis {
            clusters,
            links: vec![
                Link {
                    query: reference(0, 0),
                    target: reference(1, 1),
                    identity: 1.0,
                    similarity: 1.0,
                },
                Link {
                    query: reference(0, 1),
                    target: reference(1, 0),
                    identity: 1.0,
                    similarity: 1.0,
                },
            ],
        };
        let data = analysis.to_auto_arranged_plot_data(&[0, 1]).unwrap();
        let reversed = &data.clusters[1].loci[0];
        assert_eq!(reversed.offset, Some(0.0));
        assert_eq!(
            (
                reversed.genes[0].start,
                reversed.genes[0].end,
                reversed.genes[0].strand
            ),
            (origin, origin + 9, -1)
        );
        let session = AnalysisSession::new(analysis, 0.3, Vec::new()).unwrap();
        let project = ClinkerProject::new(session, &data).unwrap();
        assert_eq!(
            ClinkerProject::from_json(&project.to_json().unwrap()).unwrap(),
            project
        );
    }
    for end in [MAX_SEQUENCE_COORDINATE + 1, u64::MAX] {
        let annotation = format!("##sequence-region seq 1 {end}\n");
        assert!(
            parse_gff3(
                "large.gff3",
                annotation.as_bytes(),
                Some(b">seq\nATGGCTTAA\n"),
                &GenBankOptions::default()
            )
            .is_err()
        );
    }
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
fn rejects_invalid_gff_coordinates_without_panicking_or_silently_dropping_the_gene() {
    for (start, end) in [(0, 9), (1, 12), (8, 7)] {
        let input = format!("##gff-version 3\nseq\t.\tCDS\t{start}\t{end}\t.\t+\t0\tID=invalid\n");
        let result = parse_gff3(
            "invalid.gff3",
            input.as_bytes(),
            Some(b">seq\nATGGCTTAA\n"),
            &GenBankOptions::default(),
        );
        assert!(matches!(result, Err(GffError::InvalidCoordinates { .. })));
    }
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
