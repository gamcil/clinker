use std::{
    collections::{BTreeMap, BTreeSet},
    io::Cursor,
};

use bio::{
    alphabets::dna,
    io::{
        fasta,
        gff::{self, GffType},
    },
};
use thiserror::Error;

use crate::protein::{normalise_translation, translate_standard};
use crate::{Cluster, Gene, Locus, ParseOptions, SequenceRange};

#[derive(Debug, Error)]
pub enum GffError {
    #[error("failed to parse GFF3 input: {0}")]
    Gff(String),
    #[error("failed to parse FASTA input: {0}")]
    Fasta(String),
    #[error("GFF3 input has no companion or embedded FASTA sequence")]
    MissingFasta,
    #[error("FASTA input contains no records")]
    MissingFastaRecords,
    #[error("GFF3 record {record} has no matching FASTA sequence")]
    MissingSequence { record: String },
    #[error("CDS {label} has invalid coordinates {start}..{end} for record {record}")]
    InvalidCoordinates {
        record: String,
        label: String,
        start: usize,
        end: usize,
    },
    #[error(
        "range {start}..{end} for record {record} is outside its available sequence coordinates {available_start}..{available_end}"
    )]
    InvalidRange {
        record: String,
        start: usize,
        end: usize,
        available_start: usize,
        available_end: usize,
    },
    #[error("CDS {label} contains segments on different strands")]
    MixedStrands { label: String },
}

#[derive(Debug, Clone)]
struct FeatureGroup {
    label: String,
    records: Vec<gff::Record>,
}

#[derive(Debug, Clone, Copy)]
struct Region {
    /// Zero-based coordinate represented by the first base in the FASTA record.
    sequence_origin: usize,
    /// Zero-based, end-exclusive coordinates selected for this locus.
    selected: SequenceRange,
}

/// Parse one GFF3 annotation and its companion FASTA. When `fasta_bytes` is
/// absent, an embedded `##FASTA` section is used.
pub fn parse_gff3(
    file_name: &str,
    gff_bytes: &[u8],
    fasta_bytes: Option<&[u8]>,
    options: &ParseOptions,
) -> Result<Cluster, GffError> {
    let (annotation, embedded_fasta) = split_embedded_fasta(gff_bytes);
    let fasta_bytes = fasta_bytes
        .or(embedded_fasta)
        .ok_or(GffError::MissingFasta)?;
    let sequences = parse_fasta(fasta_bytes)?;
    let directives = sequence_regions(annotation);
    let features = parse_features(annotation)?;
    let mut loci = Vec::new();
    for (record_name, sequence) in &sequences {
        let record_features = features.get(record_name);
        let region = region_for(
            record_name,
            sequence.len(),
            directives.get(record_name).copied(),
            options.ranges.get(record_name).copied(),
        )?;
        let mut locus_genes = Vec::new();
        if let Some(groups) = record_features {
            for group in groups {
                if let Some(gene) = gene_from_group(record_name, sequence, region, group)? {
                    locus_genes.push(gene);
                }
            }
        }
        locus_genes.sort_by_key(|gene| (gene.start, gene.end, gene.label.clone()));
        loci.push(Locus {
            name: record_name.clone(),
            start: region.selected.start,
            end: region.selected.end,
            genes: locus_genes,
        });
    }

    for record_name in features.keys() {
        if !sequences.contains_key(record_name) {
            return Err(GffError::MissingSequence {
                record: record_name.clone(),
            });
        }
    }

    Ok(Cluster {
        name: cluster_name(file_name),
        source_file: Some(file_name.to_owned()),
        organism: None,
        loci,
    })
}

type FeatureMap = BTreeMap<String, Vec<FeatureGroup>>;

fn parse_features(annotation: &[u8]) -> Result<FeatureMap, GffError> {
    let mut reader = gff::Reader::new(Cursor::new(annotation), GffType::GFF3);
    let mut grouped = BTreeMap::<(String, String), Vec<gff::Record>>::new();
    let mut order = Vec::<(String, String)>::new();

    for (index, result) in reader.records().enumerate() {
        let record = result.map_err(|error| GffError::Gff(error.to_string()))?;
        let identifier = attribute(&record, "ID");
        if !record.feature_type().eq_ignore_ascii_case("CDS") {
            continue;
        }
        let key = identifier
            .or_else(|| attribute(&record, "Parent"))
            .unwrap_or_else(|| format!("CDS_{index}"));
        let compound_key = (record.seqname().to_owned(), key);
        if !grouped.contains_key(&compound_key) {
            order.push(compound_key.clone());
        }
        grouped.entry(compound_key).or_default().push(record);
    }

    let mut features = FeatureMap::new();
    for (record, key) in order {
        let records = grouped
            .remove(&(record.clone(), key.clone()))
            .unwrap_or_default();
        let label = records
            .first()
            .map(|record| label_for(record, &key))
            .unwrap_or(key);
        features
            .entry(record)
            .or_default()
            .push(FeatureGroup { label, records });
    }
    Ok(features)
}

fn gene_from_group(
    record_name: &str,
    sequence: &[u8],
    region: Region,
    group: &FeatureGroup,
) -> Result<Option<Gene>, GffError> {
    let raw_cds_start = group
        .records
        .iter()
        .map(|record| *record.start() as usize)
        .min()
        .unwrap_or(0);
    let cds_end = group
        .records
        .iter()
        .map(|record| *record.end() as usize)
        .max()
        .unwrap_or(0);
    if raw_cds_start == 0 || raw_cds_start > cds_end {
        return Err(GffError::InvalidCoordinates {
            record: record_name.to_owned(),
            label: group.label.clone(),
            start: raw_cds_start,
            end: cds_end,
        });
    }
    let cds_start = raw_cds_start - 1;
    if cds_start < region.selected.start || cds_end > region.selected.end {
        return Ok(None);
    }

    let strands = group
        .records
        .iter()
        .filter_map(|record| record.strand())
        .map(|strand| strand.strand_symbol().to_owned())
        .collect::<BTreeSet<_>>();
    if strands.len() > 1 {
        return Err(GffError::MixedStrands {
            label: group.label.clone(),
        });
    }
    let strand = if strands.first().is_some_and(|strand| strand == "-") {
        -1
    } else {
        1
    };

    let mut ordered = group.records.iter().collect::<Vec<_>>();
    ordered.sort_by_key(|record| *record.start());
    if strand == -1 {
        ordered.reverse();
    }
    let mut coding_sequence = Vec::new();
    for (segment_index, record) in ordered.into_iter().enumerate() {
        let raw_start = *record.start() as usize;
        let end = *record.end() as usize;
        if raw_start == 0 || raw_start > end {
            return Err(GffError::InvalidCoordinates {
                record: record_name.to_owned(),
                label: group.label.clone(),
                start: raw_start,
                end,
            });
        }
        let start = raw_start - 1;
        if start < region.sequence_origin
            || end > region.sequence_origin + sequence.len()
            || start >= end
        {
            return Err(GffError::InvalidCoordinates {
                record: record_name.to_owned(),
                label: group.label.clone(),
                start,
                end,
            });
        }
        let mut part =
            sequence[start - region.sequence_origin..end - region.sequence_origin].to_vec();
        if strand == -1 {
            part = dna::revcomp(part);
        }
        // Only the first biological segment's phase trims an incomplete codon
        // at the edge of the annotated CDS. Later phases describe the frame at
        // exon boundaries; trimming those bases would delete split codons.
        let phase: Option<u8> = record.phase().clone().try_into().unwrap_or(None);
        let skip = if segment_index == 0 {
            usize::from(phase.unwrap_or(0))
        } else {
            0
        };
        if skip < part.len() {
            coding_sequence.extend_from_slice(&part[skip..]);
        }
    }

    let first = &group.records[0];
    let names = first
        .attributes()
        .iter_all()
        .flat_map(|(key, values)| values.iter().map(move |value| (key.clone(), value.clone())))
        .collect::<Vec<_>>();
    // Protein arrows represent the CDS span, as in the GenBank parser. Parent
    // genes can include UTRs that overlap neighbouring coding features.
    Ok(Some(Gene {
        label: group.label.clone(),
        names,
        start: cds_start,
        end: cds_end,
        strand,
        translation: attribute(first, "translation")
            .map(|translation| normalise_translation(&translation))
            .unwrap_or_else(|| translate_standard(&coding_sequence)),
    }))
}

fn parse_fasta(bytes: &[u8]) -> Result<BTreeMap<String, Vec<u8>>, GffError> {
    let mut sequences = BTreeMap::new();
    for result in fasta::Reader::new(Cursor::new(bytes)).records() {
        let record = result.map_err(|error| GffError::Fasta(error.to_string()))?;
        sequences.insert(record.id().to_owned(), record.seq().to_vec());
    }
    if sequences.is_empty() {
        Err(GffError::MissingFastaRecords)
    } else {
        Ok(sequences)
    }
}

fn region_for(
    record: &str,
    sequence_length: usize,
    directive: Option<(usize, usize)>,
    selected: Option<SequenceRange>,
) -> Result<Region, GffError> {
    let sequence_origin = directive
        .filter(|(start, end)| end - start + 1 == sequence_length)
        .map_or(0, |(start, _)| start - 1);
    let available = SequenceRange {
        start: sequence_origin,
        end: sequence_origin + sequence_length,
    };
    let directive_range = directive.map(|(start, end)| SequenceRange {
        start: start - 1,
        end,
    });
    let selected = selected.or(directive_range).unwrap_or(available);
    if selected.start < available.start
        || selected.start >= selected.end
        || selected.end > available.end
    {
        return Err(GffError::InvalidRange {
            record: record.to_owned(),
            start: selected.start,
            end: selected.end,
            available_start: available.start,
            available_end: available.end,
        });
    }
    Ok(Region {
        sequence_origin,
        selected,
    })
}

fn split_embedded_fasta(bytes: &[u8]) -> (&[u8], Option<&[u8]>) {
    let marker = b"##FASTA";
    let Some(start) = bytes
        .windows(marker.len())
        .position(|window| window == marker)
    else {
        return (bytes, None);
    };
    let fasta_start = bytes[start..]
        .iter()
        .position(|byte| *byte == b'\n')
        .map_or(bytes.len(), |offset| start + offset + 1);
    (&bytes[..start], Some(&bytes[fasta_start..]))
}

fn sequence_regions(bytes: &[u8]) -> BTreeMap<String, (usize, usize)> {
    String::from_utf8_lossy(bytes)
        .lines()
        .filter_map(|line| line.strip_prefix("##sequence-region"))
        .filter_map(|value| {
            let mut fields = value.split_whitespace();
            let name = fields.next()?.to_owned();
            let start = fields.next()?.parse().ok()?;
            let end = fields.next()?.parse().ok()?;
            (start > 0 && end >= start).then_some((name, (start, end)))
        })
        .collect()
}

fn attribute(record: &gff::Record, key: &str) -> Option<String> {
    record.attributes().get(key).cloned()
}

fn label_for(record: &gff::Record, fallback: &str) -> String {
    ["protein_id", "locus_tag", "ID", "gene", "Name", "name"]
        .iter()
        .find_map(|key| attribute(record, key))
        .unwrap_or_else(|| fallback.to_owned())
}

fn cluster_name(file_name: &str) -> String {
    file_name
        .rsplit_once('.')
        .map_or_else(|| file_name.to_owned(), |(stem, _)| stem.to_owned())
}
