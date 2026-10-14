use std::{collections::BTreeMap, io::Cursor};

use gb_io::{
    reader::SeqReader,
    seq::{Feature, Location, LocationError, Seq},
};
use thiserror::Error;

use crate::{Cluster, Gene, Locus};

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("failed to parse GenBank input: {0}")]
    GenBank(#[from] gb_io::reader::GbParserError),
    #[error("GenBank record has no LOCUS name")]
    MissingRecordName,
    #[error("GenBank input contains no records")]
    MissingRecords,
    #[error("could not determine coordinates for CDS feature {label}: {source}")]
    FeatureLocation {
        label: String,
        #[source]
        source: LocationError,
    },
    #[error("range {start}..{end} for record {record} is outside its sequence length {length}")]
    InvalidRange {
        record: String,
        start: usize,
        end: usize,
        length: usize,
    },
}

/// A zero-based, end-exclusive interval within a GenBank record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SequenceRange {
    pub start: usize,
    pub end: usize,
}

/// Input parsing behavior shared by native and browser frontends.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseOptions {
    /// Record name to zero-based, end-exclusive extraction interval.
    pub ranges: BTreeMap<String, SequenceRange>,
    /// Emit each record as its own cluster instead of one cluster per file.
    pub as_separate_clusters: bool,
    /// Rotate circular records when a CDS or gene crosses the origin.
    pub set_origin: bool,
}

impl Default for ParseOptions {
    fn default() -> Self {
        Self {
            ranges: BTreeMap::new(),
            as_separate_clusters: false,
            set_origin: true,
        }
    }
}

/// Backwards-compatible name retained for callers that only parse GenBank.
pub type GenBankOptions = ParseOptions;

pub fn parse_genbank(file_name: &str, bytes: &[u8]) -> Result<Cluster, ParseError> {
    let mut clusters = parse_genbank_with_options(file_name, bytes, &ParseOptions::default())?;
    clusters.pop().ok_or(ParseError::MissingRecords)
}

/// Parse a GenBank file with explicit range, record grouping, and circular
/// origin behavior. A file normally produces one cluster, while
/// `as_separate_clusters` produces one cluster per record.
pub fn parse_genbank_with_options(
    file_name: &str,
    bytes: &[u8],
    options: &ParseOptions,
) -> Result<Vec<Cluster>, ParseError> {
    let mut loci = Vec::new();
    for record in SeqReader::new(Cursor::new(bytes)) {
        let mut record = record?;
        let name = record.name.clone().ok_or(ParseError::MissingRecordName)?;
        let organism = record
            .source
            .as_ref()
            .and_then(|source| {
                source
                    .organism
                    .as_deref()
                    .and_then(|organism| organism.lines().next())
                    .map(str::trim)
                    .filter(|organism| !organism.is_empty())
                    .map(str::to_owned)
            })
            .or_else(|| record.source.as_ref().map(|source| source.source.clone()));
        if options.set_origin
            && let Some(origin) = circular_origin_shift(&record)
        {
            record = record.set_origin(origin as i64);
        }
        let locus = locus_from_record(record, name.clone(), options.ranges.get(&name))?;
        loci.push((name, organism, locus));
    }

    if loci.is_empty() {
        return Err(ParseError::MissingRecords);
    }
    if options.as_separate_clusters {
        Ok(loci
            .into_iter()
            .map(|(record_name, organism, locus)| Cluster {
                name: organism.clone().unwrap_or(record_name),
                source_file: Some(file_name.to_owned()),
                organism,
                loci: vec![locus],
            })
            .collect())
    } else {
        let organism = common_organism(&loci);
        Ok(vec![Cluster {
            name: cluster_name(file_name),
            source_file: Some(file_name.to_owned()),
            organism,
            loci: loci.into_iter().map(|(_, _, locus)| locus).collect(),
        }])
    }
}

fn common_organism(loci: &[(String, Option<String>, Locus)]) -> Option<String> {
    let first = loci.first()?.1.as_ref()?;
    loci.iter()
        .all(|(_, organism, _)| organism.as_ref() == Some(first))
        .then(|| first.clone())
}

fn locus_from_record(
    record: Seq,
    name: String,
    range: Option<&SequenceRange>,
) -> Result<Locus, ParseError> {
    let length = record.len() as usize;
    let range = range.copied().unwrap_or(SequenceRange {
        start: 0,
        end: length,
    });
    if range.start >= range.end || range.end > length {
        return Err(ParseError::InvalidRange {
            record: name,
            start: range.start,
            end: range.end,
            length,
        });
    }

    let contained = |feature: &&Feature| {
        location_bounds(&feature.location)
            .is_ok_and(|(start, end)| start >= range.start as i64 && end <= range.end as i64)
    };
    let genes = record
        .features
        .iter()
        .filter(|feature| feature.kind.eq_ignore_ascii_case("CDS"))
        .filter(contained)
        .enumerate()
        .map(|(index, feature)| gene_from_feature(&record, feature, index))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Locus {
        name,
        start: range.start,
        end: range.end,
        genes,
    })
}

fn circular_origin_shift(record: &Seq) -> Option<usize> {
    if !record.is_circular() {
        return None;
    }
    let length = record.len();
    record
        .features
        .iter()
        .filter(|feature| {
            feature.kind.eq_ignore_ascii_case("gene") || feature.kind.eq_ignore_ascii_case("CDS")
        })
        .find_map(|feature| {
            let parts = location_parts(&feature.location);
            let touches_start = parts.iter().any(|&(start, _)| start == 0);
            let end_part = parts.iter().find(|&&(_, end)| end == length);
            (touches_start && end_part.is_some()).then(|| end_part.unwrap().0 as usize)
        })
}

fn location_parts(location: &Location) -> Vec<(i64, i64)> {
    match location {
        Location::Range((start, _), (end, _)) => vec![(*start, *end)],
        Location::Between(start, end) => vec![(*start, end + 1)],
        Location::Complement(inner) => location_parts(inner),
        Location::Join(parts)
        | Location::Order(parts)
        | Location::Bond(parts)
        | Location::OneOf(parts) => parts.iter().flat_map(location_parts).collect(),
        Location::External(_, Some(inner)) => location_parts(inner),
        Location::External(_, None) | Location::Gap(_) => Vec::new(),
    }
}

fn gene_from_feature(record: &Seq, feature: &Feature, index: usize) -> Result<Gene, ParseError> {
    let names = feature
        .qualifiers
        .iter()
        .filter_map(|(key, value)| {
            (key.as_ref() != "translation")
                .then(|| value.as_ref().map(|value| (key.to_string(), value.clone())))
                .flatten()
        })
        .collect::<Vec<_>>();
    let label = label_for(feature, index);
    let cds_bounds =
        location_bounds(&feature.location).map_err(|source| ParseError::FeatureLocation {
            label: label.clone(),
            source,
        })?;
    // A CDS may cover only part of its parent gene (for example when UTRs are
    // annotated). Protein glyphs represent the translated feature, so using
    // the enclosing gene can introduce apparent overlaps with neighbouring
    // CDS features. Keep their geometry tied to the CDS itself.
    let (start, end) = cds_bounds;
    let translation = if let Some(translation) = feature.qualifier_values("translation").next() {
        translation.split_whitespace().collect()
    } else {
        let coding_sequence = record
            .extract_location(&feature.location)
            .map_err(|source| ParseError::FeatureLocation {
                label: label.clone(),
                source,
            })?;
        translate_standard(&coding_sequence)
    };
    let strand = location_strand(&feature.location).unwrap_or(1);
    Ok(Gene {
        label,
        names,
        start: start as usize,
        end: end as usize,
        strand,
        translation: normalize_translation(translation),
    })
}

/// GenBank translations commonly include the terminal stop as `*`, but that
/// codon is not part of the translated protein sequence.
fn normalize_translation(translation: String) -> String {
    translation.trim_end_matches('*').to_owned()
}

fn location_bounds(location: &Location) -> Result<(i64, i64), LocationError> {
    match location {
        Location::Range((start, _), (end, _)) => Ok((*start, *end)),
        Location::Between(start, end) => Ok((*start, end + 1)),
        Location::Complement(inner) => location_bounds(inner),
        Location::Join(parts)
        | Location::Order(parts)
        | Location::Bond(parts)
        | Location::OneOf(parts) => bounds_for_parts(parts),
        Location::External(_, Some(inner)) => location_bounds(inner),
        _ => location.find_bounds(),
    }
}

fn bounds_for_parts(parts: &[Location]) -> Result<(i64, i64), LocationError> {
    let mut bounds = parts.iter().map(location_bounds);
    let (mut start, mut end) = bounds.next().ok_or(LocationError::Empty)??;
    for part in bounds {
        let (part_start, part_end) = part?;
        start = start.min(part_start);
        end = end.max(part_end);
    }
    Ok((start, end))
}

fn location_strand(location: &Location) -> Option<i8> {
    fn visit(location: &Location, orientation: i8) -> Option<i8> {
        match location {
            Location::Range(..) | Location::Between(..) => Some(orientation),
            Location::Complement(inner) => visit(inner, -orientation),
            Location::Join(parts)
            | Location::Order(parts)
            | Location::Bond(parts)
            | Location::OneOf(parts) => {
                let mut strands = parts.iter().filter_map(|part| visit(part, orientation));
                let strand = strands.next()?;
                strands.all(|other| other == strand).then_some(strand)
            }
            Location::External(_, Some(inner)) => visit(inner, orientation),
            Location::External(_, None) | Location::Gap(_) => None,
        }
    }
    visit(location, 1)
}

/// Translate an extracted CDS using the standard genetic code.
fn translate_standard(coding_sequence: &[u8]) -> String {
    coding_sequence
        .as_chunks::<3>()
        .0
        .iter()
        .map(|codon| amino_acid(codon).unwrap_or('X'))
        .collect()
}

fn amino_acid(codon: &[u8]) -> Option<char> {
    const TABLE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let first = base_index(codon[0])?;
    let second = base_index(codon[1])?;
    let third = base_index(codon[2])?;
    Some(TABLE[first * 16 + second * 4 + third] as char)
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

fn label_for(feature: &Feature, index: usize) -> String {
    const LABEL_KEYS: [&str; 7] = [
        "protein_id",
        "locus_tag",
        "id",
        "ID",
        "gene",
        "label",
        "name",
    ];

    LABEL_KEYS
        .iter()
        .find_map(|key| feature.qualifier_values(key).next())
        .map_or_else(|| format!("CDS_{index}"), ToOwned::to_owned)
}

fn cluster_name(file_name: &str) -> String {
    file_name
        .rsplit_once('.')
        .map_or_else(|| file_name.to_owned(), |(stem, _)| stem.to_owned())
}
