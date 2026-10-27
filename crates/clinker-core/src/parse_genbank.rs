use std::{collections::BTreeMap, io::Cursor};

use gb_io::{
    reader::SeqReader,
    seq::{Feature, Location, LocationError, Seq},
};
use thiserror::Error;

use crate::protein::{normalise_translation, translate_cds};
use crate::{Cluster, Gene, Locus, MAX_SEQUENCE_COORDINATE};

#[derive(Debug, Error)]
pub enum ParseError {
    #[error("failed to parse GenBank input: {0}")]
    GenBank(#[from] gb_io::reader::GbParserError),
    #[error("GenBank record has no LOCUS name")]
    MissingRecordName,
    #[error("GenBank input contains no records")]
    MissingRecords,
    #[error("record {record} has invalid sequence length {length}")]
    InvalidRecordLength { record: String, length: i64 },
    #[error("could not determine coordinates for CDS feature {label}: {source}")]
    FeatureLocation {
        label: String,
        #[source]
        source: LocationError,
    },
    #[error("could not translate CDS {label}: {message}")]
    Translation { label: String, message: String },
    #[error("range {start}..{end} for record {record} is outside its sequence length {length}")]
    InvalidRange {
        record: String,
        start: u64,
        end: u64,
        length: u64,
    },
}

/// A zero-based, end-exclusive interval within a GenBank record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SequenceRange {
    pub start: u64,
    pub end: u64,
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
        // Explicit trim coordinates refer to the original record, before rotation.
        if options.set_origin
            && !options.ranges.contains_key(&name)
            && let Some(origin) = circular_origin_shift(&record)
        {
            record = record.set_origin(origin);
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
    let length = u64::try_from(record.len())
        .ok()
        .filter(|length| *length > 0 && *length <= MAX_SEQUENCE_COORDINATE)
        .ok_or_else(|| ParseError::InvalidRecordLength {
            record: name.clone(),
            length: record.len(),
        })?;
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

    let mut genes = Vec::new();
    for feature in record
        .features
        .iter()
        .filter(|feature| feature.kind.eq_ignore_ascii_case("CDS"))
    {
        let index = genes.len();
        let bounds = location_interval(&feature.location, length).map_err(|source| {
            ParseError::FeatureLocation {
                label: label_for(feature, index),
                source,
            }
        })?;
        if bounds.0 >= range.start && bounds.1 <= range.end {
            genes.push(gene_from_feature(&record, feature, index, bounds)?);
        }
    }
    Ok(Locus {
        name,
        start: range.start,
        end: range.end,
        genes,
    })
}

fn circular_origin_shift(record: &Seq) -> Option<i64> {
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
            let end_part = parts
                .iter()
                .find(|&&(start, end)| end == length && start >= 0 && start < length);
            (touches_start && end_part.is_some()).then(|| end_part.unwrap().0)
        })
}

fn location_parts(location: &Location) -> Vec<(i64, i64)> {
    match location {
        Location::Range((start, _), (end, _)) => vec![(*start, *end)],
        Location::Between(..) => location_bounds(location).into_iter().collect(),
        Location::Complement(inner) => location_parts(inner),
        Location::Join(parts)
        | Location::Order(parts)
        | Location::Bond(parts)
        | Location::OneOf(parts) => parts.iter().flat_map(location_parts).collect(),
        Location::External(_, Some(inner)) => location_parts(inner),
        Location::External(_, None) | Location::Gap(_) => Vec::new(),
    }
}

fn gene_from_feature(
    record: &Seq,
    feature: &Feature,
    index: usize,
    (start, end): (u64, u64),
) -> Result<Gene, ParseError> {
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
    // A CDS may cover only part of its parent gene (for example when UTRs are
    // annotated). Protein glyphs represent the translated feature, so using
    // the enclosing gene can introduce apparent overlaps with neighbouring
    // CDS features. Keep their geometry tied to the CDS itself.
    let translation = if let Some(translation) = feature.qualifier_values("translation").next() {
        normalise_translation(translation)
    } else {
        let coding_sequence = record
            .extract_location(&feature.location)
            .map_err(|source| ParseError::FeatureLocation {
                label: label.clone(),
                source,
            })?;
        let codon_start = feature
            .qualifier_values("codon_start")
            .next()
            .unwrap_or("1");
        let offset = match codon_start.trim() {
            "1" => 0,
            "2" => 1,
            "3" => 2,
            _ => {
                return Err(ParseError::Translation {
                    label: label.clone(),
                    message: format!("codon_start must be 1, 2 or 3, got {codon_start:?}"),
                });
            }
        };
        translate_cds(
            coding_sequence.get(offset..).unwrap_or_default(),
            feature.qualifier_values("transl_table").next(),
        )
        .map_err(|message| ParseError::Translation {
            label: label.clone(),
            message,
        })?
    };
    let strand = location_strand(&feature.location).unwrap_or(1);
    Ok(Gene {
        label,
        names,
        start,
        end,
        strand,
        translation,
    })
}

fn location_interval(location: &Location, length: u64) -> Result<(u64, u64), LocationError> {
    let (start, end) = location_bounds(location)?;
    if let (Ok(start), Ok(end)) = (u64::try_from(start), u64::try_from(end))
        && start < end
        && end <= length
    {
        return Ok((start, end));
    }
    Err(LocationError::OutOfBounds(location.clone()))
}

fn location_bounds(location: &Location) -> Result<(i64, i64), LocationError> {
    match location {
        Location::Range((start, _), (end, _)) => Ok((*start, *end)),
        Location::Between(start, end) => end
            .checked_add(1)
            .map(|end| (*start, end))
            .ok_or_else(|| LocationError::OutOfBounds(location.clone())),
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
