use std::io::Cursor;

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
    #[error("could not determine coordinates for CDS feature {label}: {source}")]
    FeatureLocation {
        label: String,
        #[source]
        source: LocationError,
    },
}

/// Parse all records in one GenBank file into a clinker cluster.
///
pub fn parse_genbank(file_name: &str, bytes: &[u8]) -> Result<Cluster, ParseError> {
    let mut loci = Vec::new();

    for record in SeqReader::new(Cursor::new(bytes)) {
        let record = record?;
        let end = record.len() as usize;
        let name = record.name.clone().ok_or(ParseError::MissingRecordName)?;
        let genes = record
            .features
            .iter()
            .filter(|feature| feature.kind.eq_ignore_ascii_case("CDS"))
            .enumerate()
            .map(|(index, feature)| gene_from_feature(&record, feature, index))
            .collect::<Result<Vec<_>, _>>()?;

        loci.push(Locus {
            name,
            start: 0,
            end,
            genes,
        });
    }

    Ok(Cluster {
        name: cluster_name(file_name),
        loci,
    })
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
    let (start, end) =
        location_bounds(&feature.location).map_err(|source| ParseError::FeatureLocation {
            label: label.clone(),
            source,
        })?;
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
        translation,
    })
}

/// Return the enclosing interval for a GenBank location.
///
/// `Location::find_bounds` follows the order of `Join` members. That is wrong
/// for valid reverse-strand spliced CDS written as
/// `join(complement(high..high), complement(low..low))`: its first start is
/// greater than its last end. Plot coordinates need the genomic envelope, so
/// compound locations use the minimum start and maximum end instead.
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

/// Determine a feature's orientation even when each exon carries its own
/// `complement`, as seen in several GenBank submissions.
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
///
/// Unknown codons become `X`, and trailing incomplete codons are ignored. This
/// matches the useful behavior needed for drawing protein-homology links while
/// keeping translation independent of Python or a browser runtime.
fn translate_standard(coding_sequence: &[u8]) -> String {
    coding_sequence
        .chunks_exact(3)
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
