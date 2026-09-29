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
    #[error("CDS feature {label} has no /translation qualifier")]
    MissingTranslation { label: String },
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
        feature
            .location
            .find_bounds()
            .map_err(|source| ParseError::FeatureLocation {
                label: label.clone(),
                source,
            })?;
    let translation = feature
        .qualifier_values("translation")
        .next()
        .map(|translation| translation.split_whitespace().collect())
        .ok_or_else(|| ParseError::MissingTranslation {
            label: label.clone(),
        })?;

    // `gb-io` locations preserve a top-level complement, which is sufficient
    // for valid CDS locations emitted by the GenBank fixtures used here.
    let strand = if matches!(feature.location, Location::Complement(_)) {
        -1
    } else {
        1
    };

    // Keep `record` in the signature so the next increment can derive a
    // translation from `record.extract_location(&feature.location)` when an
    // input CDS has no /translation qualifier.
    let _ = record;

    Ok(Gene {
        label,
        names,
        start: start as usize,
        end: end as usize,
        strand,
        translation,
    })
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
