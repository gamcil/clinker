use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{Analysis, GeneGroupDefinition, GeneRef};

pub const ANALYSIS_SESSION_FORMAT: &str = "clinker-analysis";
pub const ANALYSIS_SESSION_VERSION: u32 = 1;

/// Versioned analysis state used to resume work or add more input clusters.
///
/// Unlike plot data, a session deliberately retains protein translations.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AnalysisSession {
    format: String,
    version: u32,
    pub identity_cutoff: f32,
    pub analysis: Analysis,
    #[serde(default)]
    pub group_definitions: Vec<GeneGroupDefinition>,
}

impl AnalysisSession {
    pub fn new(
        analysis: Analysis,
        identity_cutoff: f32,
        group_definitions: Vec<GeneGroupDefinition>,
    ) -> Result<Self, SessionError> {
        let session = Self {
            format: ANALYSIS_SESSION_FORMAT.to_owned(),
            version: ANALYSIS_SESSION_VERSION,
            identity_cutoff,
            analysis,
            group_definitions,
        };
        session.validate()?;
        Ok(session)
    }

    pub fn from_json(bytes: &[u8]) -> Result<Self, SessionError> {
        let session: Self = serde_json::from_slice(bytes)?;
        session.validate()?;
        Ok(session)
    }

    pub fn to_json(&self) -> Result<Vec<u8>, SessionError> {
        self.validate()?;
        Ok(serde_json::to_vec_pretty(self)?)
    }

    pub fn format(&self) -> &str {
        &self.format
    }

    pub fn version(&self) -> u32 {
        self.version
    }

    pub fn validate(&self) -> Result<(), SessionError> {
        if self.format != ANALYSIS_SESSION_FORMAT {
            return Err(SessionError::UnsupportedFormat(self.format.clone()));
        }
        if self.version != ANALYSIS_SESSION_VERSION {
            return Err(SessionError::UnsupportedVersion(self.version));
        }
        if !self.identity_cutoff.is_finite() || !(0.0..=1.0).contains(&self.identity_cutoff) {
            return Err(SessionError::InvalidIdentity(self.identity_cutoff));
        }
        for (index, link) in self.analysis.links.iter().enumerate() {
            validate_gene_ref(&self.analysis, index, "query", link.query)?;
            validate_gene_ref(&self.analysis, index, "target", link.target)?;
            if link.query.cluster == link.target.cluster {
                return Err(SessionError::WithinClusterLink(index));
            }
            if !valid_fraction(link.identity) || !valid_fraction(link.similarity) {
                return Err(SessionError::InvalidLinkScore(index));
            }
        }
        Ok(())
    }
}

fn valid_fraction(value: f32) -> bool {
    value.is_finite() && (0.0..=1.0).contains(&value)
}

fn validate_gene_ref(
    analysis: &Analysis,
    link: usize,
    endpoint: &'static str,
    reference: GeneRef,
) -> Result<(), SessionError> {
    let valid = analysis
        .clusters
        .get(reference.cluster)
        .and_then(|cluster| cluster.loci.get(reference.locus))
        .and_then(|locus| locus.genes.get(reference.gene))
        .is_some();
    if valid {
        Ok(())
    } else {
        Err(SessionError::InvalidGeneReference {
            link,
            endpoint,
            cluster: reference.cluster,
            locus: reference.locus,
            gene: reference.gene,
        })
    }
}

#[derive(Debug, Error)]
pub enum SessionError {
    #[error("invalid analysis-session JSON: {0}")]
    Json(#[from] serde_json::Error),
    #[error("unsupported session format {0:?}")]
    UnsupportedFormat(String),
    #[error("unsupported analysis-session version {0}")]
    UnsupportedVersion(u32),
    #[error("session identity cutoff must be between 0 and 1, got {0}")]
    InvalidIdentity(f32),
    #[error("session link {link} has an invalid {endpoint} reference ({cluster}, {locus}, {gene})")]
    InvalidGeneReference {
        link: usize,
        endpoint: &'static str,
        cluster: usize,
        locus: usize,
        gene: usize,
    },
    #[error("session link {0} connects genes in the same cluster")]
    WithinClusterLink(usize),
    #[error("session link {0} has an identity or similarity outside 0..=1")]
    InvalidLinkScore(usize),
}

#[cfg(test)]
mod tests {
    use crate::{
        AnalysisOptions, GeneGroupDefinition, InputFile, analyse_genbank, parse_input_files,
    };

    use super::{ANALYSIS_SESSION_FORMAT, AnalysisSession, SessionError};

    const FIRST: &[u8] =
        br#"LOCUS       FIRST                      9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="first"
                     /translation="MKT"
ORIGIN
        1 atgaaaact
//
"#;
    const SECOND: &[u8] =
        br#"LOCUS       SECOND                     9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="second"
                     /translation="MKT"
ORIGIN
        1 atgaaaact
//
"#;
    const THIRD: &[u8] =
        br#"LOCUS       THIRD                      9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="third"
                     /translation="MKT"
ORIGIN
        1 atgaaaact
//
"#;

    #[test]
    fn round_trips_full_analysis_state() {
        let analysis = analyse_genbank(
            &[
                InputFile {
                    name: "first.gbk",
                    bytes: FIRST,
                },
                InputFile {
                    name: "second.gbk",
                    bytes: SECOND,
                },
            ],
            AnalysisOptions::default(),
        )
        .unwrap();
        let expected = AnalysisSession::new(
            analysis,
            0.3,
            vec![GeneGroupDefinition {
                label: "Example".into(),
                genes: vec!["first".into()],
                colour: Some("#123456".into()),
            }],
        )
        .unwrap();

        let encoded = expected.to_json().unwrap();
        let actual = AnalysisSession::from_json(&encoded).unwrap();

        assert_eq!(actual, expected);
        assert_eq!(actual.format(), ANALYSIS_SESSION_FORMAT);
        assert_eq!(
            actual.analysis.clusters[0].loci[0].genes[0].translation,
            "MKT"
        );
    }

    #[test]
    fn appending_clusters_preserves_old_links_and_adds_only_new_pairs() {
        let mut analysis = analyse_genbank(
            &[
                InputFile {
                    name: "first.gbk",
                    bytes: FIRST,
                },
                InputFile {
                    name: "second.gbk",
                    bytes: SECOND,
                },
            ],
            AnalysisOptions::default(),
        )
        .unwrap();
        let original_links = analysis.links.clone();
        let new_clusters = parse_input_files(&[InputFile {
            name: "third.gbk",
            bytes: THIRD,
        }])
        .unwrap();

        analysis.append_clusters(new_clusters, 0.3);

        assert_eq!(&analysis.links[..original_links.len()], original_links);
        assert_eq!(analysis.links.len(), 3);
        assert_eq!(analysis.links[1].target.cluster, 2);
        assert_eq!(analysis.links[2].target.cluster, 2);
    }

    #[test]
    fn rejects_invalid_versions_and_gene_references() {
        let mut value = serde_json::json!({
            "format": ANALYSIS_SESSION_FORMAT,
            "version": 2,
            "identityCutoff": 0.3,
            "analysis": { "clusters": [], "links": [] },
            "groupDefinitions": []
        });
        assert!(matches!(
            AnalysisSession::from_json(&serde_json::to_vec(&value).unwrap()),
            Err(SessionError::UnsupportedVersion(2))
        ));

        value["version"] = serde_json::json!(1);
        value["analysis"]["links"] = serde_json::json!([{
            "query": { "cluster": 0, "locus": 0, "gene": 0 },
            "target": { "cluster": 1, "locus": 0, "gene": 0 },
            "identity": 1.0,
            "similarity": 1.0
        }]);
        assert!(matches!(
            AnalysisSession::from_json(&serde_json::to_vec(&value).unwrap()),
            Err(SessionError::InvalidGeneReference { .. })
        ));
    }
}
