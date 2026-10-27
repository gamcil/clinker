use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{Analysis, GeneGroupDefinition, ValidationError};

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
        self.analysis.validate()?;
        Ok(())
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
    #[error(transparent)]
    Analysis(#[from] ValidationError),
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
    fn invalid_analysis_is_rejected_on_load_and_before_plotting() {
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
        let original = AnalysisSession::new(analysis, 0.3, Vec::new()).unwrap();
        let edits: &[fn(&mut crate::Analysis)] = &[
            |a| a.clusters[0].loci.clear(),
            |a| a.clusters[0].loci[0].end = 0,
            |a| a.clusters[0].loci[0].genes[0].end = 10,
            |a| a.clusters[0].loci[0].genes[0].strand = 0,
            |a| a.links[0].query.gene = usize::MAX,
            |a| {
                let link = &mut a.links[0];
                std::mem::swap(&mut link.query, &mut link.target);
            },
            |a| a.links[0].identity = f32::NAN,
            |a| a.links[0].similarity = 0.0,
            |a| a.links.push(a.links[0].clone()),
        ];
        for (case, edit) in edits.iter().enumerate() {
            let mut invalid = original.clone();
            edit(&mut invalid.analysis);
            assert!(invalid.validate().is_err(), "case {case}");
            assert!(
                invalid
                    .analysis
                    .to_auto_arranged_plot_data(&[0, 1])
                    .is_err(),
                "case {case}"
            );
            let bytes = serde_json::to_vec(&invalid).unwrap();
            assert!(AnalysisSession::from_json(&bytes).is_err(), "case {case}");
        }
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
            Err(SessionError::Analysis(
                crate::ValidationError::GeneReference { .. }
            ))
        ));
    }
}
