use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{
    Analysis, AnalysisOptions, AnalysisOptionsError, GeneGroupDefinition, KmerPrefilter,
    ValidationError,
};

pub const ANALYSIS_SESSION_FORMAT: &str = "clinker-analysis";
pub const ANALYSIS_SESSION_VERSION: u32 = 1;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(tag = "mode", rename_all = "camelCase", deny_unknown_fields)]
pub enum AlignmentSettings {
    Skipped,
    Aligned { prefilter: Option<KmerPrefilter> },
}

impl AlignmentSettings {
    pub fn options(self, identity_cutoff: f32) -> Option<AnalysisOptions> {
        match self {
            Self::Skipped => None,
            Self::Aligned { prefilter } => Some(AnalysisOptions {
                identity_cutoff,
                prefilter,
            }),
        }
    }

    pub fn validate(self, identity_cutoff: f32) -> Result<(), AnalysisOptionsError> {
        AnalysisOptions {
            identity_cutoff,
            prefilter: match self {
                Self::Skipped => None,
                Self::Aligned { prefilter } => prefilter,
            },
        }
        .validate()
    }
}

/// Versioned analysis state used to resume work or add more input clusters.
///
/// Unlike plot data, a session deliberately retains protein translations.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AnalysisSession {
    format: String,
    version: u32,
    identity_cutoff: f32,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    alignment: Option<AlignmentSettings>,
    pub analysis: Analysis,
    #[serde(default)]
    pub group_definitions: Vec<GeneGroupDefinition>,
}

impl AnalysisSession {
    pub fn new(
        analysis: Analysis,
        identity_cutoff: f32,
        alignment: AlignmentSettings,
        group_definitions: Vec<GeneGroupDefinition>,
    ) -> Result<Self, SessionError> {
        let session = Self {
            format: ANALYSIS_SESSION_FORMAT.to_owned(),
            version: ANALYSIS_SESSION_VERSION,
            identity_cutoff,
            alignment: Some(alignment),
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

    pub fn identity_cutoff(&self) -> f32 {
        self.identity_cutoff
    }

    pub fn alignment(&self) -> Option<AlignmentSettings> {
        self.alignment
    }

    pub fn ensure_settings(
        &self,
        identity_cutoff: f32,
        alignment: AlignmentSettings,
    ) -> Result<(), SessionError> {
        alignment.validate(identity_cutoff)?;
        let saved = self.alignment.ok_or(SessionError::UnknownSettings)?;
        if saved != alignment || identity_cutoff != self.identity_cutoff {
            return Err(SessionError::SettingsChanged);
        }
        Ok(())
    }

    pub fn validate(&self) -> Result<(), SessionError> {
        if self.format != ANALYSIS_SESSION_FORMAT {
            return Err(SessionError::UnsupportedFormat(self.format.clone()));
        }
        if self.version != ANALYSIS_SESSION_VERSION {
            return Err(SessionError::UnsupportedVersion(self.version));
        }
        self.alignment
            .unwrap_or(AlignmentSettings::Skipped)
            .validate(self.identity_cutoff)?;
        match self.alignment {
            Some(AlignmentSettings::Skipped) if !self.analysis.links.is_empty() => {
                return Err(SessionError::InconsistentAlignment);
            }
            Some(AlignmentSettings::Aligned { .. })
                if self
                    .analysis
                    .links
                    .iter()
                    .any(|link| link.identity < self.identity_cutoff) =>
            {
                return Err(SessionError::InconsistentAlignment);
            }
            _ => {}
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
    #[error(transparent)]
    Options(#[from] AnalysisOptionsError),
    #[error(
        "session alignment settings are unknown; use --reanalyse before extending it or changing settings"
    )]
    UnknownSettings,
    #[error(
        "requested analysis settings differ from the session; use --reanalyse to rebuild all comparisons"
    )]
    SettingsChanged,
    #[error("session links do not agree with its saved alignment settings")]
    InconsistentAlignment,
    #[error(transparent)]
    Analysis(#[from] ValidationError),
}

#[cfg(test)]
mod tests {
    use crate::{AnalysisOptions, GeneGroupDefinition, InputFile, analyse, parse_input_files};

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
        let analysis = analyse(
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
            crate::AlignmentSettings::Aligned { prefilter: None },
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
        let mut inconsistent = actual.clone();
        inconsistent.alignment = Some(crate::AlignmentSettings::Skipped);
        assert!(matches!(
            inconsistent.validate(),
            Err(SessionError::InconsistentAlignment)
        ));
        inconsistent.alignment = Some(crate::AlignmentSettings::Aligned {
            prefilter: Some(crate::KmerPrefilter {
                kmer_size: 0,
                min_shared_kmers: 3,
            }),
        });
        assert!(matches!(
            inconsistent.validate(),
            Err(SessionError::Options(
                crate::AnalysisOptionsError::Prefilter
            ))
        ));

        assert_eq!(
            actual.analysis.clusters[0].loci[0].genes[0].translation,
            "MKT"
        );
    }

    #[test]
    fn appending_clusters_preserves_old_links_and_adds_only_new_pairs() {
        let mut analysis = analyse(
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

        let before = analysis.clone();
        let mut invalid = new_clusters.clone();
        invalid[0].loci[0].genes[0].strand = 0;
        assert!(
            analysis
                .append_clusters(invalid.clone(), AnalysisOptions::default())
                .is_err()
        );
        assert!(
            analysis
                .append_clusters_without_alignment(invalid.clone())
                .is_err()
        );
        #[cfg(not(target_arch = "wasm32"))]
        assert!(
            analysis
                .append_clusters_with_jobs(invalid, AnalysisOptions::default(), 1)
                .is_err()
        );
        assert_eq!(analysis, before);
        analysis
            .append_clusters(new_clusters, AnalysisOptions::default())
            .unwrap();

        assert_eq!(&analysis.links[..original_links.len()], original_links);
        assert_eq!(analysis.links.len(), 3);
        assert_eq!(analysis.links[1].target.cluster, 2);
        assert_eq!(analysis.links[2].target.cluster, 2);
    }

    #[test]
    fn invalid_analysis_is_rejected_on_load_and_before_plotting() {
        let analysis = analyse(
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
        let original = AnalysisSession::new(
            analysis,
            0.3,
            crate::AlignmentSettings::Aligned { prefilter: None },
            Vec::new(),
        )
        .unwrap();
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
                invalid.analysis.cluster_similarity_matrix().is_err(),
                "case {case}"
            );
            assert!(invalid.analysis.cluster_order().is_err(), "case {case}");
            assert!(
                invalid
                    .analysis
                    .prepare_figure(crate::RowOrder::Similarity, &[], |_| {})
                    .is_err(),
                "case {case}"
            );
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
