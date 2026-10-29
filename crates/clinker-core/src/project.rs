//! Shared project envelope. Figure and editor state belong to clustermap and
//! remain JSON values so native round-trips preserve its extension fields.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use thiserror::Error;

use crate::{AnalysisSession, PlotData, SessionError};

pub const PROJECT_FORMAT: &str = "clinker-project";
pub const PROJECT_VERSION: u32 = 2;

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct ClinkerProject {
    format: String,
    version: u32,
    pub analysis: AnalysisSession,
    pub data: Value,
    pub config: Value,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub state: Option<Value>,
    #[serde(flatten)]
    pub extensions: BTreeMap<String, Value>,
}

impl ClinkerProject {
    pub fn new(analysis: AnalysisSession, data: &PlotData) -> Result<Self, ProjectError> {
        let project = Self {
            format: PROJECT_FORMAT.to_owned(),
            version: PROJECT_VERSION,
            analysis,
            data: serde_json::to_value(data)?,
            config: json!({}),
            state: None,
            extensions: BTreeMap::new(),
        };
        project.validate()?;
        Ok(project)
    }

    pub fn from_json(bytes: &[u8]) -> Result<Self, ProjectError> {
        let project: Self = serde_json::from_slice(bytes)?;
        project.validate()?;
        Ok(project)
    }

    pub fn to_json(&self) -> Result<Vec<u8>, ProjectError> {
        self.validate()?;
        Ok(serde_json::to_vec_pretty(self)?)
    }

    /// Check the envelope and analysis. Renderer-specific IDs, trees,
    /// annotations and layout are validated by the shared JavaScript viewer.
    pub fn validate(&self) -> Result<(), ProjectError> {
        if self.format != PROJECT_FORMAT || self.version != PROJECT_VERSION {
            return Err(ProjectError::UnsupportedProject);
        }
        self.analysis.validate()?;
        if !self.config.is_object()
            || !self.data.is_object()
            || ["clusters", "links", "groups"]
                .iter()
                .any(|key| !self.data[*key].is_array())
        {
            return Err(ProjectError::InvalidFigure);
        }
        if let Some(state) = &self.state
            && (!state.is_object()
                || !matches!(state["version"].as_u64(), Some(1 | 2))
                || !state["clusterOrder"].is_array())
        {
            return Err(ProjectError::InvalidFigure);
        }
        if self.extensions.keys().any(|key| {
            matches!(
                key.as_str(),
                "format" | "version" | "analysis" | "data" | "config" | "state"
            )
        }) {
            return Err(ProjectError::InvalidFigure);
        }
        Ok(())
    }
}

#[derive(Debug, Error)]
pub enum ProjectError {
    #[error("not a supported clinker project (expected version 2)")]
    UnsupportedProject,
    #[error("project contains invalid figure settings, data or layout state")]
    InvalidFigure,
    #[error(transparent)]
    Analysis(#[from] SessionError),
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}
