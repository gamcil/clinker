//! Browser-facing adapter for `clinker-core`.
//!
//! This crate owns JavaScript value conversion only. Parsing and analysis stay
//! in the platform-independent core crate so the CLI and browser agree.

use clinker_core::{
    Analysis, AnalysisError, GeneRef, InputFile, Link, PlotData, analyse_cluster_pair,
    parse_input_files,
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[derive(Debug, Deserialize)]
struct BrowserFile {
    name: String,
    bytes: Vec<u8>,
}

/// A link returned by one cluster-pair worker.
#[derive(Debug, Serialize)]
struct PairLink {
    query: PairGeneRef,
    target: PairGeneRef,
    identity: f32,
    similarity: f32,
}

#[derive(Debug, Serialize)]
struct PairGeneRef {
    locus: usize,
    gene: usize,
}

impl From<GeneRef> for PairGeneRef {
    fn from(reference: GeneRef) -> Self {
        Self {
            locus: reference.locus,
            gene: reference.gene,
        }
    }
}

impl From<Link> for PairLink {
    fn from(link: Link) -> Self {
        Self {
            query: link.query.into(),
            target: link.target.into(),
            identity: link.identity,
            similarity: link.similarity,
        }
    }
}

/// Parse browser-uploaded GenBank files and return plot-ready cluster metadata.
///
/// `files` is an array of `{ name, bytes }` objects. In JavaScript, `bytes`
/// should be a `Uint8Array`; it is copied into Rust-owned memory before
/// analysis so the caller may release the original file buffers afterwards.
#[wasm_bindgen]
pub fn parse_files(files: JsValue) -> Result<JsValue, JsValue> {
    let files = serde_wasm_bindgen::from_value::<Vec<BrowserFile>>(files)
        .map_err(|error| JsValue::from_str(&format!("invalid browser input: {error}")))?;
    let plot_data =
        parse_browser_files(&files).map_err(|error| JsValue::from_str(&error.to_string()))?;
    serde_wasm_bindgen::to_value(&plot_data)
        .map_err(|error| JsValue::from_str(&format!("could not encode plot data: {error}")))
}

/// Align exactly two browser-uploaded GenBank files and return their links.
#[wasm_bindgen]
pub fn analyse_pair(files: JsValue, identity_cutoff: f32) -> Result<JsValue, JsValue> {
    let files = serde_wasm_bindgen::from_value::<Vec<BrowserFile>>(files)
        .map_err(|error| JsValue::from_str(&format!("invalid browser input: {error}")))?;
    if files.len() != 2 {
        return Err(JsValue::from_str(
            "a pairwise analysis requires exactly two files",
        ));
    }
    let links = analyse_browser_pair(&files, identity_cutoff)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    serde_wasm_bindgen::to_value(&links)
        .map_err(|error| JsValue::from_str(&format!("could not encode pair links: {error}")))
}

fn input_files(files: &[BrowserFile]) -> Vec<InputFile<'_>> {
    files
        .iter()
        .map(|file| InputFile {
            name: &file.name,
            bytes: &file.bytes,
        })
        .collect()
}

fn parse_browser_files(files: &[BrowserFile]) -> Result<PlotData, AnalysisError> {
    let clusters = parse_input_files(&input_files(files))?;
    Ok(Analysis {
        clusters,
        links: Vec::new(),
    }
    .to_plot_data())
}

fn analyse_browser_pair(
    files: &[BrowserFile],
    identity_cutoff: f32,
) -> Result<Vec<PairLink>, AnalysisError> {
    let clusters = parse_input_files(&input_files(files))?;
    Ok(
        analyse_cluster_pair(&clusters[0], &clusters[1], identity_cutoff)
            .into_iter()
            .map(PairLink::from)
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::{BrowserFile, analyse_browser_pair, parse_browser_files};

    #[test]
    fn adapter_separates_parsing_from_pairwise_alignment() {
        let record = b"LOCUS       TEST                       9 bp    DNA     linear   UNA 01-JAN-2000\nFEATURES             Location/Qualifiers\n     CDS             1..9\n                     /locus_tag=\"test\"\nORIGIN\n        1 atggcttaa\n//\n";
        let files = vec![
            BrowserFile {
                name: "first.gbk".into(),
                bytes: record.to_vec(),
            },
            BrowserFile {
                name: "second.gbk".into(),
                bytes: record.to_vec(),
            },
        ];

        let data = parse_browser_files(&files).unwrap();
        let links = analyse_browser_pair(&files, 0.30).unwrap();
        assert_eq!(data.clusters.len(), 2);
        assert!(data.links.is_empty());
        assert_eq!(links.len(), 1);
        assert_eq!(links[0].query.locus, 0);
        assert_eq!(links[0].target.gene, 0);
    }
}
