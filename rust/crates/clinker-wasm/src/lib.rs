//! Browser-facing adapter for `clinker-core`.
//!
//! This crate owns JavaScript value conversion only. Parsing and analysis stay
//! in the platform-independent core crate so the CLI and browser agree.

use clinker_core::{AnalysisError, AnalysisOptions, InputFile, PlotData, analyse_genbank};
use serde::Deserialize;
use wasm_bindgen::prelude::*;

#[derive(Debug, Deserialize)]
struct BrowserFile {
    name: String,
    bytes: Vec<u8>,
}

/// Analyse browser-uploaded GenBank files and return clustermap.js plot data.
///
/// `files` is an array of `{ name, bytes }` objects. In JavaScript, `bytes`
/// should be a `Uint8Array`; it is copied into Rust-owned memory before
/// analysis so the caller may release the original file buffers afterwards.
#[wasm_bindgen]
pub fn analyse(files: JsValue, identity_cutoff: f32) -> Result<JsValue, JsValue> {
    let files = serde_wasm_bindgen::from_value::<Vec<BrowserFile>>(files)
        .map_err(|error| JsValue::from_str(&format!("invalid browser input: {error}")))?;
    let plot_data = analyse_files(&files, identity_cutoff)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;

    serde_wasm_bindgen::to_value(&plot_data)
        .map_err(|error| JsValue::from_str(&format!("could not encode plot data: {error}")))
}

fn analyse_files(files: &[BrowserFile], identity_cutoff: f32) -> Result<PlotData, AnalysisError> {
    let inputs = files
        .iter()
        .map(|file| InputFile {
            name: &file.name,
            bytes: &file.bytes,
        })
        .collect::<Vec<_>>();
    analyse_genbank(&inputs, AnalysisOptions { identity_cutoff })
        .map(|analysis| analysis.to_plot_data())
}

#[cfg(test)]
mod tests {
    use super::{BrowserFile, analyse_files};

    #[test]
    fn adapter_delegates_to_the_shared_core() {
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

        let data = analyse_files(&files, 0.30).unwrap();
        assert_eq!(data.clusters.len(), 2);
        assert_eq!(data.links.len(), 1);
    }
}
