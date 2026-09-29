//! Browser-facing adapter for `clinker-core`.
//!
//! This crate owns JavaScript value conversion only. Parsing and analysis stay
//! in the platform-independent core crate so the CLI and browser agree.

use clinker_core::{
    Analysis, AnalysisError, InputFile, PlotData, analyse_protein_tile, parse_input_files,
};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[derive(Debug, Deserialize)]
struct BrowserFile {
    name: String,
    bytes: Vec<u8>,
}

/// The compact alignment input for one parsed CDS gene.
#[derive(Debug, Clone, Deserialize, Serialize)]
struct BrowserProtein {
    cluster: usize,
    locus: usize,
    gene: usize,
    translation: String,
}

/// Metadata and alignment inputs produced by the one-time parsing stage.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct ParsedBrowserFiles {
    plot_data: PlotData,
    proteins: Vec<BrowserProtein>,
}

/// A retained link returned by an alignment-only worker tile.
#[derive(Debug, Serialize)]
struct TileLink {
    query: BrowserProteinRef,
    target: BrowserProteinRef,
    identity: f32,
    similarity: f32,
}

#[derive(Debug, Serialize)]
struct BrowserProteinRef {
    cluster: usize,
    locus: usize,
    gene: usize,
}

impl From<&BrowserProtein> for BrowserProteinRef {
    fn from(protein: &BrowserProtein) -> Self {
        Self {
            cluster: protein.cluster,
            locus: protein.locus,
            gene: protein.gene,
        }
    }
}

/// Parse browser-uploaded GenBank files once and return plot metadata plus
/// compact protein records for alignment workers.
///
/// `files` is an array of `{ name, bytes }` objects. In JavaScript, `bytes`
/// should be a `Uint8Array`; it is copied into Rust-owned memory before
/// analysis so the caller may release the original file buffers afterwards.
#[wasm_bindgen]
pub fn parse_files(files: JsValue) -> Result<JsValue, JsValue> {
    let files = serde_wasm_bindgen::from_value::<Vec<BrowserFile>>(files)
        .map_err(|error| JsValue::from_str(&format!("invalid browser input: {error}")))?;
    let parsed =
        parse_browser_files(&files).map_err(|error| JsValue::from_str(&error.to_string()))?;
    serde_wasm_bindgen::to_value(&parsed)
        .map_err(|error| JsValue::from_str(&format!("could not encode parsed files: {error}")))
}

/// Align the Cartesian product of two compact protein blocks.
#[wasm_bindgen]
pub fn analyse_tile(
    query: JsValue,
    target: JsValue,
    identity_cutoff: f32,
) -> Result<JsValue, JsValue> {
    let query = serde_wasm_bindgen::from_value::<Vec<BrowserProtein>>(query)
        .map_err(|error| JsValue::from_str(&format!("invalid query proteins: {error}")))?;
    let target = serde_wasm_bindgen::from_value::<Vec<BrowserProtein>>(target)
        .map_err(|error| JsValue::from_str(&format!("invalid target proteins: {error}")))?;
    let query_sequences = query
        .iter()
        .map(|protein| protein.translation.as_bytes())
        .collect::<Vec<_>>();
    let target_sequences = target
        .iter()
        .map(|protein| protein.translation.as_bytes())
        .collect::<Vec<_>>();
    let links = analyse_protein_tile(&query_sequences, &target_sequences, identity_cutoff)
        .into_iter()
        .map(|alignment| TileLink {
            query: (&query[alignment.query_index]).into(),
            target: (&target[alignment.target_index]).into(),
            identity: alignment.identity,
            similarity: alignment.similarity,
        })
        .collect::<Vec<_>>();

    serde_wasm_bindgen::to_value(&links)
        .map_err(|error| JsValue::from_str(&format!("could not encode tile links: {error}")))
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

fn parse_browser_files(files: &[BrowserFile]) -> Result<ParsedBrowserFiles, AnalysisError> {
    let clusters = parse_input_files(&input_files(files))?;
    let proteins = clusters
        .iter()
        .enumerate()
        .flat_map(|(cluster_index, cluster)| {
            cluster
                .loci
                .iter()
                .enumerate()
                .flat_map(move |(locus_index, locus)| {
                    locus
                        .genes
                        .iter()
                        .enumerate()
                        .map(move |(gene_index, gene)| BrowserProtein {
                            cluster: cluster_index,
                            locus: locus_index,
                            gene: gene_index,
                            translation: gene.translation.clone(),
                        })
                })
        })
        .collect();
    let plot_data = Analysis {
        clusters,
        links: Vec::new(),
    }
    .to_plot_data();
    Ok(ParsedBrowserFiles {
        plot_data,
        proteins,
    })
}

#[cfg(test)]
mod tests {
    use super::{BrowserFile, analyse_tile, parse_browser_files};
    use wasm_bindgen::JsValue;

    #[test]
    fn adapter_separates_parsing_from_tile_alignment() {
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
        assert_eq!(data.plot_data.clusters.len(), 2);
        assert!(data.plot_data.links.is_empty());
        assert_eq!(data.proteins.len(), 2);

        // The exported function is exercised by browser/WASM integration; the
        // parsing test above keeps native unit tests independent of JsValue.
        let _ = analyse_tile as fn(JsValue, JsValue, f32) -> Result<JsValue, JsValue>;
        assert_eq!(data.proteins[0].translation, "MA*");
    }
}
