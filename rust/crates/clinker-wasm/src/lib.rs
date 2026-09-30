//! Browser-facing adapter for `clinker-core`.
//!
//! This crate owns JavaScript value conversion only. Parsing and analysis stay
//! in the platform-independent core crate so the CLI and browser agree.

use clinker_core::{
    Analysis, AnalysisError, Cluster, Gene, GeneRef, InputFile, KmerPrefilter, Link, Locus,
    PlotData, ProteinPair, analyse_protein_pairs_with_progress, analyse_protein_tile_with_progress,
    kmer_candidate_pairs, parse_input_files,
};
use js_sys::Function;
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

/// Indices into `ParsedBrowserFiles.proteins`, not sequence copies.
#[derive(Debug, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
struct BrowserProteinPair {
    query_index: usize,
    target_index: usize,
}

/// Metadata and alignment inputs produced by the one-time parsing stage.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct ParsedBrowserFiles {
    layout: BrowserLayout,
    proteins: Vec<BrowserProtein>,
    candidate_pairs: Option<Vec<BrowserProteinPair>>,
}

/// Optional approximate candidate selection for browser alignment work.
#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct BrowserPrefilter {
    enabled: bool,
    kmer_size: usize,
    min_shared_kmers: usize,
    identity_cutoff: f32,
}

/// Coordinate and annotation metadata retained between parsing and the final
/// layout pass. Protein sequences stay only in `BrowserProtein` tile inputs.
#[derive(Debug, Deserialize, Serialize)]
struct BrowserLayout {
    clusters: Vec<BrowserLayoutCluster>,
}

#[derive(Debug, Deserialize, Serialize)]
struct BrowserLayoutCluster {
    name: String,
    loci: Vec<BrowserLayoutLocus>,
}

#[derive(Debug, Deserialize, Serialize)]
struct BrowserLayoutLocus {
    name: String,
    start: usize,
    end: usize,
    genes: Vec<BrowserLayoutGene>,
}

#[derive(Debug, Deserialize, Serialize)]
struct BrowserLayoutGene {
    label: String,
    names: Vec<(String, String)>,
    start: usize,
    end: usize,
    strand: i8,
}

/// A retained link returned by an alignment-only worker tile.
#[derive(Debug, Serialize)]
struct TileLink {
    query: BrowserProteinRef,
    target: BrowserProteinRef,
    identity: f32,
    similarity: f32,
}

/// Link shape returned by an alignment worker and later sent to the grouping
/// worker. It uses positions, rather than renderer UIDs, so Rust can form
/// groups without reparsing strings.
#[derive(Debug, Deserialize)]
struct BrowserLink {
    query: BrowserProteinRef,
    target: BrowserProteinRef,
    identity: f32,
    similarity: f32,
}

#[derive(Debug, Deserialize, Serialize)]
struct BrowserProteinRef {
    cluster: usize,
    locus: usize,
    gene: usize,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct BrowserPostProcess {
    plot_data: PlotData,
}

impl From<BrowserProteinRef> for GeneRef {
    fn from(reference: BrowserProteinRef) -> Self {
        Self {
            cluster: reference.cluster,
            locus: reference.locus,
            gene: reference.gene,
        }
    }
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
pub fn parse_files(files: JsValue, prefilter: JsValue) -> Result<JsValue, JsValue> {
    let files = serde_wasm_bindgen::from_value::<Vec<BrowserFile>>(files)
        .map_err(|error| JsValue::from_str(&format!("invalid browser input: {error}")))?;
    let prefilter = serde_wasm_bindgen::from_value::<BrowserPrefilter>(prefilter)
        .map_err(|error| JsValue::from_str(&format!("invalid prefilter settings: {error}")))?;
    let parsed = parse_browser_files(&files, prefilter)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    serde_wasm_bindgen::to_value(&parsed)
        .map_err(|error| JsValue::from_str(&format!("could not encode parsed files: {error}")))
}

/// Align the Cartesian product of two compact protein blocks.
#[wasm_bindgen]
pub fn analyse_tile(
    query: JsValue,
    target: JsValue,
    identity_cutoff: f32,
    progress: &Function,
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
    let links = analyse_protein_tile_with_progress(
        &query_sequences,
        &target_sequences,
        identity_cutoff,
        |processed| report_progress(progress, processed),
    )
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

/// Align selected pairs from one compact protein tile.
#[wasm_bindgen]
pub fn analyse_pairs(
    proteins: JsValue,
    pairs: JsValue,
    identity_cutoff: f32,
    progress: &Function,
) -> Result<JsValue, JsValue> {
    let proteins = serde_wasm_bindgen::from_value::<Vec<BrowserProtein>>(proteins)
        .map_err(|error| JsValue::from_str(&format!("invalid tile proteins: {error}")))?;
    let pairs = serde_wasm_bindgen::from_value::<Vec<BrowserProteinPair>>(pairs)
        .map_err(|error| JsValue::from_str(&format!("invalid tile pairs: {error}")))?;
    let sequences = proteins
        .iter()
        .map(|protein| protein.translation.as_bytes())
        .collect::<Vec<_>>();
    let pairs = pairs
        .iter()
        .map(|pair| ProteinPair {
            query_index: pair.query_index,
            target_index: pair.target_index,
        })
        .collect::<Vec<_>>();
    let links =
        analyse_protein_pairs_with_progress(&sequences, &pairs, identity_cutoff, |processed| {
            report_progress(progress, processed)
        })
        .into_iter()
        .map(|alignment| TileLink {
            query: (&proteins[alignment.query_index]).into(),
            target: (&proteins[alignment.target_index]).into(),
            identity: alignment.identity,
            similarity: alignment.similarity,
        })
        .collect::<Vec<_>>();

    serde_wasm_bindgen::to_value(&links)
        .map_err(|error| JsValue::from_str(&format!("could not encode tile links: {error}")))
}

fn report_progress(progress: &Function, processed: usize) {
    // `postMessage` only enqueues a parent-worker event; it does not yield or
    // interrupt the synchronous Wasm alignment loop.
    let _ = progress.call1(&JsValue::NULL, &JsValue::from_f64(processed as f64));
}

/// Build groups and the default synteny ordering after browser tile work.
#[wasm_bindgen]
pub fn post_process(layout: JsValue, links: JsValue) -> Result<JsValue, JsValue> {
    let layout = serde_wasm_bindgen::from_value::<BrowserLayout>(layout)
        .map_err(|error| JsValue::from_str(&format!("invalid layout data: {error}")))?;
    let links = browser_links(links)?;
    let analysis = Analysis {
        clusters: clusters_from_layout(layout),
        links,
    };
    let order = analysis.cluster_order(clinker_core::DEFAULT_CONTIGUITY_WEIGHT);
    let arranged = analysis.to_auto_arranged_plot_data(&order);
    let result = BrowserPostProcess {
        plot_data: arranged,
    };
    serde_wasm_bindgen::to_value(&result)
        .map_err(|error| JsValue::from_str(&format!("could not encode post-processing: {error}")))
}

fn browser_links(links: JsValue) -> Result<Vec<Link>, JsValue> {
    serde_wasm_bindgen::from_value::<Vec<BrowserLink>>(links)
        .map_err(|error| JsValue::from_str(&format!("invalid browser links: {error}")))
        .map(|links| {
            links
                .into_iter()
                .map(|link| Link {
                    query: link.query.into(),
                    target: link.target.into(),
                    identity: link.identity,
                    similarity: link.similarity,
                })
                .collect()
        })
}

fn clusters_from_layout(layout: BrowserLayout) -> Vec<Cluster> {
    layout
        .clusters
        .into_iter()
        .map(|cluster| Cluster {
            name: cluster.name,
            loci: cluster
                .loci
                .into_iter()
                .map(|locus| Locus {
                    name: locus.name,
                    start: locus.start,
                    end: locus.end,
                    genes: locus
                        .genes
                        .into_iter()
                        .map(|gene| Gene {
                            label: gene.label,
                            names: gene.names.into_iter().collect(),
                            start: gene.start,
                            end: gene.end,
                            strand: gene.strand,
                            translation: String::new(),
                        })
                        .collect(),
                })
                .collect(),
        })
        .collect()
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

fn parse_browser_files(
    files: &[BrowserFile],
    prefilter: BrowserPrefilter,
) -> Result<ParsedBrowserFiles, AnalysisError> {
    let clusters = parse_input_files(&input_files(files))?;
    let layout = layout_from_clusters(&clusters);
    let proteins: Vec<BrowserProtein> = clusters
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
    let candidate_pairs = prefilter.enabled.then(|| {
        let by_cluster = proteins_by_cluster(&proteins, clusters.len());
        let settings = KmerPrefilter {
            kmer_size: prefilter.kmer_size,
            min_shared_kmers: prefilter.min_shared_kmers,
        };
        let mut pairs = Vec::new();
        for query_cluster in 0..by_cluster.len() {
            for target_cluster in query_cluster + 1..by_cluster.len() {
                let query = &by_cluster[query_cluster];
                let target = &by_cluster[target_cluster];
                let query_sequences = query
                    .iter()
                    .map(|&index| proteins[index].translation.as_bytes())
                    .collect::<Vec<_>>();
                let target_sequences = target
                    .iter()
                    .map(|&index| proteins[index].translation.as_bytes())
                    .collect::<Vec<_>>();
                pairs.extend(
                    kmer_candidate_pairs(
                        &query_sequences,
                        &target_sequences,
                        prefilter.identity_cutoff,
                        settings,
                    )
                    .into_iter()
                    .map(|pair| BrowserProteinPair {
                        query_index: query[pair.query_index],
                        target_index: target[pair.target_index],
                    }),
                );
            }
        }
        pairs
    });
    Ok(ParsedBrowserFiles {
        layout,
        proteins,
        candidate_pairs,
    })
}

fn proteins_by_cluster(proteins: &[BrowserProtein], cluster_count: usize) -> Vec<Vec<usize>> {
    let mut clusters = vec![Vec::new(); cluster_count];
    for (index, protein) in proteins.iter().enumerate() {
        clusters[protein.cluster].push(index);
    }
    clusters
}

fn layout_from_clusters(clusters: &[Cluster]) -> BrowserLayout {
    BrowserLayout {
        clusters: clusters
            .iter()
            .map(|cluster| BrowserLayoutCluster {
                name: cluster.name.clone(),
                loci: cluster
                    .loci
                    .iter()
                    .map(|locus| BrowserLayoutLocus {
                        name: locus.name.clone(),
                        start: locus.start,
                        end: locus.end,
                        genes: locus
                            .genes
                            .iter()
                            .map(|gene| BrowserLayoutGene {
                                label: gene.label.clone(),
                                names: gene.names.clone(),
                                start: gene.start,
                                end: gene.end,
                                strand: gene.strand,
                            })
                            .collect(),
                    })
                    .collect(),
            })
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BrowserFile, BrowserPrefilter, Function, analyse_pairs, analyse_tile, parse_browser_files,
        post_process,
    };
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

        let data = parse_browser_files(
            &files,
            BrowserPrefilter {
                enabled: false,
                kmer_size: 3,
                min_shared_kmers: 3,
                identity_cutoff: 0.30,
            },
        )
        .unwrap();
        assert_eq!(data.layout.clusters.len(), 2);
        assert_eq!(data.proteins.len(), 2);

        // The exported function is exercised by browser/WASM integration; the
        // parsing test above keeps native unit tests independent of JsValue.
        let _ = analyse_tile as fn(JsValue, JsValue, f32, &Function) -> Result<JsValue, JsValue>;
        let _ = analyse_pairs as fn(JsValue, JsValue, f32, &Function) -> Result<JsValue, JsValue>;
        let _ = post_process as fn(JsValue, JsValue) -> Result<JsValue, JsValue>;
        assert_eq!(data.proteins[0].translation, "MA*");
    }
}
