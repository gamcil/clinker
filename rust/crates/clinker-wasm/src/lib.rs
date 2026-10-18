//! Browser-facing adapter for `clinker-core`.
//!
//! This crate owns JavaScript value conversion only. Parsing and analysis stay
//! in the platform-independent core crate so the CLI and browser agree.

use std::collections::BTreeMap;

use clinker_core::{
    Analysis, AnalysisSession, Cluster, ClusterPairSimilarity, GenBankOptions, Gene,
    GeneGroupDefinition, GeneRef, InputFile, KmerPrefilter, Link, Locus, PlotData, ProteinPair,
    SequenceRange, analyse_protein_pairs_with_progress, analyse_protein_tile_with_progress,
    can_reach_identity, cluster_kmer_candidate_pairs_with_progress, parse_input_files_with_options,
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
    #[serde(default)]
    genbank: BrowserGenBankOptions,
}

#[derive(Debug, Deserialize)]
#[serde(default, rename_all = "camelCase")]
struct BrowserGenBankOptions {
    ranges: Vec<BrowserSequenceRange>,
    as_separate_clusters: bool,
    #[serde(default = "default_set_origin")]
    set_origin: bool,
}

#[derive(Debug, Deserialize)]
struct BrowserSequenceRange {
    record: String,
    start: usize,
    end: usize,
}

fn default_set_origin() -> bool {
    true
}

impl Default for BrowserGenBankOptions {
    fn default() -> Self {
        Self {
            ranges: Vec::new(),
            as_separate_clusters: false,
            set_origin: true,
        }
    }
}

impl From<&BrowserGenBankOptions> for GenBankOptions {
    fn from(options: &BrowserGenBankOptions) -> Self {
        Self {
            ranges: options
                .ranges
                .iter()
                .map(|range| {
                    (
                        range.record.clone(),
                        SequenceRange {
                            start: range.start,
                            end: range.end,
                        },
                    )
                })
                .collect::<BTreeMap<_, _>>(),
            as_separate_clusters: options.as_separate_clusters,
            set_origin: options.set_origin,
        }
    }
}

/// Incremental parser used by the browser worker so JavaScript regains control
/// between files and can deliver and paint progress updates.
#[wasm_bindgen]
pub struct BrowserFileParser {
    clusters: Vec<Cluster>,
    file_positions: BTreeMap<String, usize>,
    pending_files: Vec<BrowserFile>,
    prefilter: BrowserPrefilter,
    genbank_options: GenBankOptions,
}

#[wasm_bindgen]
impl BrowserFileParser {
    #[wasm_bindgen(constructor)]
    pub fn new(prefilter: JsValue) -> Result<BrowserFileParser, JsValue> {
        let prefilter = serde_wasm_bindgen::from_value::<BrowserPrefilter>(prefilter)
            .map_err(|error| JsValue::from_str(&format!("invalid prefilter settings: {error}")))?;
        let genbank_options = GenBankOptions::from(&prefilter.genbank);
        Ok(Self {
            clusters: Vec::new(),
            file_positions: BTreeMap::new(),
            pending_files: Vec::new(),
            prefilter,
            genbank_options,
        })
    }

    pub fn parse_file(&mut self, file: JsValue) -> Result<(), JsValue> {
        let file = serde_wasm_bindgen::from_value::<BrowserFile>(file)
            .map_err(|error| JsValue::from_str(&format!("invalid browser input: {error}")))?;
        let position = self.file_positions.len();
        self.file_positions
            .entry(file.name.clone())
            .or_insert(position);
        if !is_genbank_name(&file.name) {
            self.pending_files.push(file);
            return Ok(());
        }
        let input = InputFile {
            name: &file.name,
            bytes: &file.bytes,
        };
        let mut clusters = parse_input_files_with_options(&[input], &self.genbank_options)
            .map_err(|error| JsValue::from_str(&error.to_string()))?;
        self.clusters.append(&mut clusters);
        Ok(())
    }

    pub fn finish(&mut self, progress: &Function) -> Result<JsValue, JsValue> {
        if !self.pending_files.is_empty() {
            let pending = std::mem::take(&mut self.pending_files);
            let mut clusters =
                parse_input_files_with_options(&input_files(&pending), &self.genbank_options)
                    .map_err(|error| JsValue::from_str(&error.to_string()))?;
            self.clusters.append(&mut clusters);
        }
        self.clusters.sort_by_key(|cluster| {
            cluster
                .source_file
                .as_ref()
                .and_then(|name| self.file_positions.get(name))
                .copied()
        });
        let clusters = std::mem::take(&mut self.clusters);
        let mut last_report = 0;
        let parsed = prepare_browser_files_with_filter_progress(
            clusters,
            &self.prefilter,
            |completed, total| {
                report_progress_event_if_due(progress, "filter", completed, total, &mut last_report)
            },
        );
        serde_wasm_bindgen::to_value(&parsed)
            .map_err(|error| JsValue::from_str(&format!("could not encode parsed files: {error}")))
    }
}

fn is_genbank_name(name: &str) -> bool {
    matches!(
        name.rsplit_once('.')
            .map(|(_, extension)| extension.to_ascii_lowercase())
            .as_deref(),
        Some("gb" | "gbk" | "genbank" | "gbf" | "gbff")
    )
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
    source_file: Option<String>,
    organism: Option<String>,
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

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct BrowserTileResult {
    links: Vec<TileLink>,
    alignment_count: usize,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct BrowserProgress<'a> {
    stage: &'a str,
    completed: usize,
    total: usize,
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

#[derive(Debug, Clone, Deserialize, Serialize)]
struct BrowserGroupDefinition {
    label: String,
    genes: Vec<String>,
    colour: Option<String>,
}

impl From<BrowserGroupDefinition> for GeneGroupDefinition {
    fn from(group: BrowserGroupDefinition) -> Self {
        Self {
            label: group.label,
            genes: group.genes,
            colour: group.colour,
        }
    }
}

impl From<&BrowserGroupDefinition> for GeneGroupDefinition {
    fn from(group: &BrowserGroupDefinition) -> Self {
        Self {
            label: group.label.clone(),
            genes: group.genes.clone(),
            colour: group.colour.clone(),
        }
    }
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
    cluster_names: Vec<String>,
    cluster_order: Vec<usize>,
    clustering_merges: Vec<[usize; 2]>,
    clustering_heights: Vec<f64>,
    similarity_matrix: Vec<Vec<ClusterPairSimilarity>>,
    analysis_session: AnalysisSession,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct BrowserSavedAnalysis {
    similarity_matrix: Vec<Vec<ClusterPairSimilarity>>,
    analysis_session: AnalysisSession,
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
    let alignment_count = query_sequences
        .iter()
        .flat_map(|query| {
            target_sequences.iter().filter(move |target| {
                can_reach_identity(query.len(), target.len(), identity_cutoff)
            })
        })
        .count();
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

    serde_wasm_bindgen::to_value(&BrowserTileResult {
        links,
        alignment_count,
    })
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
    let alignment_count = pairs
        .iter()
        .filter(|pair| {
            let (Some(query), Some(target)) = (
                sequences.get(pair.query_index),
                sequences.get(pair.target_index),
            ) else {
                return false;
            };
            can_reach_identity(query.len(), target.len(), identity_cutoff)
        })
        .count();
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

    serde_wasm_bindgen::to_value(&BrowserTileResult {
        links,
        alignment_count,
    })
    .map_err(|error| JsValue::from_str(&format!("could not encode tile links: {error}")))
}

fn report_progress(progress: &Function, processed: usize) {
    // `postMessage` only enqueues a parent-worker event; it does not yield or
    // interrupt the synchronous Wasm alignment loop.
    let _ = progress.call1(&JsValue::NULL, &JsValue::from_f64(processed as f64));
}

fn report_progress_event(progress: &Function, stage: &str, completed: usize, total: usize) {
    if let Ok(value) = serde_wasm_bindgen::to_value(&BrowserProgress {
        stage,
        completed,
        total,
    }) {
        let _ = progress.call1(&JsValue::NULL, &value);
    }
}

fn report_progress_event_if_due(
    progress: &Function,
    stage: &str,
    completed: usize,
    total: usize,
    last_report: &mut usize,
) {
    let interval = (total / 100).max(1);
    if completed == 0 || completed == total || completed.saturating_sub(*last_report) >= interval {
        report_progress_event(progress, stage, completed, total);
        *last_report = completed;
    }
}

/// Build groups and the default synteny ordering after browser tile work.
#[wasm_bindgen]
pub fn post_process(
    layout: JsValue,
    proteins: JsValue,
    links: JsValue,
    group_definitions: JsValue,
    identity_cutoff: f32,
    use_file_order: bool,
    progress: &Function,
) -> Result<JsValue, JsValue> {
    let layout = serde_wasm_bindgen::from_value::<BrowserLayout>(layout)
        .map_err(|error| JsValue::from_str(&format!("invalid layout data: {error}")))?;
    let proteins = serde_wasm_bindgen::from_value::<Vec<BrowserProtein>>(proteins)
        .map_err(|error| JsValue::from_str(&format!("invalid browser proteins: {error}")))?;
    let links = browser_links(links)?;
    let browser_group_definitions =
        serde_wasm_bindgen::from_value::<Vec<BrowserGroupDefinition>>(group_definitions)
            .map_err(|error| JsValue::from_str(&format!("invalid gene groups: {error}")))?;
    let group_definitions = browser_group_definitions
        .iter()
        .map(GeneGroupDefinition::from)
        .collect::<Vec<_>>();
    let analysis = Analysis {
        clusters: clusters_from_layout(layout, proteins)?,
        links,
    };
    post_process_analysis(
        analysis,
        group_definitions,
        identity_cutoff,
        use_file_order,
        progress,
    )
}

/// Validate saved analysis and recover export scores without rebuilding its figure.
#[wasm_bindgen]
pub fn open_session(session: JsValue) -> Result<JsValue, JsValue> {
    let session = serde_wasm_bindgen::from_value::<AnalysisSession>(session)
        .map_err(|error| JsValue::from_str(&format!("invalid analysis session: {error}")))?;
    session
        .validate()
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let result = BrowserSavedAnalysis {
        similarity_matrix: session.analysis.cluster_similarity_matrix(),
        analysis_session: session,
    };
    serde_wasm_bindgen::to_value(&result)
        .map_err(|error| JsValue::from_str(&format!("could not encode saved analysis: {error}")))
}

fn post_process_analysis(
    analysis: Analysis,
    group_definitions: Vec<GeneGroupDefinition>,
    identity_cutoff: f32,
    use_file_order: bool,
    progress: &Function,
) -> Result<JsValue, JsValue> {
    let cluster_count = analysis.clusters.len();
    let similarity_total = cluster_count.saturating_mul(cluster_count.saturating_sub(1)) / 2;
    let mut last_similarity_report = 0;
    let similarity_matrix = analysis.cluster_similarity_matrix_with_progress(|completed| {
        report_progress_event_if_due(
            progress,
            "similarity",
            completed,
            similarity_total,
            &mut last_similarity_report,
        )
    });
    report_progress_event(progress, "ordering", 0, cluster_count);
    let ordering = if use_file_order {
        None
    } else {
        Some(analysis.cluster_ordering_from_similarity_matrix(&similarity_matrix))
    };
    let (order, clustering_merges, clustering_heights) = ordering.map_or_else(
        || ((0..cluster_count).collect(), Vec::new(), Vec::new()),
        |ordering| (ordering.order, ordering.merges, ordering.heights),
    );
    report_progress_event(progress, "ordering", cluster_count, cluster_count);
    let group_total = analysis.links.len() * 2;
    let mut last_group_report = 0;
    let mut last_layout_report = 0;
    let arranged = analysis.to_auto_arranged_plot_data_with_group_definitions_and_progress(
        &order,
        &group_definitions,
        |completed| {
            report_progress_event_if_due(
                progress,
                "groups",
                completed,
                group_total,
                &mut last_group_report,
            )
        },
        |completed| {
            report_progress_event_if_due(
                progress,
                "layout",
                completed,
                cluster_count,
                &mut last_layout_report,
            )
        },
    );
    let cluster_names = analysis
        .clusters
        .iter()
        .map(|cluster| cluster.name.clone())
        .collect();
    let analysis_session = AnalysisSession::new(analysis, identity_cutoff, group_definitions)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let result = BrowserPostProcess {
        plot_data: arranged,
        cluster_names,
        cluster_order: order,
        clustering_merges,
        clustering_heights,
        similarity_matrix,
        analysis_session,
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

fn clusters_from_layout(
    layout: BrowserLayout,
    proteins: Vec<BrowserProtein>,
) -> Result<Vec<Cluster>, JsValue> {
    let mut clusters = layout
        .clusters
        .into_iter()
        .map(|cluster| Cluster {
            name: cluster.name,
            source_file: cluster.source_file,
            organism: cluster.organism,
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
        .collect::<Vec<_>>();
    let expected = clusters
        .iter()
        .flat_map(|cluster| &cluster.loci)
        .map(|locus| locus.genes.len())
        .sum::<usize>();
    let mut assigned = std::collections::BTreeSet::new();
    for protein in proteins {
        let reference = (protein.cluster, protein.locus, protein.gene);
        if !assigned.insert(reference) {
            return Err(JsValue::from_str(&format!(
                "duplicate browser protein reference {reference:?}"
            )));
        }
        let gene = clusters
            .get_mut(protein.cluster)
            .and_then(|cluster| cluster.loci.get_mut(protein.locus))
            .and_then(|locus| locus.genes.get_mut(protein.gene))
            .ok_or_else(|| {
                JsValue::from_str(&format!("invalid browser protein reference {reference:?}"))
            })?;
        gene.translation = protein.translation;
    }
    if assigned.len() != expected {
        return Err(JsValue::from_str(
            "browser analysis is missing one or more protein sequences",
        ));
    }
    Ok(clusters)
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

fn prepare_browser_files_with_filter_progress<F>(
    clusters: Vec<Cluster>,
    prefilter: &BrowserPrefilter,
    mut progress: F,
) -> ParsedBrowserFiles
where
    F: FnMut(usize, usize),
{
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
        let sequences = proteins
            .iter()
            .map(|protein| protein.translation.as_bytes())
            .collect::<Vec<_>>();
        let clusters = proteins
            .iter()
            .map(|protein| protein.cluster)
            .collect::<Vec<_>>();
        cluster_kmer_candidate_pairs_with_progress(
            &sequences,
            &clusters,
            prefilter.identity_cutoff,
            KmerPrefilter {
                kmer_size: prefilter.kmer_size,
                min_shared_kmers: prefilter.min_shared_kmers,
            },
            &mut progress,
        )
        .into_iter()
        .map(|pair| BrowserProteinPair {
            query_index: pair.query_index,
            target_index: pair.target_index,
        })
        .collect()
    });
    ParsedBrowserFiles {
        layout,
        proteins,
        candidate_pairs,
    }
}

fn layout_from_clusters(clusters: &[Cluster]) -> BrowserLayout {
    BrowserLayout {
        clusters: clusters
            .iter()
            .map(|cluster| BrowserLayoutCluster {
                name: cluster.name.clone(),
                source_file: cluster.source_file.clone(),
                organism: cluster.organism.clone(),
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
