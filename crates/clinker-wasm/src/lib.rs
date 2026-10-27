//! Browser-facing adapter for `clinker-core`.
//!
//! This crate owns JavaScript value conversion only. Parsing and analysis stay
//! in the platform-independent core crate so the CLI and browser agree.

use std::collections::BTreeMap;

use clinker_core::{
    AlignmentCounts, AlignmentSettings, Analysis, AnalysisOptions, AnalysisSession, Cluster,
    ClusterPairSimilarity, GenBankOptions, Gene, GeneGroupDefinition, GeneRef, InputFile,
    KmerPrefilter, Link, Locus, PlotData, PreparedAlignment, ProteinBatch, ProteinBatchAligner,
    ProteinBatchResult, ProteinPair, SequenceRange, parse_input_files_with_options,
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

/// Metadata and alignment inputs produced by the one-time parsing stage.
#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct ParsedBrowserFiles {
    layout: BrowserLayout,
    proteins: Vec<BrowserProtein>,
    candidate_pairs: Option<Vec<ProteinPair>>,
    alignment_batches: Vec<ProteinBatch>,
    alignment_counts: AlignmentCounts,
}

/// Analysis and parsing settings supplied by the browser host.
#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct BrowserAnalysisSettings {
    enabled: bool,
    #[serde(default)]
    no_align: bool,
    kmer_size: usize,
    min_shared_kmers: usize,
    identity_cutoff: f32,
    #[serde(default)]
    genbank: BrowserGenBankOptions,
}

impl BrowserAnalysisSettings {
    fn alignment(&self) -> AlignmentSettings {
        if self.no_align {
            AlignmentSettings::Skipped
        } else {
            AlignmentSettings::Aligned {
                prefilter: self.enabled.then_some(KmerPrefilter {
                    kmer_size: self.kmer_size,
                    min_shared_kmers: self.min_shared_kmers,
                }),
            }
        }
    }

    fn validate(&self) -> Result<(), JsValue> {
        self.alignment()
            .validate(self.identity_cutoff)
            .map_err(|error| JsValue::from_str(&error.to_string()))
    }
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
    start: u64,
    end: u64,
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
    prefilter: BrowserAnalysisSettings,
    genbank_options: GenBankOptions,
}

#[wasm_bindgen]
impl BrowserFileParser {
    #[wasm_bindgen(constructor)]
    pub fn new(prefilter: JsValue) -> Result<BrowserFileParser, JsValue> {
        let prefilter = serde_wasm_bindgen::from_value::<BrowserAnalysisSettings>(prefilter)
            .map_err(|error| JsValue::from_str(&format!("invalid prefilter settings: {error}")))?;
        prefilter.validate()?;
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
    start: u64,
    end: u64,
    genes: Vec<BrowserLayoutGene>,
}

#[derive(Debug, Deserialize, Serialize)]
struct BrowserLayoutGene {
    label: String,
    names: Vec<(String, String)>,
    start: u64,
    end: u64,
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
    pair_count: usize,
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

/// Alignment buffers retained between batches in one browser worker.
#[wasm_bindgen]
#[derive(Default)]
pub struct AlignmentWorker {
    aligner: ProteinBatchAligner,
}

#[wasm_bindgen]
impl AlignmentWorker {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self::default()
    }

    pub fn analyse_batch(
        &mut self,
        proteins: JsValue,
        batch: JsValue,
        pairs: JsValue,
        identity_cutoff: f32,
        progress: &Function,
    ) -> Result<JsValue, JsValue> {
        let proteins = serde_wasm_bindgen::from_value::<Vec<BrowserProtein>>(proteins)
            .map_err(|error| JsValue::from_str(&format!("invalid batch proteins: {error}")))?;
        let batch = serde_wasm_bindgen::from_value::<ProteinBatch>(batch)
            .map_err(|error| JsValue::from_str(&format!("invalid batch: {error}")))?;
        let pairs = serde_wasm_bindgen::from_value::<Vec<ProteinPair>>(pairs)
            .map_err(|error| JsValue::from_str(&format!("invalid batch pairs: {error}")))?;
        if !batch.is_valid(proteins.len(), pairs.len()) {
            return Err(JsValue::from_str(
                "alignment batch indices are out of bounds",
            ));
        }
        let sequences = proteins
            .iter()
            .map(|p| p.translation.as_bytes())
            .collect::<Vec<_>>();
        let result = self.aligner.analyse_pairs(
            &sequences,
            &sequences,
            batch.pairs(&pairs),
            identity_cutoff,
            |processed| report_progress(progress, processed),
        );
        encode_batch(result, &proteins, batch.pair_count())
    }
}

fn encode_batch(
    result: ProteinBatchResult,
    proteins: &[BrowserProtein],
    pair_count: usize,
) -> Result<JsValue, JsValue> {
    let links = result
        .matches
        .into_iter()
        .map(|matched| TileLink {
            query: (&proteins[matched.query_index]).into(),
            target: (&proteins[matched.target_index]).into(),
            identity: matched.identity,
            similarity: matched.similarity,
        })
        .collect();
    serde_wasm_bindgen::to_value(&BrowserTileResult {
        links,
        alignment_count: result.aligned_pairs,
        pair_count,
    })
    .map_err(|error| JsValue::from_str(&format!("could not encode batch links: {error}")))
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
    settings: JsValue,
    use_file_order: bool,
    progress: &Function,
) -> Result<JsValue, JsValue> {
    let settings = serde_wasm_bindgen::from_value::<BrowserAnalysisSettings>(settings)
        .map_err(|error| JsValue::from_str(&format!("invalid analysis settings: {error}")))?;
    settings.validate()?;
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
        settings,
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
    settings: BrowserAnalysisSettings,
    use_file_order: bool,
    progress: &Function,
) -> Result<JsValue, JsValue> {
    let analysis_session = AnalysisSession::new(
        analysis,
        settings.identity_cutoff,
        settings.alignment(),
        group_definitions,
    )
    .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let analysis = &analysis_session.analysis;
    let group_definitions = &analysis_session.group_definitions;
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
    let arranged = analysis
        .to_auto_arranged_plot_data_with_group_definitions_and_progress(
            &order,
            group_definitions,
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
        )
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let cluster_names = analysis
        .clusters
        .iter()
        .map(|cluster| cluster.name.clone())
        .collect();
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
    prefilter: &BrowserAnalysisSettings,
    mut progress: F,
) -> ParsedBrowserFiles
where
    F: FnMut(usize, usize),
{
    let layout = layout_from_clusters(&clusters);
    let prepared = PreparedAlignment::new(
        &clusters,
        0,
        AnalysisOptions {
            identity_cutoff: prefilter.identity_cutoff,
            prefilter: prefilter.enabled.then_some(KmerPrefilter {
                kmer_size: prefilter.kmer_size,
                min_shared_kmers: prefilter.min_shared_kmers,
            }),
        },
        &mut progress,
    );
    let proteins = prepared
        .references
        .iter()
        .map(|r| BrowserProtein {
            cluster: r.cluster,
            locus: r.locus,
            gene: r.gene,
            translation: clusters[r.cluster].loci[r.locus].genes[r.gene]
                .translation
                .clone(),
        })
        .collect();
    ParsedBrowserFiles {
        layout,
        proteins,
        candidate_pairs: prepared.candidates,
        alignment_batches: prepared.batches,
        alignment_counts: prepared.counts,
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
