use std::{
    collections::HashSet,
    fs, io,
    path::{Path, PathBuf},
};

use clinker_core::{Cluster, InputFile, ParseOptions, parse_input_files_with_options};

pub struct InputSource {
    pub path: PathBuf,
    companion: Option<PathBuf>,
}

impl InputSource {
    pub fn paths(&self) -> impl Iterator<Item = &Path> {
        std::iter::once(self.path.as_path()).chain(self.companion.as_deref())
    }
}

pub fn discover(inputs: &[PathBuf]) -> io::Result<Vec<InputSource>> {
    let mut paths = Vec::new();
    for input in inputs {
        collect(input, true, &mut paths)?;
    }
    let mut seen = HashSet::new();
    let mut unique = Vec::new();
    for path in paths {
        if seen.insert(fs::canonicalize(&path)?) {
            unique.push(path);
        }
    }
    let explicit_fasta = unique
        .iter()
        .filter(|path| is_fasta(path))
        .map(|path| Ok(fs::canonicalize(parent_directory(path))?.join(path.file_name().unwrap())))
        .collect::<io::Result<Vec<_>>>()?;
    let mut used_fasta = HashSet::new();
    let mut sources = Vec::new();
    for path in unique.iter().filter(|path| !is_fasta(path)) {
        let companion = if is_gff(path) {
            companion(path, &explicit_fasta)?
        } else {
            None
        };
        if let Some(path) = &companion {
            used_fasta.insert(fs::canonicalize(path)?);
        }
        sources.push(InputSource {
            path: path.clone(),
            companion,
        });
    }
    for path in explicit_fasta {
        if !used_fasta.contains(&fs::canonicalize(&path)?) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{} has no GFF3/GFF input with the same directory and stem",
                    path.display()
                ),
            ));
        }
    }
    Ok(sources)
}

pub fn parse(
    sources: &[InputSource],
    options: &ParseOptions,
) -> Result<Vec<Cluster>, Box<dyn std::error::Error>> {
    let mut clusters = Vec::new();
    for source in sources {
        let files = source
            .paths()
            .map(|path| {
                let name = path
                    .file_name()
                    .unwrap_or(path.as_os_str())
                    .to_string_lossy()
                    .into_owned();
                fs::read(path).map(|bytes| (name, bytes)).map_err(|error| {
                    io::Error::new(error.kind(), format!("{}: {error}", path.display()))
                })
            })
            .collect::<io::Result<Vec<_>>>()?;
        let inputs = files
            .iter()
            .map(|(name, bytes)| InputFile { name, bytes })
            .collect::<Vec<_>>();
        // Pairing is resolved by path; display filenames are scoped to this source.
        clusters.extend(parse_input_files_with_options(&inputs, options)?);
    }
    Ok(clusters)
}

fn collect(path: &Path, explicit: bool, paths: &mut Vec<PathBuf>) -> io::Result<()> {
    let metadata = fs::symlink_metadata(path)?;
    let file_type = metadata.file_type();
    if file_type.is_symlink() && path.is_dir() {
        return Ok(());
    }
    if file_type.is_dir() {
        let mut entries = fs::read_dir(path)?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<io::Result<Vec<_>>>()?;
        entries.sort();
        for entry in entries {
            collect(&entry, false, paths)?;
        }
    } else if path.is_file() && (is_primary(path) || (explicit && is_fasta(path))) {
        paths.push(path.to_owned());
    } else if explicit {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "unsupported input {}; expected GenBank, GFF3/GFF, or FASTA",
                path.display()
            ),
        ));
    }
    Ok(())
}

fn companion(gff: &Path, explicit_fasta: &[PathBuf]) -> io::Result<Option<PathBuf>> {
    let parent = fs::canonicalize(parent_directory(gff))?;
    let matches = |path: &Path| {
        is_fasta(path)
            && path
                .file_stem()
                .zip(gff.file_stem())
                .is_some_and(|(one, two)| {
                    one.to_string_lossy()
                        .eq_ignore_ascii_case(&two.to_string_lossy())
                })
    };
    let mut candidates = explicit_fasta
        .iter()
        .filter(|path| path.parent() == Some(parent.as_path()) && matches(path))
        .cloned()
        .collect::<Vec<_>>();
    if candidates.is_empty() {
        let mut seen = HashSet::new();
        for entry in fs::read_dir(&parent)? {
            let path = entry?.path();
            if path.is_file() && matches(&path) && seen.insert(fs::canonicalize(&path)?) {
                candidates.push(path);
            }
        }
    }
    candidates.sort();
    match candidates.as_slice() {
        [] => Ok(None),
        [path] => Ok(Some(path.clone())),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "ambiguous FASTA companions for {}: {}; select one explicitly",
                gff.display(),
                candidates
                    .iter()
                    .map(|path| path.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        )),
    }
}

fn parent_directory(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn extension(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase)
}

fn is_primary(path: &Path) -> bool {
    matches!(
        extension(path).as_deref(),
        Some("gb" | "gbk" | "genbank" | "gbf" | "gbff" | "gff" | "gff3")
    )
}

fn is_gff(path: &Path) -> bool {
    matches!(extension(path).as_deref(), Some("gff" | "gff3"))
}

fn is_fasta(path: &Path) -> bool {
    matches!(
        extension(path).as_deref(),
        Some("fa" | "fsa" | "fna" | "fasta" | "faa")
    )
}
