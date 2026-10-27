use std::{
    fs::{self, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
};

use clinker_core::{Analysis, ClusterPairSimilarity};
use tempfile::NamedTempFile;

use crate::{Args, inputs::InputSource};

pub enum PlotFormat {
    Json,
    Html,
}

impl PlotFormat {
    pub fn from_path(path: &Path) -> io::Result<Self> {
        match path
            .extension()
            .and_then(|extension| extension.to_str())
            .map(str::to_ascii_lowercase)
            .as_deref()
        {
            Some("json") => Ok(Self::Json),
            Some("html" | "htm") => Ok(Self::Html),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "--plot output must end in .json, .html, or .htm",
            )),
        }
    }
}

pub struct Outputs {
    pub summary: Option<OutputFile>,
    pub plot: Option<OutputFile>,
    pub project: Option<OutputFile>,
    pub matrix: Option<OutputFile>,
    session: Option<OutputFile>,
    pub plot_format: Option<PlotFormat>,
}

impl Outputs {
    pub fn prepare(args: &Args, sources: &[InputSource], session_exists: bool) -> io::Result<Self> {
        let plot_path = args.plot.as_ref().and_then(Option::as_deref);
        let plot_format = plot_path.map(PlotFormat::from_path).transpose()?;
        let save_session = args.session.is_some()
            && (!session_exists
                || args.reanalyse
                || !sources.is_empty()
                || args.gene_functions.is_some());
        let protected = sources
            .iter()
            .flat_map(InputSource::paths)
            .chain(args.gene_functions.as_deref())
            .chain(args.colour_map.as_deref())
            .chain(args.session.as_deref().filter(|_| !save_session))
            .collect::<Vec<_>>();
        let mut directories = Vec::<(PathBuf, tempfile::TempDir)>::new();
        let mut prepare = |path: Option<&Path>, overwrite| -> io::Result<Option<OutputFile>> {
            let Some(path) = path else { return Ok(None) };
            let path = resolve(path)?;
            if path.try_exists()? {
                for input in &protected {
                    if same_file::is_same_file(input, &path)? {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!("{} is also an input file", path.display()),
                        ));
                    }
                }
            }
            check_destination(&path, overwrite)?;
            let parent = parent_directory(&path);
            let mut index = 0;
            while index < directories.len()
                && !same_file::is_same_file(parent, &directories[index].0)?
            {
                index += 1;
            }
            if index == directories.len() {
                directories.push((parent.to_owned(), tempfile::tempdir_in(parent)?));
            }
            // Use the actual filesystem's filename rules, including case and Unicode aliases.
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(directories[index].1.path().join(path.file_name().unwrap()))
                .map_err(|error| {
                    io::Error::new(
                        error.kind(),
                        format!("conflicting or invalid output {}: {error}", path.display()),
                    )
                })?;
            Ok(Some(OutputFile::new(path, overwrite)?))
        };
        Ok(Self {
            summary: prepare(args.output.as_deref(), args.force)?,
            plot: prepare(plot_path, args.force)?,
            project: prepare(args.project.as_deref(), args.force)?,
            matrix: prepare(args.matrix_out.as_deref(), args.force)?,
            session: prepare(
                args.session.as_deref().filter(|_| save_session),
                session_exists,
            )?,
            plot_format,
        })
    }

    pub fn commit(
        mut self,
        session: &clinker_core::AnalysisSession,
    ) -> Result<(), Box<dyn std::error::Error>> {
        if let Some(output) = &mut self.session {
            output.write(&session.to_json()?)?;
        }
        for output in [
            self.summary,
            self.plot,
            self.project,
            self.matrix,
            self.session,
        ]
        .into_iter()
        .flatten()
        {
            output.commit()?;
        }
        Ok(())
    }
}

pub struct OutputFile {
    path: PathBuf,
    temporary: NamedTempFile,
    overwrite: bool,
}

impl OutputFile {
    fn new(path: PathBuf, overwrite: bool) -> io::Result<Self> {
        let temporary = NamedTempFile::new_in(parent_directory(&path))?;
        if let Ok(metadata) = fs::metadata(&path) {
            temporary
                .as_file()
                .set_permissions(metadata.permissions())?;
        }
        Ok(Self {
            path,
            temporary,
            overwrite,
        })
    }

    pub fn write(&mut self, contents: &[u8]) -> io::Result<()> {
        self.temporary.write_all(contents).map_err(|error| {
            io::Error::new(
                error.kind(),
                format!("could not prepare {}: {error}", self.path.display()),
            )
        })
    }

    fn commit(self) -> io::Result<()> {
        self.temporary.as_file().sync_all()?;
        let result = if self.overwrite {
            self.temporary.persist(&self.path)
        } else {
            self.temporary.persist_noclobber(&self.path)
        };
        result.map(|_| ()).map_err(|error| {
            io::Error::new(
                error.error.kind(),
                format!("could not save {}: {}", self.path.display(), error.error),
            )
        })
    }
}

fn resolve(path: &Path) -> io::Result<PathBuf> {
    let name = path.file_name().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{} is not an output filename", path.display()),
        )
    })?;
    let parent = parent_directory(path);
    Ok(fs::canonicalize(parent)
        .map_err(|error| {
            io::Error::new(
                error.kind(),
                format!("output directory {}: {error}", parent.display()),
            )
        })?
        .join(name))
}

fn parent_directory(path: &Path) -> &Path {
    path.parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn check_destination(path: &Path, overwrite: bool) -> io::Result<()> {
    let metadata = match fs::symlink_metadata(path) {
        Ok(metadata) => metadata,
        Err(error) if error.kind() == io::ErrorKind::NotFound => return Ok(()),
        Err(error) => return Err(error),
    };
    if !metadata.is_file() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "{} must be a regular file, not a directory or symlink",
                path.display()
            ),
        ));
    }
    if !overwrite {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            format!(
                "{} already exists; pass --force to overwrite it",
                path.display()
            ),
        ));
    }
    if metadata.permissions().readonly() {
        return Err(io::Error::new(
            io::ErrorKind::PermissionDenied,
            format!("{} is read-only", path.display()),
        ));
    }
    Ok(())
}

pub fn distance_matrix(
    analysis: &Analysis,
    similarities: &[Vec<ClusterPairSimilarity>],
) -> Result<Vec<u8>, csv::Error> {
    let mut writer = csv::Writer::from_writer(Vec::new());
    writer.write_record(
        std::iter::once("").chain(
            analysis
                .clusters
                .iter()
                .map(|cluster| cluster.name.as_str()),
        ),
    )?;
    for (query, (cluster, row)) in analysis.clusters.iter().zip(similarities).enumerate() {
        writer.write_record(std::iter::once(cluster.name.clone()).chain(
            row.iter().enumerate().map(|(target, pair)| {
                let distance = if query == target {
                    0.0
                } else {
                    1.0 - pair.similarity
                };
                distance.to_string()
            }),
        ))?;
    }
    writer
        .into_inner()
        .map_err(|error| error.into_error().into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_failed_write_leaves_the_previous_file_intact() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("session.json");
        fs::write(&path, b"previous session").unwrap();
        let mut output = OutputFile::new(path.clone(), true).unwrap();
        output.write(b"partial output").unwrap();
        // Reopen the staged file read-only to force a real write error.
        let (file, temporary_path) = output.temporary.into_parts();
        drop(file);
        let read_only = fs::File::open(&temporary_path).unwrap();
        output.temporary = NamedTempFile::from_parts(read_only, temporary_path);
        assert!(output.write(b"remaining output").is_err());
        drop(output);
        assert_eq!(fs::read(path).unwrap(), b"previous session");
    }

    #[test]
    fn exclusive_commit_preserves_a_file_created_during_analysis() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("figure.json");
        let mut output = OutputFile::new(path.clone(), false).unwrap();
        output.write(b"new output").unwrap();
        fs::write(&path, b"other writer").unwrap();
        assert_eq!(
            output.commit().unwrap_err().kind(),
            io::ErrorKind::AlreadyExists
        );
        assert_eq!(fs::read(&path).unwrap(), b"other writer");
    }
}
