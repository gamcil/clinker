use std::{
    fs::{self, File, OpenOptions},
    io::{self, Write},
    path::{Path, PathBuf},
};

use clinker_core::Analysis;
use tempfile::NamedTempFile;

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

pub fn validate<'a>(
    destinations: impl IntoIterator<Item = (&'a Path, bool)>,
    inputs: impl IntoIterator<Item = &'a Path>,
) -> io::Result<()> {
    let protected = inputs
        .into_iter()
        .map(fs::canonicalize)
        .collect::<io::Result<Vec<_>>>()?;
    let mut directories = Vec::<(PathBuf, tempfile::TempDir)>::new();
    for (path, overwrite) in destinations {
        let destination = resolve(path)?;
        if destination.try_exists()? {
            for input in &protected {
                if same_file::is_same_file(input, &destination)? {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{} is also an input file", path.display()),
                    ));
                }
            }
        }
        check_destination(&destination, overwrite)?;
        let parent = parent_directory(&destination);
        let mut probe_index = 0;
        while probe_index < directories.len()
            && !same_file::is_same_file(parent, &directories[probe_index].0)?
        {
            probe_index += 1;
        }
        if probe_index == directories.len() {
            directories.push((
                parent.to_owned(),
                tempfile::Builder::new()
                    .prefix(".clinker-check-")
                    .tempdir_in(parent)?,
            ));
        }
        let directory = &directories[probe_index].1;
        // Probe names on the destination filesystem to catch case and Unicode aliases.
        let probe = directory.path().join(destination.file_name().unwrap());
        OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(probe)
            .map_err(|error| {
                if error.kind() == io::ErrorKind::AlreadyExists {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("{} is used by more than one output", path.display()),
                    )
                } else {
                    io::Error::new(
                        error.kind(),
                        format!("invalid output {}: {error}", path.display()),
                    )
                }
            })?;
    }
    Ok(())
}

pub struct PendingOutput {
    path: PathBuf,
    temporary: NamedTempFile,
    overwrite: bool,
}

impl PendingOutput {
    pub fn stage(path: &Path, contents: &[u8], overwrite: bool) -> io::Result<Self> {
        Self::stage_with(path, overwrite, |file| file.write_all(contents)).map_err(|error| {
            io::Error::new(
                error.kind(),
                format!("could not prepare {}: {error}", path.display()),
            )
        })
    }

    fn stage_with(
        path: &Path,
        overwrite: bool,
        write: impl FnOnce(&mut File) -> io::Result<()>,
    ) -> io::Result<Self> {
        let path = resolve(path)?;
        check_destination(&path, overwrite)?;
        let mut temporary = NamedTempFile::new_in(parent_directory(&path))?;
        write(temporary.as_file_mut())?;
        if let Ok(metadata) = fs::metadata(&path) {
            temporary
                .as_file()
                .set_permissions(metadata.permissions())?;
        }
        temporary.as_file().sync_all()?;
        Ok(Self {
            path,
            temporary,
            overwrite,
        })
    }

    pub fn commit(self) -> io::Result<()> {
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

pub fn distance_matrix(analysis: &Analysis) -> Result<Vec<u8>, csv::Error> {
    let mut writer = csv::Writer::from_writer(Vec::new());
    writer.write_record(
        std::iter::once("").chain(
            analysis
                .clusters
                .iter()
                .map(|cluster| cluster.name.as_str()),
        ),
    )?;
    for (cluster, row) in analysis
        .clusters
        .iter()
        .zip(analysis.synteny_distance_matrix())
    {
        writer.write_record(
            std::iter::once(cluster.name.clone()).chain(row.iter().map(ToString::to_string)),
        )?;
    }
    writer
        .into_inner()
        .map_err(|error| error.into_error().into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_round_trips_cluster_names_with_commas_quotes_and_newlines() {
        let names = ["name, with comma", "a \"quoted\" name\nsecond line"];
        let analysis = Analysis {
            clusters: names
                .iter()
                .map(|name| clinker_core::Cluster {
                    name: (*name).to_owned(),
                    source_file: None,
                    organism: None,
                    loci: Vec::new(),
                })
                .collect(),
            links: Vec::new(),
        };
        let bytes = distance_matrix(&analysis).unwrap();
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(false)
            .from_reader(bytes.as_slice());
        let rows = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(
            rows[0],
            csv::StringRecord::from(vec!["", names[0], names[1]])
        );
        assert_eq!(rows[1], csv::StringRecord::from(vec![names[0], "0", "1"]));
        assert_eq!(rows[2], csv::StringRecord::from(vec![names[1], "1", "0"]));
    }

    #[test]
    fn failed_staging_preserves_the_previous_file_and_cleans_up() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("session.json");
        fs::write(&path, b"previous").unwrap();
        let result = PendingOutput::stage_with(&path, true, |file| {
            file.write_all(b"partial")?;
            Err(io::Error::other("simulated write failure"))
        });
        assert!(result.is_err());
        assert_eq!(fs::read(&path).unwrap(), b"previous");
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 1);
    }

    #[test]
    fn dropping_staged_output_leaves_the_destination_unchanged() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("session.json");
        fs::write(&path, b"previous").unwrap();
        let pending = PendingOutput::stage(&path, b"replacement", true).unwrap();
        assert_eq!(fs::read(&path).unwrap(), b"previous");
        drop(pending);
        assert_eq!(fs::read(&path).unwrap(), b"previous");
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 1);
    }

    #[test]
    fn exclusive_commit_preserves_a_file_created_after_staging() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("figure.json");
        let pending = PendingOutput::stage(&path, b"new output", false).unwrap();
        fs::write(&path, b"other writer").unwrap();
        assert_eq!(
            pending.commit().unwrap_err().kind(),
            io::ErrorKind::AlreadyExists
        );
        assert_eq!(fs::read(&path).unwrap(), b"other writer");
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 1);
    }

    #[test]
    fn commit_replaces_the_previous_file_completely() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("session.json");
        fs::write(&path, b"previous longer content").unwrap();
        PendingOutput::stage(&path, b"new", true)
            .unwrap()
            .commit()
            .unwrap();
        assert_eq!(fs::read(&path).unwrap(), b"new");
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 1);
    }
}
