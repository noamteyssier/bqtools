use std::{
    ffi::CString,
    fs, io,
    os::unix::{ffi::OsStrExt, fs::FileTypeExt},
    path::Path,
};

use anyhow::Result;
use log::{trace, warn};

use super::{PairedChannels, RecordPair};
use crate::cli::FileFormat;

/// The record pairs that get a FIFO and writer thread.
pub fn pairs(paired: bool, channels: PairedChannels) -> &'static [RecordPair] {
    match (paired, channels) {
        (false, _) => &[RecordPair::Unpaired],
        (true, PairedChannels::Both) => &[RecordPair::R1, RecordPair::R2],
        (true, PairedChannels::R1Only) => &[RecordPair::R1],
        (true, PairedChannels::R2Only) => &[RecordPair::R2],
    }
}

/// Creates many FIFOs (named-pipes) at the given basepath.
///
/// For paired files, `channels` controls which channels are created. For
/// unpaired files, `channels` is ignored and a single unlabelled FIFO per
/// thread is created.
///
/// Note: this does not open the FIFOs for writing.
pub fn create_fifos(
    basepath: &str,
    paired: bool,
    num_threads: usize,
    format: FileFormat,
    channels: PairedChannels,
) -> Result<Vec<String>> {
    let pairs = pairs(paired, channels);
    (0..num_threads)
        .flat_map(|idx| pairs.iter().map(move |&pair| (idx, pair)))
        .map(|(idx, pair)| {
            let path = name_fifo(basepath, idx, pair, format);
            trace!("Creating FIFO at path: {path}");
            mkfifo(Path::new(&path))?;
            Ok(path)
        })
        .collect()
}

fn is_fifo(path: &Path) -> bool {
    fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_fifo())
}

/// Creates a FIFO with mode 0600. An existing FIFO is reused; any other
/// existing file is an error.
fn mkfifo(path: &Path) -> io::Result<()> {
    let c_path = CString::new(path.as_os_str().as_bytes())?;
    // SAFETY: `c_path` is a valid NUL-terminated string for the call's duration.
    if unsafe { libc::mkfifo(c_path.as_ptr(), libc::S_IRUSR | libc::S_IWUSR) } == 0 {
        return Ok(());
    }
    match io::Error::last_os_error() {
        err if err.kind() == io::ErrorKind::AlreadyExists && is_fifo(path) => {
            trace!("FIFO already exists at {}, reconnecting...", path.display());
            Ok(())
        }
        err => Err(err),
    }
}

/// RAII guard that unlinks a set of FIFOs when dropped.
///
/// Cleanup is tied to the guard's lifetime rather than the happy path, so the
/// FIFOs are removed from disk on any early return, `?` propagation, or panic
/// (via stack unwinding) — not just on successful completion.
pub struct FifoGuard(pub Vec<String>);

impl Drop for FifoGuard {
    /// Unlink each FIFO. A missing path (`ENOENT`) is treated as success so
    /// cleanup is idempotent; other errors are logged but not propagated, since
    /// this runs during teardown.
    fn drop(&mut self) {
        for path in &self.0 {
            trace!("Closing FIFO at path: {path}");
            match fs::remove_file(path) {
                Ok(()) => {}
                Err(err) if err.kind() == io::ErrorKind::NotFound => {}
                Err(err) => warn!("Failed to unlink FIFO at {path}: {err}"),
            }
        }
    }
}

pub fn name_fifo(basepath: &str, pid: usize, pair: RecordPair, format: FileFormat) -> String {
    let suffix = match pair {
        RecordPair::R1 => "_R1",
        RecordPair::R2 => "_R2",
        RecordPair::Unpaired => "",
    };
    format!("{basepath}_{pid}{suffix}.{}", format.extension())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    fn fifo_at(dir: &tempfile::TempDir) -> std::path::PathBuf {
        dir.path().join("f")
    }

    #[test]
    fn creates_fifo_with_mode_0600() {
        let dir = tempfile::tempdir().unwrap();
        let path = fifo_at(&dir);
        mkfifo(&path).unwrap();
        let meta = fs::metadata(&path).unwrap();
        assert!(meta.file_type().is_fifo());
        assert_eq!(meta.permissions().mode() & 0o777, 0o600);
    }

    #[test]
    fn existing_fifo_is_reused() {
        let dir = tempfile::tempdir().unwrap();
        let base = dir.path().join("base").to_string_lossy().into_owned();
        let args = |b: &str| create_fifos(b, false, 2, FileFormat::Fastq, PairedChannels::Both);
        let first = args(&base).unwrap();
        assert_eq!(args(&base).unwrap(), first);
    }

    #[test]
    fn existing_regular_file_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = fifo_at(&dir);
        fs::write(&path, b"data").unwrap();
        let err = mkfifo(&path).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::AlreadyExists);
        // the file must be left untouched
        assert_eq!(fs::read(&path).unwrap(), b"data");
    }

    #[test]
    fn missing_parent_dir_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let err = mkfifo(&dir.path().join("nope/f")).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
    }

    #[test]
    fn interior_nul_is_an_error() {
        let err = mkfifo(Path::new("a\0b")).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
    }

    #[test]
    fn guard_removes_fifos_and_tolerates_missing() {
        let dir = tempfile::tempdir().unwrap();
        let (a, b) = (fifo_at(&dir), dir.path().join("g"));
        mkfifo(&a).unwrap();
        mkfifo(&b).unwrap();
        fs::remove_file(&b).unwrap(); // already gone: must not panic
        drop(FifoGuard(vec![
            a.to_string_lossy().into_owned(),
            b.to_string_lossy().into_owned(),
        ]));
        assert!(!a.exists());
    }
}
