use nix::errno::Errno;
use std::path::Path;

use anyhow::Result;
use log::{trace, warn};
use nix::sys::stat;
use nix::unistd;

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
            unistd::mkfifo(Path::new(&path), stat::Mode::S_IRUSR | stat::Mode::S_IWUSR).or_else(
                |err| match err {
                    Errno::EEXIST => {
                        trace!("FIFO already exists at {path}, reconnecting...");
                        Ok(())
                    }
                    err => Err(err),
                },
            )?;
            Ok(path)
        })
        .collect()
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
            match unistd::unlink(Path::new(path)) {
                Ok(()) | Err(Errno::ENOENT) => {}
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
