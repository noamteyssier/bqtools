use std::{
    fs::File,
    io::BufWriter,
    io::{stderr, Write},
    path::{Path, PathBuf},
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
};

use anyhow::{bail, Result};
use binseq::{BinseqWriter, BinseqWriterBuilder, ParallelProcessor, SequencingRecordBuilder};

use crate::{cli::SplitCommand, commands::split::splitter::Splitter, types::BoxedWriter};

#[derive(Clone)]
pub struct SplitProcessor {
    /// Thread-local matcher
    matcher: Splitter,

    /// Index of the undetermined writer (always last), if active
    undetermined_idx: Option<usize>,

    /// Thread-local writers for the split processor.
    t_writer: Vec<BinseqWriter<Vec<u8>>>,
    t_counts: Vec<usize>,

    /// Global writers for the split processor.
    writer: Arc<Vec<Mutex<BinseqWriter<BoxedWriter>>>>,
    counts: Arc<Vec<AtomicUsize>>,

    /// Name of each output bin, with the undetermined bin last when active.
    aliases: Vec<String>,

    /// Output file path for each writer (parallel to `writer`/`counts`).
    paths: Vec<PathBuf>,
}
impl SplitProcessor {
    pub fn new(
        matcher: Splitter,
        builder: &BinseqWriterBuilder,
        args: &SplitCommand,
    ) -> Result<Self> {
        let output_mode = args.input.mode()?;
        let write_undetermined = !args.split.skip_unmatched;
        let undetermined_basepath = args.split.unmatched_basename.as_str();

        let mut aliases = matcher.aliases().to_vec();
        if write_undetermined && aliases.iter().any(|a| a == undetermined_basepath) {
            bail!(
                "Pattern alias '{undetermined_basepath}' collides with the unmatched output name; set --unmatched-basename to something else"
            );
        }
        let undetermined_idx = write_undetermined.then_some(aliases.len());
        if write_undetermined {
            aliases.push(undetermined_basepath.to_string());
        }

        let mut t_writer = Vec::default();
        let mut writer = Vec::default();
        let mut paths = Vec::default();
        for basename in &aliases {
            let output_path = Path::new(&args.split.basepath).join(format!(
                "{}{}",
                basename,
                output_mode.extension()
            ));
            let output_handle: BoxedWriter = Box::new(BufWriter::new(File::create(&output_path)?));

            let gw = builder.clone().build(output_handle)?;
            t_writer.push(gw.new_headless_buffer()?);
            writer.push(Mutex::new(gw));
            paths.push(output_path);
        }

        Ok(Self {
            matcher,
            undetermined_idx,
            t_counts: vec![0; writer.len()],
            counts: Arc::new((0..writer.len()).map(|_| AtomicUsize::new(0)).collect()),
            t_writer,
            writer: Arc::new(writer),
            aliases,
            paths,
        })
    }

    pub fn finish(&mut self) -> binseq::Result<()> {
        self.writer
            .iter()
            .try_for_each(|w| w.lock().unwrap().finish())
    }

    /// Removes any output files that received fewer than `min_records` records.
    ///
    /// Must be called after [`finish`](Self::finish) so all writers are flushed.
    /// Returns the number of files removed.
    pub fn prune_below(&self, min_records: usize) -> Result<usize> {
        let mut removed = 0;
        for (path, count) in self
            .paths
            .iter()
            .zip(self.counts.iter().map(|c| c.load(Ordering::Relaxed)))
        {
            if count < min_records {
                log::debug!(
                    "Removing {} ({count} records, below threshold of {min_records})",
                    path.display(),
                );
                std::fs::remove_file(path)?;
                removed += 1;
            }
        }
        Ok(removed)
    }

    pub fn pprint_counts(&self) -> Result<()> {
        let mut handle = stderr();
        self.aliases
            .iter()
            .zip(self.counts.iter().map(|c| c.load(Ordering::Relaxed)))
            .try_for_each(|(alias, count)| writeln!(&mut handle, "{alias}\t{count}"))?;
        handle.flush().map_err(Into::into)
    }

    /// Moves each thread-local buffer into its global writer with `ingest`.
    fn ingest_all(
        &mut self,
        ingest: impl Fn(
            &mut BinseqWriter<BoxedWriter>,
            &mut BinseqWriter<Vec<u8>>,
        ) -> binseq::Result<()>,
    ) -> binseq::Result<()> {
        self.writer
            .iter()
            .zip(self.t_writer.iter_mut())
            .try_for_each(|(global, local)| ingest(&mut global.lock().unwrap(), local))
    }
}
impl ParallelProcessor for SplitProcessor {
    fn process_record<R: binseq::prelude::BinseqRecord>(
        &mut self,
        record: R,
    ) -> binseq::Result<()> {
        let paired = record.is_paired();
        let has_qual = record.has_quality();
        let (sseq, xseq) = (record.sseq(), record.xseq());
        let rec = SequencingRecordBuilder::default()
            .s_seq(sseq)
            .opt_s_qual(has_qual.then(|| record.squal()))
            .s_header(record.sheader())
            .opt_x_seq(paired.then_some(xseq))
            .opt_x_qual((paired && has_qual).then(|| record.xqual()))
            .opt_x_header(paired.then(|| record.xheader()))
            .build()?;

        // matched bin, else the undetermined bin (if active)
        if let Some(idx) = self.matcher.split_idx(sseq, xseq).or(self.undetermined_idx) {
            self.t_writer[idx].push(rec)?;
            self.t_counts[idx] += 1;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        for (global, local) in self.counts.iter().zip(self.t_counts.iter_mut()) {
            global.fetch_add(std::mem::take(local), Ordering::Relaxed);
        }
        self.ingest_all(BinseqWriter::ingest_completed)
    }

    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.ingest_all(BinseqWriter::ingest)
    }
}
