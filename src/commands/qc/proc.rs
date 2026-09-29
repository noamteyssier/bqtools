use std::{io::Write, path::PathBuf};

use crate::{
    cli::QcOptions,
    commands::qc::{
        base_content::PerBaseSequenceContent, base_quality::PerBaseSequenceQuality,
        dup_levels::SequenceDuplicationLevels, gc_content::PerSequenceGcContent,
        modules::QcModuleType, seq_length::SequenceLengthDistribution,
        seq_quality::PerSequenceQuality,
    },
};

use super::report::table;

use anyhow::{bail, Result};
use binseq::ParallelProcessor;
use log::trace;

use crate::commands::match_output;

const SUMMARY_PATH: &str = "summary.md";

/// TODO: adapter content
#[derive(Clone)]
pub struct QcProcessor {
    outdir: PathBuf,
    modules: Vec<QcModuleType>,
    input_path: String,
    num_records: usize,
    paired: bool,
}
impl QcProcessor {
    /// `span_start` is the first record index processed (non-zero with `--span`).
    pub fn new(
        opts: &QcOptions,
        span_start: usize,
        input_path: String,
        num_records: usize,
        paired: bool,
    ) -> Result<Self> {
        let mut modules = Vec::default();
        trace!("Loading QC modules...");
        if !opts.skip_base_qual {
            modules.push(QcModuleType::BaseQuality(PerBaseSequenceQuality::default()));
        }
        if !opts.skip_seq_qual {
            modules.push(QcModuleType::SeqQuality(PerSequenceQuality::default()));
        }
        if !opts.skip_base_content {
            modules.push(QcModuleType::BaseContent(PerBaseSequenceContent::default()));
        }
        if !opts.skip_seq_gc {
            modules.push(QcModuleType::GcContent(PerSequenceGcContent::default()));
        }
        if !opts.skip_seq_length {
            modules.push(QcModuleType::SeqLength(
                SequenceLengthDistribution::default(),
            ));
        }
        if !opts.skip_dup_levels || !opts.skip_overrepresented {
            modules.push(QcModuleType::Duplication(SequenceDuplicationLevels::new(
                opts, span_start,
            )));
        }
        trace!("{} modules loaded", modules.len());

        if modules.is_empty() {
            bail!("Must provide at least one QC module to process")
        }
        Ok(Self {
            outdir: PathBuf::from(&opts.outdir),
            modules,
            input_path,
            num_records,
            paired,
        })
    }

    pub fn finish(&mut self) -> Result<()> {
        if !self.outdir.exists() {
            std::fs::create_dir_all(&self.outdir)?;
        }
        self.modules
            .iter_mut()
            .try_for_each(|m| m.finish(&self.outdir))?;
        self.write_summary()
    }

    /// Writes the high-level `summary.md` report: an overview table followed
    /// by each module's headline stats (the full data still lives in each
    /// module's own TSV).
    fn write_summary(&self) -> Result<()> {
        let mut handle = match_output(Some(self.outdir.join(SUMMARY_PATH)))?;

        writeln!(handle, "# BQtools QC Report\n")?;
        write!(
            handle,
            "{}",
            table(
                &["Metric", "Value"],
                &[
                    vec!["Input".into(), self.input_path.clone()],
                    vec!["Reads".into(), self.num_records.to_string()],
                    vec!["Paired".into(), self.paired.to_string()],
                ],
            )
        )?;
        writeln!(handle)?;

        for module in &self.modules {
            let section = module.summarize();
            if !section.is_empty() {
                writeln!(handle, "{section}")?;
            }
        }

        Ok(())
    }
}
impl ParallelProcessor for QcProcessor {
    fn process_record<R: binseq::prelude::BinseqRecord>(
        &mut self,
        record: R,
    ) -> binseq::Result<()> {
        self.modules.iter_mut().for_each(|m| m.push(&record));
        Ok(())
    }

    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.modules.iter_mut().for_each(QcModuleType::sync_final);
        Ok(())
    }
}
