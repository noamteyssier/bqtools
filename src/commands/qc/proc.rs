use std::{io::Write, path::PathBuf};

use crate::{cli::QcOptions, commands::qc::modules::QcModuleType};

use super::{report::table, QcModule};

use anyhow::{bail, Result};
use binseq::ParallelProcessor;
use log::trace;

use crate::commands::match_output;

const SUMMARY_PATH: &str = "summary.md";

/// TODO: adapter content
#[derive(Clone, Default)]
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
        let mut add_module = |module: QcModuleType| {
            trace!("Loaded: {}", module.desc());
            modules.push(module);
        };

        trace!("Loading QC modules...");
        if !opts.skip_base_qual {
            add_module(QcModuleType::new_base_quality());
        }
        if !opts.skip_seq_qual {
            add_module(QcModuleType::new_seq_quality());
        }
        if !opts.skip_base_content {
            add_module(QcModuleType::new_base_content());
        }
        if !opts.skip_seq_gc {
            add_module(QcModuleType::new_gc_content());
        }
        if !opts.skip_seq_length {
            add_module(QcModuleType::new_seq_length());
        }
        if !opts.skip_dup_levels || !opts.skip_overrepresented {
            add_module(QcModuleType::new_duplication(opts, span_start));
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
        self.modules
            .iter_mut()
            .try_for_each(|m| m.finish(&self.outdir))?;
        self.write_summary()
    }

    /// Writes the high-level `summary.md` report: an overview table followed
    /// by each module's headline stats (the full data still lives in each
    /// module's own TSV).
    fn write_summary(&self) -> Result<()> {
        if !self.outdir.exists() {
            std::fs::create_dir_all(&self.outdir)?;
        }

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

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        self.modules
            .iter_mut()
            .for_each(super::modules::QcModule::sync_batch);
        Ok(())
    }

    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.modules
            .iter_mut()
            .for_each(super::modules::QcModule::sync_final);
        Ok(())
    }
}
