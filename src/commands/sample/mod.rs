use crate::cli::{BinseqMode, OutputBinseqInherited, SampleCommand};
use anyhow::Result;
use binseq::{BinseqReader, BinseqRecord};
use std::io::Write;

use super::{
    decode::{keep, run_with},
    encode::processor::Encoder,
    utils::builder_from_reader,
};

/// Encoder that only keeps the sampled records.
struct Sampler<W: Write + Send> {
    inner: Encoder<W>,
    fraction: f64,
    seed: u64,
}
impl<W: Write + Send> Clone for Sampler<W> {
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            fraction: self.fraction,
            seed: self.seed,
        }
    }
}
impl<W: Write + Send> binseq::ParallelProcessor for Sampler<W> {
    fn process_record<R: BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        if keep(record.index(), self.fraction, self.seed) {
            self.inner.process_record(record)?;
        }
        Ok(())
    }
    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        self.inner.on_batch_complete()
    }
    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.inner.on_thread_complete()
    }
}

fn run_binseq(args: &SampleCommand) -> Result<()> {
    let reader = BinseqReader::new(args.input.path())?;
    let out = OutputBinseqInherited {
        output: args.output.output.clone(),
        pipe: false,
        threads: args.output.threads,
    };
    let writer = builder_from_reader(&reader).build(out.as_writer(args.input.mode()?)?)?;
    let mut proc = Sampler {
        inner: Encoder::new(writer)?,
        fraction: args.sample.fraction,
        seed: args.sample.seed,
    };
    let range = args.input.range(reader.num_records()?)?;
    reader.process_parallel_range(proc.clone(), out.threads(), range)?;
    proc.inner.finish()?;
    Ok(())
}

pub fn run(args: &SampleCommand) -> Result<()> {
    if args
        .output
        .output
        .as_deref()
        .is_some_and(|p| BinseqMode::determine(p).is_ok())
    {
        return run_binseq(args);
    }
    run_with(
        &args.input,
        &args.output,
        Some((args.sample.fraction, args.sample.seed)),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use clap::Parser;
    use itertools::iproduct;
    use tempfile::NamedTempFile;

    use crate::cli::{BinseqMode, FileFormat};
    use crate::testutils::{count_fastx_records, encode, write_fastx};

    fn sample(
        bq_path: &std::path::Path,
        out_path: &std::path::Path,
        fraction: f64,
        seed: u64,
        threads: u32,
    ) -> Result<()> {
        let cmd = crate::cli::SampleCommand::try_parse_from([
            "sample",
            bq_path.to_str().unwrap(),
            "-F",
            &fraction.to_string(),
            "-S",
            &seed.to_string(),
            "-T",
            &threads.to_string(),
            "-o",
            out_path.to_str().unwrap(),
        ])?;
        super::run(&cmd)
    }

    /// Sampling at 0.5 should produce approximately half the records (±20%).
    #[test]
    fn test_sample_half() -> Result<()> {
        let nrec = 1000;
        let fraction = 0.5_f64;
        #[allow(clippy::cast_sign_loss)]
        let expected = (nrec as f64 * fraction) as usize;
        let tolerance = nrec / 5;

        for (mode, fmt) in iproduct!(BinseqMode::enum_iter(), FileFormat::fastx_iter()) {
            let in_tmp = write_fastx().format(fmt).nrec(nrec).call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let out_tmp = NamedTempFile::with_suffix(fmt.fastx_suffix())?;
            sample(bq_tmp.path(), out_tmp.path(), fraction, 42, 1)?;

            let count = count_fastx_records(out_tmp.path())?;
            assert!(
                count.abs_diff(expected) <= tolerance,
                "sample count {count} far from expected {expected} (±{tolerance}) for {mode:?} {fmt:?}"
            );
        }
        Ok(())
    }

    /// Sampling at fraction=1.0 must return all records exactly.
    #[test]
    fn test_sample_fraction_one() -> Result<()> {
        let nrec = 200;
        for mode in BinseqMode::enum_iter() {
            let in_tmp = write_fastx().nrec(nrec).call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let out_tmp = NamedTempFile::with_suffix(".fastq")?;
            sample(bq_tmp.path(), out_tmp.path(), 1.0, 42, 1)?;

            assert_eq!(
                count_fastx_records(out_tmp.path())?,
                nrec,
                "sample fraction=1.0 should return all records for {mode:?}"
            );
        }
        Ok(())
    }

    /// `-o x.<mode>` writes a BINSEQ file with the sampled records.
    #[test]
    fn test_sample_binseq_output() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let in_tmp = write_fastx().nrec(200).call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let out_tmp = NamedTempFile::with_suffix(mode.extension())?;
            sample(bq_tmp.path(), out_tmp.path(), 1.0, 42, 1)?;
            assert_eq!(crate::testutils::count_binseq(out_tmp.path())?, 200);

            sample(bq_tmp.path(), out_tmp.path(), 0.5, 42, 1)?;
            let n = crate::testutils::count_binseq(out_tmp.path())?;
            assert!(n > 50 && n < 150, "{n} {mode:?}");
        }
        Ok(())
    }

    /// Sorted sequence lines of a FASTQ file (output order depends on batch timing).
    fn sorted_seqs(path: &std::path::Path) -> Result<Vec<String>> {
        let mut seqs: Vec<String> = std::fs::read_to_string(path)?
            .lines()
            .skip(1)
            .step_by(4)
            .map(String::from)
            .collect();
        seqs.sort_unstable();
        Ok(seqs)
    }

    /// The same seed must select the same records regardless of `-T`.
    #[test]
    fn test_sample_seed_independent_of_threads() -> Result<()> {
        let in_tmp = write_fastx().nrec(2000).call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let mut results = Vec::new();
        for threads in [1, 4] {
            let out_tmp = NamedTempFile::with_suffix(".fastq")?;
            sample(bq_tmp.path(), out_tmp.path(), 0.3, 7, threads)?;
            results.push(sorted_seqs(out_tmp.path())?);
        }
        assert_eq!(results[0], results[1]);
        Ok(())
    }

    /// Different seeds should (very likely) produce different sample sizes.
    #[test]
    fn test_sample_different_seeds_vary() -> Result<()> {
        let nrec = 1000;
        let in_tmp = write_fastx().nrec(nrec).call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let counts: Vec<usize> = [42_u64, 123, 999]
            .iter()
            .map(|&seed| {
                let out_tmp = NamedTempFile::with_suffix(".fastq")?;
                sample(bq_tmp.path(), out_tmp.path(), 0.5, seed, 1)?;
                count_fastx_records(out_tmp.path())
            })
            .collect::<Result<_>>()?;

        let unique: hashbrown::HashSet<_> = counts.iter().collect();
        assert!(
            unique.len() > 1,
            "all seeds produced the same count — suspicious: {counts:?}"
        );
        Ok(())
    }
}
