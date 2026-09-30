use clap::Parser;

use super::{InputBinseq, OutputFile};

/// Subsample a BINSEQ file and output to FASTQ, FASTA, TSV, or BINSEQ
///
/// An `-o` ending in `.bq/.vbq/.cbq` writes a BINSEQ file in the input's mode
/// and settings (the extension must match the input).
///
/// Output defaults to TSV on stdout; use `-o reads.fastq[.gz]` or `-f q` for
/// FASTQ. Record order in the output is not preserved.
#[derive(Parser)]
pub struct SampleCommand {
    #[clap(flatten)]
    pub input: InputBinseq,

    #[clap(flatten)]
    pub output: OutputFile,

    #[clap(flatten)]
    pub sample: SampleArgs,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "SAMPLE OPTIONS")]
pub struct SampleArgs {
    /// Fraction of reads to keep, in (0, 1]
    ///
    /// Each record is kept independently with this probability, so the output
    /// count is approximate. Applied within `--span` when given.
    #[clap(short = 'F', long, value_parser = parse_fraction, required_unless_present = "num", conflicts_with = "num")]
    pub fraction: Option<f64>,

    /// Exact number of reads to keep
    ///
    /// Selects exactly this many records uniformly at random, drawn from
    /// `--span` when given. If it exceeds the available records, all are kept.
    #[clap(short = 'n', long)]
    pub num: Option<usize>,

    /// Seed for random sampling
    ///
    /// The same seed selects the same records regardless of thread count.
    #[clap(short = 'S', long, default_value = "42")]
    pub seed: u64,
}

fn parse_fraction(s: &str) -> Result<f64, String> {
    match s.parse::<f64>() {
        Ok(f) if f > 0.0 && f <= 1.0 => Ok(f),
        Ok(_) => Err("Fraction must be between 0 and 1".into()),
        Err(e) => Err(e.to_string()),
    }
}
