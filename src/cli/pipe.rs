use clap::{
    builder::{PossibleValue, PossibleValuesParser, TypedValueParser},
    Parser,
};

use crate::cli::FileFormat;

use super::InputBinseq;

/// Split a BINSEQ file into multiple named pipes (FIFOs) for legacy tools (Unix only).
///
/// Without `-x`/`-X`, the command blocks until every FIFO has been opened and
/// fully read. FIFOs are removed on exit.
#[derive(Parser, Debug)]
pub struct PipeCommand {
    #[clap(flatten)]
    pub input: InputBinseq,

    #[clap(flatten)]
    pub pipe: PipeOptions,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "PIPE OPTIONS")]
pub struct PipeOptions {
    /// Number of FIFOs to create (0 = number of CPUs; capped at CPU count)
    ///
    /// For paired input this is split into p/2 R1/R2 pairs, so `{n}` ranges
    /// over 0..p/2.
    #[clap(short = 'p', long, default_value = "0")]
    num_pipes: usize,

    /// Record format written to each FIFO
    #[clap(short, long, default_value = "q", value_parser = parse_pipe_format())]
    format: FileFormat,

    /// Base path for the FIFOs
    ///
    /// FIFOs are named `{basepath}_{n}.{fa|fq}` (single-end) or
    /// `{basepath}_{n}_R1.{ext}` / `{basepath}_{n}_R2.{ext}` (paired), with `n`
    /// starting at 0. An existing FIFO at that path is reused.
    #[clap(short, long, default_value = "bqtools_fifo")]
    basepath: String,

    /// Execute a shell command once per pipe, substituting FIFO paths.
    ///
    /// Use `{}` for the FIFO path (single-end), or `{R1}` / `{R2}` for the
    /// respective paths (paired-end). Referencing only one of `{R1}` / `{R2}`
    /// processes just that mate — the other channel's FIFOs are never created.
    /// `{n}` expands to the pipe index, useful for per-shard output paths.
    /// Commands run via `sh -c`; the template must contain the placeholder for
    /// the input type, and the run exits non-zero if any command fails.
    /// Mutually exclusive with `--exec-batch`.
    #[clap(short = 'x', long, conflicts_with = "exec_batch")]
    exec: Option<String>,

    /// Execute a single shell command with all FIFO paths substituted.
    ///
    /// `{}` (single-end) or `{R1}` / `{R2}` (paired-end) each expand to a
    /// space-joined list of every matching FIFO path. Writing `{R1} {R2}`
    /// adjacently interleaves the paths as pairs (`r1_0` `r2_0` `r1_1` `r2_1` …) so
    /// positional-argument tools receive each pair together.
    /// `{n}` is not expanded. Mutually exclusive with `--exec`.
    #[clap(short = 'X', long, conflicts_with = "exec")]
    exec_batch: Option<String>,
}

/// Pipes only support FASTA (`a`) and FASTQ (`q`).
fn parse_pipe_format() -> impl TypedValueParser<Value = FileFormat> {
    PossibleValuesParser::new([
        PossibleValue::new("a").help("FASTA file format"),
        PossibleValue::new("q").help("FASTQ file format"),
    ])
    .map(|s| {
        if s == "a" {
            FileFormat::Fasta
        } else {
            FileFormat::Fastq
        }
    })
}

impl PipeCommand {
    pub fn format(&self) -> FileFormat {
        self.pipe.format
    }
    pub fn num_pipes(&self) -> usize {
        match self.pipe.num_pipes {
            0 => num_cpus::get(),
            n => n.min(num_cpus::get()),
        }
    }
    pub fn basepath(&self) -> &str {
        &self.pipe.basepath
    }
    pub fn exec(&self) -> Option<&str> {
        self.pipe.exec.as_deref()
    }
    pub fn exec_batch(&self) -> Option<&str> {
        self.pipe.exec_batch.as_deref()
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::PipeCommand;
    use crate::cli::FileFormat;

    #[test]
    fn test_format_defaults_to_fastq_and_rejects_non_fastx() {
        let cmd = PipeCommand::try_parse_from(["pipe", "x.cbq"]).unwrap();
        assert_eq!(cmd.format(), FileFormat::Fastq);
        let cmd = PipeCommand::try_parse_from(["pipe", "x.cbq", "-f", "a"]).unwrap();
        assert_eq!(cmd.format(), FileFormat::Fasta);
        for bad in ["b", "t"] {
            assert!(PipeCommand::try_parse_from(["pipe", "x.cbq", "-f", bad]).is_err());
        }
    }
}
