use clap::Parser;

use super::{InputBinseq, OutputFile};

/// Decode BINSEQ files to FASTQ, FASTA, or TSV.
///
/// The format is inferred from the `-o` extension, and defaults to TSV on
/// stdout. Records without stored qualities are written with `?` (Phred 30)
/// as FASTQ.
#[derive(Parser, Debug)]
pub struct DecodeCommand {
    #[clap(flatten)]
    pub input: InputBinseq,

    #[clap(flatten)]
    pub output: OutputFile,
}
