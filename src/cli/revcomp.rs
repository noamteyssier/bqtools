use clap::Parser;

use super::{InputBinseq, Mate, OutputBinseqInherited};

/// Reverse complement the sequences in a BINSEQ file.
///
/// The output keeps the input's BINSEQ mode and encoding settings. Records are
/// written as parallel batches complete, so output order is not preserved.
#[derive(Parser, Debug)]
pub struct RevcompCommand {
    #[clap(flatten)]
    pub input: InputBinseq,

    #[clap(flatten)]
    pub output: OutputBinseqInherited,

    /// Which mate(s) to reverse complement
    ///
    /// Only relevant for paired BINSEQ files. Defaults to reverse
    /// complementing both mates; ignored (with a warning) on single-end
    /// files, where the single read is always reverse complemented.
    // `-M` rather than `-m` for consistency with `verify --mate`.
    #[clap(short = 'M', long, default_value = "both")]
    pub mate: Mate,
}
