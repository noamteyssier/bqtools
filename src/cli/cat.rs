use clap::Parser;

use super::{MultiInputBinseq, OutputBinseqInherited};

#[derive(Parser, Debug)]
/// Concatenate BINSEQ files.
///
/// All inputs must be the same BINSEQ mode with identical file headers; the
/// output inherits that mode and its settings. bq inputs are byte-copied in
/// order, while vbq/cbq records are re-encoded in parallel, so their order
/// is not preserved.
pub struct CatCommand {
    #[clap(flatten)]
    pub input: MultiInputBinseq,

    #[clap(flatten)]
    pub output: OutputBinseqInherited,
}
