use clap::Parser;

use super::{MultiInputBinseq, OutputBinseqInherited};

#[derive(Parser, Debug)]
/// Concatenate BINSEQ files.
pub struct CatCommand {
    #[clap(flatten)]
    pub input: MultiInputBinseq,

    #[clap(flatten)]
    pub output: OutputBinseqInherited,
}
