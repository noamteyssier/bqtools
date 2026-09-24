use clap::Parser;

#[derive(Parser, Debug)]
/// Show metadata and statistics for one or more BINSEQ files (BQ/VBQ/CBQ).
///
/// Unreadable inputs are skipped with a warning; the command fails only if
/// none can be read.
pub struct InfoCommand {
    /// One or more BINSEQ files (.bq/.vbq/.cbq) to inspect
    #[clap(num_args=1.., required=true)]
    pub input: Vec<String>,

    #[clap(flatten)]
    pub opts: InfoOpts,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "INFO OPTIONS")]
#[allow(clippy::struct_excessive_bools)]
pub struct InfoOpts {
    /// Print only the record count of each input, as `<count>\t<path>`
    #[clap(short, long, conflicts_with_all=["json", "show_index", "show_headers"])]
    pub num: bool,

    /// Print file metadata as a JSON array (one object per input)
    #[clap(short, long, conflicts_with_all=["show_index", "show_headers", "num"])]
    pub json: bool,

    /// Print the block index (VBQ/CBQ only; BQ files have no index)
    #[clap(long, conflicts_with_all=["json", "show_headers", "num"])]
    pub show_index: bool,

    /// Print per-block headers in debug format (CBQ only; other inputs are skipped)
    #[clap(long, conflicts_with_all=["json", "show_index", "num"])]
    pub show_headers: bool,
}
