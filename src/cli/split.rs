use clap::Parser;

#[cfg(feature = "fuzzy")]
use super::FuzzyArgs;
use super::{InputBinseq, PatternFileArgs};

/// Split a BINSEQ file into one output per pattern alias
///
/// Each record is written to the file of the single alias it matches, named
/// `<basepath>/<alias>.<ext>` in the input's BINSEQ mode. The alias is the FASTA
/// header or TSV alias of the pattern, or the pattern itself; patterns sharing
/// an alias share a file. Records matching no alias, or more than one, go to
/// the unmatched file.
#[derive(Parser, Debug)]
#[clap(group = clap::ArgGroup::new("pattern_files").args(["file", "sfile", "xfile"]).required(true).multiple(true))]
pub struct SplitCommand {
    #[clap(flatten)]
    pub input: InputBinseq,

    #[clap(flatten)]
    pub split: SplitOptions,

    #[clap(flatten)]
    pub patterns: PatternFileArgs,

    #[cfg(feature = "fuzzy")]
    #[clap(flatten)]
    pub fuzzy_args: FuzzyArgs,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "SPLIT OPTIONS")]
#[allow(clippy::struct_excessive_bools)]
pub struct SplitOptions {
    /// Output directory for split files (created if missing)
    #[clap(long, default_value = "./split_outs")]
    pub basepath: String,

    /// Skip writing records that match no pattern or more than one alias
    #[clap(long)]
    pub skip_unmatched: bool,

    /// File stem of the unmatched output (`<basepath>/<name>.<ext>`)
    #[clap(long, default_value = "unmatched")]
    pub unmatched_basename: String,

    /// Remove output files with fewer than this many records.
    ///
    /// Defaults to 1, which removes empty output files. Set to 0 to keep all files.
    #[clap(long, default_value_t = 1)]
    pub min_records: usize,

    /// Denotes patterns are fixed strings (non-regex)
    ///
    /// Allows usage of Aho-Corasick algorithm for efficient matching.
    /// Auto-detected when all patterns are uppercase ACGT; forcing it makes
    /// regex-looking patterns match literally. Ignored with fuzzy matching.
    #[clap(short = 'x', long)]
    pub fixed: bool,

    /// Reverse complement all patterns before matching
    ///
    /// Only supported for fixed ACGT patterns; regex patterns are rejected
    /// since reverse complementing a regex is undefined.
    #[clap(long)]
    pub rc: bool,

    /// Don't use Aho-Corasick DFA (slower, but lower memory)
    ///
    /// Only affects fixed-string (Aho-Corasick) matching.
    #[clap(long)]
    pub no_dfa: bool,

    /// Number of threads to use (0 = all CPUs)
    #[clap(short = 'T', long, default_value_t = 0)]
    pub threads: usize,

    /// Suppress the `alias<TAB>count` summary written to stderr
    #[clap(long)]
    pub quiet: bool,
}
