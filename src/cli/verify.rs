use clap::Parser;

use super::{InputBinseq, Mate};

/// Compute an order-independent checksum over a BINSEQ file.
///
/// BINSEQ files are frequently produced by parallel encoders, which make no
/// guarantee that record order matches the input FASTQ/FASTA. `verify`
/// accounts for this by hashing each record independently and combining the
/// per-record hashes with a commutative operation (wrapping sum), so the
/// resulting checksum is identical regardless of record order.
///
/// Use this to confirm that two BINSEQ files carry the same data even if a
/// parallel encoder wrote them in different record orders. Caveats: bq/vbq
/// encodes of input containing `N` differ between runs under the default
/// random N policy (use `-p a`), headers are never stored in bq files, and
/// `--span` selects records by file position, so spans are order-dependent.
///
/// Prints `<16-hex checksum>\t<num_records>\t<path>`.
#[derive(Parser, Debug)]
pub struct VerifyCommand {
    #[clap(flatten)]
    pub input: InputBinseq,

    #[clap(flatten)]
    pub opts: VerifyOptions,
}

// Each `skip_*` flag independently toggles one field out of the checksum -
// they're orthogonal CLI switches, not states in a state machine, so
// collapsing them into an enum wouldn't fit clap's flag model here.
#[allow(clippy::struct_excessive_bools)]
#[derive(Parser, Debug)]
#[clap(next_help_heading = "VERIFY OPTIONS")]
pub struct VerifyOptions {
    /// Exclude sequence data from the checksum
    #[clap(long)]
    pub skip_seq: bool,

    /// Exclude quality scores from the checksum (no effect on files without qualities)
    #[clap(long)]
    pub skip_qual: bool,

    /// Exclude sequence/record headers from the checksum
    ///
    /// Headers are automatically excluded (with a warning) for files that store
    /// none, such as all bq files.
    #[clap(long)]
    pub skip_headers: bool,

    /// Exclude the per-record flag from the checksum (no effect on files without flags)
    ///
    /// At least one of sequence, quality, headers, or flags must remain included.
    #[clap(long)]
    pub skip_flags: bool,

    /// Which mate(s) to include in the checksum for paired records
    ///
    /// `1` and `both` work on single-end files (both resolve to the primary
    /// channel, which always exists). `2` errors on single-end files, since
    /// there is no extended/mate-2 channel to checksum.
    // `-M` rather than `-m` for consistency with `revcomp --mate`.
    #[clap(short = 'M', long, default_value = "both")]
    pub mate: Mate,

    /// Number of threads to use (0 = all CPUs)
    #[clap(short = 'T', long, default_value_t = 0)]
    pub threads: usize,

    /// Print the checksum report as JSON
    #[clap(short, long)]
    pub json: bool,
}
