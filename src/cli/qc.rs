use clap::Parser;

use super::InputBinseq;

/// Run FastQC-style quality control on a BINSEQ file.
///
/// Writes `summary.md` and per-module TSVs to the output directory. At least
/// one module must remain enabled.
#[derive(Parser, Debug)]
pub struct QcCommand {
    #[clap(flatten)]
    pub input: InputBinseq,

    #[clap(flatten)]
    pub qc: QcOptions,
}

// Each `skip_*` flag independently toggles one QC module on/off - they're
// orthogonal CLI switches, not states in a state machine, so collapsing them
// into an enum wouldn't fit clap's flag model or improve readability here.
#[allow(clippy::struct_excessive_bools)]
#[derive(Parser, Debug)]
#[clap(next_help_heading = "QC OPTIONS")]
pub struct QcOptions {
    /// Number of threads to use (0 = all CPUs)
    #[clap(short = 'T', long, default_value_t = 0)]
    pub threads: usize,

    /// Skip the per-base sequence quality module
    #[clap(long)]
    pub skip_base_qual: bool,

    /// Skip the per-sequence quality module
    #[clap(long)]
    pub skip_seq_qual: bool,

    /// Skip the per-base sequence content module
    #[clap(long)]
    pub skip_base_content: bool,

    /// Skip the per-sequence GC content module
    #[clap(long)]
    pub skip_seq_gc: bool,

    /// Skip the sequence length distribution module
    #[clap(long)]
    pub skip_seq_length: bool,

    /// Skip the sequence duplication levels module
    #[clap(long)]
    pub skip_dup_levels: bool,

    /// Skip the overrepresented sequences module
    #[clap(long)]
    pub skip_overrepresented: bool,

    /// Number of leading records (of the processed span) to sample for
    /// duplication level and overrepresented-sequence estimation
    ///
    /// 0 uses all records; memory then grows with the number of distinct sequences.
    #[clap(long, default_value_t = 100_000)]
    pub dup_sample_size: usize,

    /// Minimum percent (0-100) of sampled reads a sequence must represent to be
    /// flagged as overrepresented (0.1 = 0.1%)
    #[clap(long, default_value_t = 0.1, value_parser = parse_percent)]
    pub overrepresented_threshold: f64,

    /// Output directory for the report and TSVs (created if missing)
    #[clap(short, long, default_value = "./bqtools-qc")]
    pub outdir: String,
}

fn parse_percent(input: &str) -> Result<f64, String> {
    let value: f64 = input
        .parse()
        .map_err(|_| format!("Invalid percentage: {input}"))?;
    if !(0.0..=100.0).contains(&value) {
        return Err(format!(
            "expected a percentage between 0 and 100, got {value}"
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::parse_percent;

    #[test]
    fn test_parse_percent_bounds() {
        assert_eq!(parse_percent("0.1"), Ok(0.1));
        assert_eq!(parse_percent("100"), Ok(100.0));
        assert!(parse_percent("-1").is_err());
        assert!(parse_percent("101").is_err());
        assert!(parse_percent("x").is_err());
    }
}
