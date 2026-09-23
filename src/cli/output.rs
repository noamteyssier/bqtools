use anyhow::{bail, Result};
use binseq::{BitSize, Policy};
use clap::{
    builder::{PossibleValue, PossibleValuesParser, TypedValueParser},
    Parser, ValueEnum,
};
use std::{io::Write, path::Path};

use crate::{
    cli::FileFormat,
    commands::{compress_passthrough, match_output, CompressionType},
};

#[derive(Parser, Debug, Clone)]
#[clap(next_help_heading = "OUTPUT FILE OPTIONS")]
pub struct OutputFile {
    /// Output file [default: stdout]
    #[clap(short = 'o', long)]
    pub output: Option<String>,

    /// Write paired output to `{PREFIX}_R1.{ext}` and `{PREFIX}_R2.{ext}`
    ///
    /// Only valid for paired input with `--mate both`; a `.gz`/`.zst` suffix is
    /// added when compressing. Without it, paired records are interleaved.
    #[clap(short, long, conflicts_with = "output")]
    pub prefix: Option<String>,

    /// Which mate(s) to output for paired BINSEQ files
    ///
    /// Ignored (with a warning) for single-end files.
    #[clap(short = 'm', long, default_value = "both")]
    pub mate: Mate,

    /// Output file format [default: inferred from the output extension, else TSV]
    #[clap(short, long, value_parser = parse_record_format())]
    pub format: Option<FileFormat>,

    /// Compress output file [default: inferred from the output extension, else uncompressed]
    #[clap(short, long)]
    pub compress: Option<CompressionType>,

    /// Number of threads for processing and compression (0 = all CPUs; capped at CPU count)
    #[clap(short = 'T', long, default_value = "0")]
    pub threads: usize,
}
impl OutputFile {
    pub fn as_writer(&self) -> Result<Box<dyn Write + Send>> {
        let writer = match_output(self.output.as_deref())?;
        compress_passthrough(writer, self.compress(), self.threads())
    }

    /// Explicit `-c` wins; otherwise compression is inferred from the output extension.
    #[allow(clippy::case_sensitive_file_extension_comparisons)]
    pub fn compress(&self) -> CompressionType {
        if let Some(compress) = self.compress {
            return compress;
        }
        self.output
            .as_ref()
            .map_or(CompressionType::Uncompressed, |path| {
                if path.ends_with(".gz") {
                    CompressionType::Gzip
                } else if path.ends_with(".zst") {
                    CompressionType::Zstd
                } else {
                    CompressionType::Uncompressed
                }
            })
    }

    pub fn mate(&self) -> Mate {
        self.mate
    }

    pub fn format(&self) -> Result<FileFormat> {
        let format = if let Some(format) = self.format {
            format
        } else if let Some(path) = self.output.as_ref() {
            FileFormat::from_path(path).ok_or_else(|| {
                anyhow::anyhow!("Could not infer file format from `{path}`; pass -f a|q|t")
            })?
        } else {
            FileFormat::Tsv
        };

        // `-f` can't select BAM, but an output path ending in `.bam` can
        if format == FileFormat::Bam {
            bail!(
                "BAM output is not supported here; use FASTA (-f a), FASTQ (-f q), or TSV (-f t) instead"
            );
        }

        Ok(format)
    }

    /// Returns the number of threads to use.
    ///
    /// The default of 0 sets the maximum; all other values are clamped to the maximum.
    pub fn threads(&self) -> usize {
        match self.threads {
            0 => num_cpus::get(),
            n => n.min(num_cpus::get()),
        }
    }

    pub fn as_paired_writer(
        &self,
        format: FileFormat,
    ) -> Result<(Box<dyn Write + Send>, Box<dyn Write + Send>)> {
        // Check for prefix
        let prefix = self.prefix.as_ref().ok_or_else(|| {
            anyhow::anyhow!("Output file format prefix is required for paired BINSEQ files")
        })?;

        // Construct the output file names
        let compress = self.compress();
        let r1_name = if let Some(ext) = compress.extension() {
            format!("{}_R1.{}.{}", prefix, format.extension(), ext)
        } else {
            format!("{}_R1.{}", prefix, format.extension())
        };
        let r2_name = if let Some(ext) = compress.extension() {
            format!("{}_R2.{}.{}", prefix, format.extension(), ext)
        } else {
            format!("{}_R2.{}", prefix, format.extension())
        };

        // Open the output files
        let r1 = match_output(Some(&r1_name))?;
        let r2 = match_output(Some(&r2_name))?;

        // Compress the output files (if necessary)
        let r1 = compress_passthrough(r1, compress, self.threads())?;
        let r2 = compress_passthrough(r2, compress, self.threads())?;

        Ok((r1, r2))
    }
}

#[derive(ValueEnum, PartialEq, Eq, Clone, Copy, Debug, Default)]
pub enum Mate {
    /// Primary (R1) mate only
    #[clap(name = "1")]
    One,
    /// Extended (R2) mate only
    #[clap(name = "2")]
    Two,
    /// Both mates
    #[default]
    Both,
}

/// Record output formats (`-f`) for commands writing FASTA/FASTQ/TSV.
fn parse_record_format() -> impl TypedValueParser<Value = FileFormat> {
    PossibleValuesParser::new([
        PossibleValue::new("a").help("FASTA file format"),
        PossibleValue::new("q").help("FASTQ file format"),
        PossibleValue::new("t").help("TSV file format"),
    ])
    .map(|s| match s.as_str() {
        "a" => FileFormat::Fasta,
        "q" => FileFormat::Fastq,
        _ => FileFormat::Tsv,
    })
}

#[derive(Parser, Debug, Clone)]
#[clap(next_help_heading = "OUTPUT BINSEQ OPTIONS")]
#[allow(clippy::struct_excessive_bools)]
pub struct OutputBinseq {
    #[clap(short = 'o', long)]
    /// Output binseq file
    ///
    /// To output to stdout, use the `--pipe` flag.
    pub output: Option<String>,

    #[clap(flatten)]
    pub options: OutputBinseqOptions,

    /// Pipe the output to stdout
    #[clap(long, conflicts_with = "output")]
    pub pipe: bool,
}
impl OutputBinseq {
    pub fn mode(&self) -> Result<BinseqMode> {
        if let Some(mode) = self.options.mode {
            Ok(mode)
        } else if let Some(ref path) = self.output {
            BinseqMode::determine(path)
        } else {
            // STDOUT
            Ok(BinseqMode::default())
        }
    }

    pub fn threads(&self) -> usize {
        self.options.threads()
    }
}

/// BINSEQ output for commands that copy the input file's mode and encoding settings.
#[derive(Parser, Debug, Clone)]
#[clap(next_help_heading = "OUTPUT BINSEQ OPTIONS")]
pub struct OutputBinseqInherited {
    /// Output binseq file
    ///
    /// The output keeps the input's BINSEQ mode and settings, so a `.bq/.vbq/.cbq`
    /// extension must match the input. To output to stdout, use the `--pipe` flag.
    #[clap(short = 'o', long)]
    pub output: Option<String>,

    /// Pipe the output to stdout
    #[clap(long, conflicts_with = "output")]
    pub pipe: bool,

    /// Number of threads to use (0 = all CPUs; clamped to CPU count)
    #[clap(short = 'T', long, default_value = "0")]
    pub threads: usize,
}
impl OutputBinseqInherited {
    /// Opens the output, refusing a known BINSEQ extension that disagrees with `mode`.
    pub fn as_writer(&self, mode: BinseqMode) -> Result<Box<dyn Write + Send>> {
        if let Some(path) = self.output.as_deref() {
            if let Ok(ext_mode) = BinseqMode::determine(path) {
                if ext_mode != mode {
                    bail!(
                        "Output extension implies {ext_mode:?} but input is {mode:?}; the output always keeps the input's BINSEQ mode"
                    );
                }
            }
        }
        if self.output.is_none() && !self.pipe {
            bail!(
                "Refusing to write binary BINSEQ data to stdout. Provide an output path with `-o/--output`, or pass `--pipe` to write to stdout explicitly."
            );
        }
        match_output(self.output.as_deref())
    }

    pub fn threads(&self) -> usize {
        match self.threads {
            0 => num_cpus::get(),
            n => n.min(num_cpus::get()),
        }
    }
}

#[derive(Parser, Debug, Clone, Copy)]
#[allow(clippy::struct_excessive_bools)]
pub struct OutputBinseqOptions {
    /// BINSEQ mode to write [default: inferred from the -o extension, else cbq]
    #[clap(short = 'm', long)]
    pub mode: Option<BinseqMode>,

    /// Policy for handling Ns in sequences
    ///
    /// Only applied with 2-bit encoding (bq/vbq); cbq stores Ns directly.
    #[clap(short = 'p', long, default_value = "r")]
    pub policy: PolicyWrapper,

    /// Encoding bitsize (2 or 4 bits per nucleotide)
    ///
    /// Used by bq+vbq; ignored by cbq.
    #[clap(short = 'S', long, default_value = "2", value_parser = parse_bitsize())]
    bitsize: u8,

    /// Exclude sequence names (headers) in the binseq file
    ///
    /// Used by vbq+cbq
    #[clap(short = 'H', long)]
    skip_headers: bool,

    /// Skip ZSTD compression of VBQ blocks (default: compressed)
    ///
    /// Only used by vbq; bq is never compressed and cbq always is.
    #[clap(short = 'u', long)]
    pub uncompressed: bool,

    /// Skip inclusion of quality scores (default: included)
    ///
    /// Used by vbq+cbq
    #[clap(short = 'Q', long)]
    pub skip_quality: bool,

    /// Virtual block size in bytes; accepts K/M/G suffixes (powers of 1024)
    ///
    /// Used by vbq+cbq
    #[clap(short = 'B', long, value_parser = parse_memory_size, default_value = "128K")]
    block_size: usize,

    /// Number of threads to use (0 = all CPUs; capped at CPU count)
    ///
    /// When batch encoding several files, threads are split across concurrently
    /// encoded files.
    #[clap(short = 'T', long, default_value = "0")]
    pub threads: usize,

    /// Zstd compression level
    ///
    /// Between 1 and 22; higher levels compress better at the cost of speed.
    /// 0 uses zstd's default level.
    ///
    /// Only used by cbq (vbq always uses level 3).
    #[clap(short, long, default_value = "3")]
    pub level: i32,

    /// Archive mode
    ///
    /// Sets 4-bit encoding, keeps headers and quality scores, uses a 200M block
    /// size, and enables zstd compression. Intended for vbq; it does not change
    /// `--mode`, so pair it with `-o out.vbq` or `-m vbq`.
    #[clap(short = 'A', long, conflicts_with_all = ["uncompressed", "skip_headers", "bitsize", "block_size", "skip_quality", "level"])]
    pub archive: bool,
}
impl OutputBinseqOptions {
    pub fn headers(&self) -> bool {
        if self.archive {
            true
        } else {
            !self.skip_headers
        }
    }

    pub fn block_size(&self) -> usize {
        if self.archive {
            200 * 1024 * 1024
        } else {
            self.block_size
        }
    }

    pub fn compress(&self) -> bool {
        if self.archive {
            true
        } else {
            !self.uncompressed
        }
    }

    pub fn quality(&self) -> bool {
        if self.archive {
            true
        } else {
            !self.skip_quality
        }
    }

    pub fn threads(&self) -> usize {
        match self.threads {
            0 => num_cpus::get(),
            n => n.min(num_cpus::get()),
        }
    }

    pub fn bitsize(&self) -> BitSize {
        if self.archive {
            BitSize::Four
        } else {
            // `parse_bitsize` restricts the value to 2 or 4
            match self.bitsize {
                4 => BitSize::Four,
                _ => BitSize::Two,
            }
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum PolicyWrapper {
    /// Ignore any sequence if it contains an N
    #[clap(name = "i")]
    IgnoreSequence,

    /// Panic if any sequence contains an N
    #[clap(name = "p")]
    BreakOnInvalid,

    /// Randomly draw a nucleotide for each N in sequences.
    #[clap(name = "r")]
    RandomDraw,

    /// Sets all Ns to A
    #[clap(name = "a")]
    SetToA,

    /// Sets all Ns to C
    #[clap(name = "c")]
    SetToC,

    /// Sets all Ns to G
    #[clap(name = "g")]
    SetToG,

    /// Sets all Ns to T
    #[clap(name = "t")]
    SetToT,
}
impl From<PolicyWrapper> for Policy {
    fn from(value: PolicyWrapper) -> Self {
        match value {
            PolicyWrapper::IgnoreSequence => Policy::IgnoreSequence,
            PolicyWrapper::BreakOnInvalid => Policy::BreakOnInvalid,
            PolicyWrapper::RandomDraw => Policy::RandomDraw,
            PolicyWrapper::SetToA => Policy::SetToA,
            PolicyWrapper::SetToC => Policy::SetToC,
            PolicyWrapper::SetToG => Policy::SetToG,
            PolicyWrapper::SetToT => Policy::SetToT,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum, Default, PartialEq)]
pub enum BinseqMode {
    /// Fixed-length, 2/4-bit encoded records
    #[clap(name = "bq")]
    Bq,
    /// Variable-length, block-compressed records
    #[clap(name = "vbq")]
    Vbq,
    /// Columnar, block-compressed records (default)
    #[clap(name = "cbq")]
    #[default]
    Cbq,
}
#[cfg(test)]
impl BinseqMode {
    pub fn enum_iter() -> impl Iterator<Item = Self> + Clone {
        [Self::Bq, Self::Vbq, Self::Cbq].into_iter()
    }
}
impl BinseqMode {
    pub fn determine(path: &str) -> Result<Self> {
        let pathbuf = Path::new(path);
        if let Some(ext) = pathbuf.extension() {
            match ext.to_str() {
                Some("bq") => Ok(Self::Bq),
                Some("vbq") => Ok(Self::Vbq),
                Some("cbq") => Ok(Self::Cbq),
                _ => bail!("Could not determine BINSEQ output mode from path: {path}"),
            }
        } else {
            bail!("Could not determine BINSEQ output mode from path: {path}")
        }
    }
    pub fn extension(&self) -> &str {
        match self {
            Self::Bq => ".bq",
            Self::Vbq => ".vbq",
            Self::Cbq => ".cbq",
        }
    }
}
impl From<BinseqMode> for binseq::write::Format {
    fn from(val: BinseqMode) -> Self {
        match val {
            BinseqMode::Bq => binseq::write::Format::Bq,
            BinseqMode::Vbq => binseq::write::Format::Vbq,
            BinseqMode::Cbq => binseq::write::Format::Cbq,
        }
    }
}

fn parse_bitsize() -> impl TypedValueParser<Value = u8> {
    PossibleValuesParser::new(["2", "4"]).map(|s| s.parse::<u8>().unwrap())
}

fn parse_memory_size(input: &str) -> Result<usize, String> {
    let input = input.trim().to_uppercase();
    let last_char = input.chars().last().unwrap_or('0');

    let (number_str, multiplier) = match last_char {
        'K' | 'k' => (&input[..input.len() - 1], 1024),
        'M' | 'm' => (&input[..input.len() - 1], 1024 * 1024),
        'G' | 'g' => (&input[..input.len() - 1], 1024 * 1024 * 1024),
        _ if last_char.is_ascii_digit() => (input.as_str(), 1),
        _ => return Err(format!("Invalid memory size format: {input}")),
    };

    match number_str.parse::<usize>() {
        Ok(number) => Ok(number * multiplier),
        Err(_) => Err(format!("Failed to parse number: {number_str}")),
    }
}

#[derive(Copy, Clone)]
pub struct BinseqConfig {
    pub compress: bool,
    pub quality: bool,
    pub block_size: usize,
    pub policy: Policy,
    pub bitsize: BitSize,
    pub headers: bool,
    pub threads: usize,
    pub compression_level: i32,
}
impl From<OutputBinseqOptions> for BinseqConfig {
    fn from(options: OutputBinseqOptions) -> Self {
        BinseqConfig {
            compress: options.compress(),
            quality: options.quality(),
            block_size: options.block_size(),
            policy: options.policy.into(),
            bitsize: options.bitsize(),
            headers: options.headers(),
            threads: options.threads(),
            compression_level: options.level,
        }
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;

    use super::{BinseqMode, OutputBinseq, OutputBinseqInherited, OutputFile};
    use crate::commands::CompressionType;

    fn compress_for(flags: &[&str]) -> CompressionType {
        let mut argv = vec!["output"];
        argv.extend_from_slice(flags);
        OutputFile::try_parse_from(argv).unwrap().compress()
    }

    #[test]
    fn test_compress_explicit_flag_wins() {
        assert!(matches!(compress_for(&["-c", "g"]), CompressionType::Gzip));
        assert!(matches!(
            compress_for(&["-o", "x.fq", "-c", "z"]),
            CompressionType::Zstd
        ));
        assert!(matches!(
            compress_for(&["-o", "x.fq.gz", "-c", "u"]),
            CompressionType::Uncompressed
        ));
    }

    #[test]
    fn test_compress_inferred_from_extension() {
        assert!(matches!(
            compress_for(&["-o", "x.fq.gz"]),
            CompressionType::Gzip
        ));
        assert!(matches!(
            compress_for(&["-o", "x.fq.zst"]),
            CompressionType::Zstd
        ));
        assert!(matches!(
            compress_for(&["-o", "x.fq"]),
            CompressionType::Uncompressed
        ));
        assert!(matches!(compress_for(&[]), CompressionType::Uncompressed));
    }

    /// Without `-o` or `--pipe`, writing binary BINSEQ data to stdout must be
    /// refused rather than silently dumping binary into the terminal.
    #[test]
    fn test_as_writer_rejects_bare_stdout() {
        let args = OutputBinseqInherited::try_parse_from(["output"]).unwrap();
        assert!(args.as_writer(BinseqMode::Cbq).is_err());
    }

    #[test]
    fn test_pipe_conflicts_with_output() {
        assert!(OutputBinseq::try_parse_from(["output", "--pipe", "-o", "x.cbq"]).is_err());
    }

    #[test]
    fn test_bitsize_rejects_invalid_values() {
        assert!(OutputBinseq::try_parse_from(["output", "-S", "3"]).is_err());
        assert!(OutputBinseq::try_parse_from(["output", "-S", "4"]).is_ok());
    }

    #[test]
    fn test_as_writer_allows_explicit_pipe() {
        let args = OutputBinseqInherited::try_parse_from(["output", "--pipe"]).unwrap();
        assert!(args.as_writer(BinseqMode::Cbq).is_ok());
    }

    #[test]
    fn test_as_writer_allows_output_path() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let args =
            OutputBinseqInherited::try_parse_from(["output", "-o", tmp.path().to_str().unwrap()])
                .unwrap();
        assert!(args.as_writer(BinseqMode::Cbq).is_ok());
    }

    #[test]
    fn test_inherited_rejects_mismatched_extension() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.bq");
        let args = OutputBinseqInherited::try_parse_from(["output", "-o", path.to_str().unwrap()])
            .unwrap();
        assert!(args.as_writer(BinseqMode::Cbq).is_err());
        assert!(args.as_writer(BinseqMode::Bq).is_ok());
    }
}
