use std::{path::PathBuf, str::FromStr};

use anyhow::{bail, Result};
use binseq::BinseqReader;
use clap::{
    builder::{PossibleValue, PossibleValuesParser, TypedValueParser},
    Parser,
};
use log::{debug, error, warn};
use paraseq::{fastx, ReaderBuilder};

use crate::{cli::BinseqMode, types::BoxedReader};

use super::FileFormat;

#[derive(Parser, Debug, Clone)]
#[clap(next_help_heading = "INPUT FILE OPTIONS")]
pub struct InputFile {
    /// Input file(s) [default: stdin]
    ///
    /// Zero (stdin), one, or two (an R1/R2 pair) files are encoded to a single
    /// output. More than two files are each encoded separately (use `--collate`
    /// to merge them); these must match `*.{fastq,fq,fasta,fa}[.gz|.zst]`
    /// (and `_R1`/`_R2` with `--paired`), others are skipped.
    #[clap(num_args = 0..)]
    pub input: Vec<String>,

    /// Input file format
    ///
    /// FASTA/FASTQ are auto-detected. Use `b` for SAM/BAM/CRAM from stdin or with
    /// a non-standard extension (a single `.sam/.bam/.cram` path is detected).
    #[clap(short, long, value_parser = parse_input_format())]
    format: Option<FileFormat>,

    /// Batch size (in records) to use in parallel processing
    ///
    /// Set this to a lower value for embedding genomes to better
    /// make use of parallelism (e.g. 2-4). Not applied to SAM/BAM/CRAM input.
    #[clap(short, long)]
    pub batch_size: Option<usize>,

    /// Input is paired-interleaved (alternating R1/R2 records, or name-sorted paired SAM/BAM/CRAM)
    #[clap(short = 'I', long, conflicts_with = "paired")]
    pub interleaved: bool,

    /// Encode all fasta/fastq files found under a single input directory
    ///
    /// Each file (or R1/R2 pair with `--paired`) is written next to its input
    /// unless `--collate` is set. Options used will be applied to all files.
    #[clap(short = 'r', long)]
    pub recursive: bool,

    /// Path to a text file listing input files to process (one per line)
    ///
    /// For R1/R2 encodings pair this with the `--paired` option. Paths go through
    /// the same extension filter as positional inputs. Options used will be
    /// applied to all files in the manifest.
    #[clap(short = 'M', long, conflicts_with_all = ["recursive", "input"])]
    pub manifest: Option<String>,

    #[clap(flatten)]
    pub recursion: RecursiveOptions,

    #[clap(flatten)]
    pub batch_encoding_options: BatchEncodingOptions,
}
impl InputFile {
    pub fn single_path(&self) -> Result<Option<&str>> {
        match self.input.len() {
            0 => Ok(None),
            1 => Ok(Some(&self.input[0])),
            _ => bail!("Requested single input file, but multiple files were provided."),
        }
    }

    pub fn paired_paths(&self) -> Result<(&str, &str)> {
        match self.input.len() {
            2 => Ok((&self.input[0], &self.input[1])),
            _ => bail!("Two input files are required."),
        }
    }

    /// Two inputs are implicitly treated as a pair unless they are collated or interleaved.
    pub fn paired(&self) -> bool {
        self.batch_encoding_options.paired
            || (self.input.len() == 2 && !self.batch_encoding_options.collate && !self.interleaved)
    }

    /// Returns the number of input files.
    pub fn num_files(&self) -> usize {
        self.input.len()
    }

    pub fn format(&self) -> Option<FileFormat> {
        if let Some(format) = self.format {
            Some(format)
        } else if self.input.len() == 1 {
            let path = &self.input[0];
            let p = std::path::Path::new(path);
            if p.extension().is_some_and(|ext| {
                ext.eq_ignore_ascii_case("bam")
                    || ext.eq_ignore_ascii_case("sam")
                    || ext.eq_ignore_ascii_case("cram")
            }) {
                Some(FileFormat::Bam)
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn is_stdin(&self) -> bool {
        self.input.is_empty()
    }

    pub fn as_directory(&self) -> Result<PathBuf> {
        if !self.recursive {
            bail!("Recursive mode is required to process a directory.");
        }
        let path = match self.input.as_slice() {
            [dir] => PathBuf::from(dir),
            [] => bail!("Recursive mode requires a directory as input."),
            _ => bail!("Recursive mode accepts exactly one directory as input."),
        };
        if !path.is_dir() {
            bail!("Input path is not a directory: {}", path.display());
        }
        Ok(path)
    }

    pub fn build_single_reader(&self) -> Result<fastx::Reader<BoxedReader>> {
        let path = self.single_path()?;
        let reader = load_reader(path, self.batch_size)?;
        Ok(reader)
    }

    /// Builds a vector of readers from the input paths.
    fn build_readers_from_paths(&self) -> Result<Vec<fastx::Reader<BoxedReader>>> {
        self.input
            .iter()
            .map(|path| load_reader(Some(path), self.batch_size))
            .collect()
    }

    pub fn build_single_collection(&self) -> Result<fastx::Collection<BoxedReader>> {
        self.build_collection_with_optional_stdin(fastx::CollectionType::Single)
    }

    pub fn build_interleaved_collection(&self) -> Result<fastx::Collection<BoxedReader>> {
        self.build_collection_with_optional_stdin(fastx::CollectionType::Interleaved)
    }

    pub fn build_paired_collection(&self) -> Result<fastx::Collection<BoxedReader>> {
        if self.input.is_empty() {
            bail!("Cannot build paired collection from stdin");
        }
        if !self.input.len().is_multiple_of(2) {
            bail!("Input must contain an even number of paths for paired collection");
        }
        let collection = fastx::Collection::new(
            self.build_readers_from_paths()?,
            fastx::CollectionType::Paired,
        )?;
        Ok(collection)
    }

    fn build_collection_with_optional_stdin(
        &self,
        collection_type: fastx::CollectionType,
    ) -> Result<fastx::Collection<BoxedReader>> {
        let collection = if self.input.is_empty() {
            fastx::Collection::new(vec![self.build_single_reader()?], collection_type)
        } else {
            fastx::Collection::new(self.build_readers_from_paths()?, collection_type)
        }?;
        Ok(collection)
    }
}

/// Encode input formats (`-f`); TSV is output-only.
fn parse_input_format() -> impl TypedValueParser<Value = FileFormat> {
    PossibleValuesParser::new([
        PossibleValue::new("a").help("FASTA file format"),
        PossibleValue::new("q").help("FASTQ file format"),
        PossibleValue::new("b").help("SAM/BAM/CRAM file format"),
    ])
    .map(|s| match s.as_str() {
        "a" => FileFormat::Fasta,
        "q" => FileFormat::Fastq,
        _ => FileFormat::Bam,
    })
}

fn load_reader(
    path: Option<&str>,
    batch_size: Option<usize>,
) -> Result<fastx::Reader<BoxedReader>> {
    if let Some(path) = path {
        if path.starts_with("gs://") {
            #[cfg(not(feature = "gcs"))]
            {
                error!("Missing feature flag - gcs. To process Google Cloud Storage files, enable the 'gcs' feature flag.");
                bail!("Missing feature flag - gcs");
            }

            #[cfg(feature = "gcs")]
            return Ok(load_gcs_reader(path, batch_size)?);
        }
        Ok(load_simple_reader(Some(path), batch_size)?)
    } else {
        Ok(load_simple_reader(None, batch_size)?)
    }
}

fn load_simple_reader(
    path: Option<&str>,
    batch_size: Option<usize>,
) -> Result<fastx::Reader<BoxedReader>, paraseq::Error> {
    let path_display = if let Some(path) = path {
        path.to_string()
    } else {
        "stdin".to_string()
    };

    debug!("building on-disk fastx reader (batch size: {batch_size:?}) from: {path_display}");
    let mut builder = ReaderBuilder::optional_path(path);
    if let Some(size) = batch_size {
        builder = builder.batch_size(size);
    }
    builder.build()
}

#[cfg(feature = "gcs")]
fn load_gcs_reader(
    path: &str,
    batch_size: Option<usize>,
) -> Result<fastx::Reader<BoxedReader>, paraseq::Error> {
    debug!("building GCS fastx reader (batch size: {batch_size:?}) from: {path}");
    let mut builder = ReaderBuilder::gcs(path);
    if let Some(size) = batch_size {
        builder = builder.batch_size(size);
    }
    builder.build()
}

#[derive(Parser, Debug, Clone, PartialEq, Eq)]
#[clap(next_help_heading = "RECURSION OPTIONS")]
pub struct RecursiveOptions {
    /// Maximum directory depth to search (1 = only the directory's own files; default: no limit)
    #[clap(long, requires = "recursive")]
    pub depth: Option<usize>,
}

#[derive(Parser, Debug, Clone, PartialEq, Eq)]
#[clap(next_help_heading = "BATCH ENCODING OPTIONS")]
pub struct BatchEncodingOptions {
    /// Treat inputs as *{_R1,_R2}* record pairs (the number of files must be even)
    ///
    /// Two positional files are paired by default unless `--collate` or `-I` is set.
    #[clap(short = 'P', long)]
    pub paired: bool,

    /// Collate all input files into a single output file. Will respect paired records if `--paired` is specified.
    #[clap(short = 'C', long)]
    pub collate: bool,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "INPUT FILE OPTIONS")]
pub struct InputBinseq {
    /// Input BINSEQ file (.bq/.vbq/.cbq)
    pub input: String,

    /// Span of records to process, as `START..END` [default: all records]
    ///
    /// 0-based and end-exclusive; either bound may be omitted (e.g. `100..`,
    /// `..500`). An end past the last record is clipped.
    #[clap(long)]
    pub span: Option<Span>,
}
impl InputBinseq {
    pub fn path(&self) -> &str {
        &self.input
    }

    pub fn mode(&self) -> Result<BinseqMode> {
        let reader = BinseqReader::new(&self.input)?;
        match reader {
            BinseqReader::Bq(_) => Ok(BinseqMode::Bq),
            BinseqReader::Vbq(_) => Ok(BinseqMode::Vbq),
            BinseqReader::Cbq(_) => Ok(BinseqMode::Cbq),
        }
    }
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "INPUT FILE OPTIONS")]
pub struct MultiInputBinseq {
    /// Input BINSEQ files; all must share the same mode and file header
    #[clap(num_args = 1..)]
    pub input: Vec<String>,
}

#[derive(Debug, Clone, Copy)]
pub struct Span {
    start: Option<usize>,
    end: Option<usize>,
}
impl Span {
    fn validate(&mut self, max_records: usize) -> Result<()> {
        if let Some(start) = self.start {
            if start > max_records {
                error!(
                    "Provided start ({start}) exceeds maximum number of records ({max_records})"
                );
                bail!("Maximum number of records exceeded")
            }
        }
        if let Some(end) = self.end {
            if end > max_records {
                warn!(
                    "Clipping provided endpoint ({end}) to maximum number of records ({max_records})"
                );
            }
            self.end = Some(end.min(max_records));
        }
        Ok(())
    }
    pub fn get_range(&mut self, max_records: usize) -> Result<std::ops::Range<usize>> {
        self.validate(max_records)?;
        match (self.start, self.end) {
            (Some(start), Some(end)) => Ok(start..end),
            (Some(start), None) => Ok(start..max_records),
            (None, Some(end)) => Ok(0..end),
            (None, None) => Ok(0..max_records),
        }
    }
}

impl FromStr for Span {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (start_str, end_str) = s
            .split_once("..")
            .ok_or_else(|| format!("expected range like '10..20', got '{s}'"))?;

        let parse_bound = |bound_str: &str, name: &str| {
            if bound_str.is_empty() {
                Ok(None)
            } else {
                bound_str
                    .parse()
                    .map(Some)
                    .map_err(|_| format!("invalid {name}: '{bound_str}'"))
            }
        };

        let start = parse_bound(start_str, "start")?;
        let end = parse_bound(end_str, "end")?;

        Ok(Self { start, end })
    }
}
