use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use aho_corasick::{AhoCorasick, AhoCorasickBuilder, AhoCorasickKind};
#[cfg(feature = "fuzzy")]
use anyhow::bail;
use anyhow::Result;
use binseq::{BinseqReader, BinseqWriterBuilder};
use gzp::{
    deflate::Gzip,
    par::compress::{ParCompress, ParCompressBuilder},
};
use log::trace;
#[cfg(feature = "fuzzy")]
use sassy::{profiles::Iupac, EncodedPatterns, Searcher};

/// Builds a case-sensitive Aho-Corasick automaton (DFA unless `no_dfa`).
pub fn corasick_builder(patterns: &[Vec<u8>], no_dfa: bool) -> Result<AhoCorasick> {
    Ok(AhoCorasickBuilder::new()
        .ascii_case_insensitive(false)
        .kind(if no_dfa {
            None
        } else {
            Some(AhoCorasickKind::DFA)
        })
        .build(patterns)?)
}

pub fn match_output<P: AsRef<Path>>(path: Option<P>) -> Result<Box<dyn Write + Send>> {
    let inner: Box<dyn Write + Send> = if let Some(path) = path {
        trace!("Opening writer handle at: {}", path.as_ref().display());
        Box::new(File::create(path)?)
    } else {
        trace!("Opening writer handle to stdout");
        Box::new(io::stdout())
    };
    Ok(Box::new(BufWriter::new(inner)))
}

/// Builds a writer that mirrors the input file's own header/configuration.
pub fn builder_from_reader(reader: &BinseqReader) -> BinseqWriterBuilder {
    match reader {
        BinseqReader::Bq(r) => BinseqWriterBuilder::from_bq_header(r.header()),
        BinseqReader::Vbq(r) => BinseqWriterBuilder::from_vbq_header(r.header()),
        BinseqReader::Cbq(r) => BinseqWriterBuilder::from_cbq_header(r.header()),
    }
}

#[derive(Clone, Copy, Default, Debug, clap::ValueEnum)]
pub enum CompressionType {
    /// Uncompressed
    #[default]
    #[value(name = "u")]
    Uncompressed,
    /// Gzip
    #[value(name = "g")]
    Gzip,
    /// Zstd
    #[value(name = "z")]
    Zstd,
}
impl CompressionType {
    pub fn extension(self) -> Option<&'static str> {
        match self {
            CompressionType::Uncompressed => None,
            CompressionType::Gzip => Some("gz"),
            CompressionType::Zstd => Some("zst"),
        }
    }
}

pub fn compress_passthrough(
    writer: Box<dyn Write + Send>,
    compression_type: CompressionType,
    num_threads: usize,
) -> Result<Box<dyn Write + Send>> {
    match compression_type {
        CompressionType::Uncompressed => Ok(writer),
        CompressionType::Gzip => {
            let encoder: ParCompress<Gzip, _> = ParCompressBuilder::new()
                .num_threads(num_threads)?
                .from_writer(writer);
            Ok(Box::new(encoder))
        }
        CompressionType::Zstd => {
            let mut encoder = zstd::Encoder::new(writer, 3)?;
            encoder.multithread(num_threads as u32)?;
            Ok(Box::new(encoder.auto_finish()))
        }
    }
}

/// Builds a fuzzy (sassy) searcher for one pattern set: validates uniform
/// pattern length, resolves the `max_n_frac` filter (defaulting to
/// `k / pattern_length` when unset), and encodes the patterns into the
/// searcher if any were provided.
#[cfg(feature = "fuzzy")]
pub fn build_fuzzy_searcher(
    patterns: &[Vec<u8>],
    k: usize,
    max_n_frac: Option<f32>,
) -> Result<(Searcher<Iupac>, Option<EncodedPatterns<Iupac>>)> {
    // sassy panics on mixed pattern lengths, so reject them up front
    let pattern_len = patterns.first().map_or(0, Vec::len);
    if let Some(bad) = patterns.iter().find(|p| p.len() != pattern_len) {
        log::error!("Multiple pattern lengths provided - currently cannot handle variable-length patterns in fuzzy matching");
        bail!(
            "Pattern length mismatch: expected length {pattern_len}, found length {}",
            bad.len()
        );
    }
    // default: k / pattern_length, or 1.0 (no restriction) without patterns
    let frac = max_n_frac.unwrap_or_else(|| {
        if pattern_len == 0 {
            1.0
        } else {
            (k as f32 / pattern_len as f32).min(1.0)
        }
    });
    let mut searcher = Searcher::new_fwd().with_max_n_frac(frac);
    let encoded = (!patterns.is_empty()).then(|| searcher.encode_patterns(patterns));
    Ok((searcher, encoded))
}
