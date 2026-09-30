use std::io::Write;

use anyhow::Result;

use crate::cli::{FileFormat, Mate};
use crate::types::BoxedWriter;

/// Placeholder quality (`?`, Phred 30) of length `len`, growing `buf` as needed.
pub fn fill_qual(buf: &mut Vec<u8>, len: usize) -> &[u8] {
    if buf.len() < len {
        buf.resize(len, b'?');
    }
    &buf[..len]
}

pub enum SplitWriter {
    Interleaved {
        inner: BoxedWriter,
    },
    Split {
        left: BoxedWriter,
        right: BoxedWriter,
    },
}
impl SplitWriter {
    pub fn is_split(&self) -> bool {
        match self {
            Self::Interleaved { .. } => false,
            Self::Split { .. } => true,
        }
    }

    /// Write one batch: `left`/`right` when split, `mixed` when interleaved.
    pub fn write_batch(&mut self, batch: &Batch) -> std::io::Result<()> {
        match self {
            Self::Interleaved { inner } => inner.write_all(&batch.mixed),
            Self::Split { left: l, right: r } => {
                l.write_all(&batch.left)?;
                r.write_all(&batch.right)
            }
        }?;
        self.flush()
    }

    pub fn flush(&mut self) -> Result<(), std::io::Error> {
        match self {
            SplitWriter::Interleaved { inner } => {
                inner.flush()?;
                Ok(())
            }
            SplitWriter::Split { left, right } => {
                left.flush()?;
                right.flush()?;
                Ok(())
            }
        }
    }
}

// ANSI red + bold, and reset
const RED_BOLD: &[u8] = b"\x1b[31;1m";
const RESET: &[u8] = b"\x1b[0m";

/// Half-open `(start, end)` span of a hit within a sequence.
type Interval = (usize, usize);

/// Sorts `matches` and merges overlapping or touching intervals in place.
fn merge_matches(matches: &mut Vec<Interval>) {
    matches.sort_unstable();
    let mut write = 0;
    for read in 1..matches.len() {
        if matches[read].0 <= matches[write].1 {
            matches[write].1 = matches[write].1.max(matches[read].1);
        } else {
            write += 1;
            matches[write] = matches[read];
        }
    }
    matches.truncate((write + 1).min(matches.len()));
}

/// Writes `buffer`, wrapping each span of `hits` (if any) in color codes.
fn write_marked<W: Write>(
    writer: &mut W,
    buffer: &[u8],
    hits: Option<&mut Vec<Interval>>,
) -> std::io::Result<()> {
    let Some(hits) = hits else {
        return writer.write_all(buffer);
    };
    merge_matches(hits);
    let mut pos = 0;
    for &(start, end) in hits.iter() {
        writer.write_all(&buffer[pos..start])?;
        writer.write_all(RED_BOLD)?;
        writer.write_all(&buffer[start..end])?;
        writer.write_all(RESET)?;
        pos = end;
    }
    writer.write_all(&buffer[pos..])
}

/// One read of a record: its header, sequence and (at least as long) quality.
#[derive(Clone, Copy)]
pub struct SeqRead<'a> {
    pub header: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

/// Writes `read`, coloring the `hits` in its sequence and quality if given.
pub fn write_record<W: Write>(
    writer: &mut W,
    read: SeqRead,
    hits: Option<&mut Vec<Interval>>,
    format: FileFormat,
) -> Result<(), std::io::Error> {
    let (head, mid): (&[u8], &[u8]) = match format {
        FileFormat::Fasta => (b">", b"\n"),
        FileFormat::Fastq => (b"@", b"\n"),
        FileFormat::Tsv => (b"", b"\t"),
        FileFormat::Bam => unimplemented!("Cannot write BAM record from here"),
    };
    writer.write_all(head)?;
    writer.write_all(read.header)?;
    writer.write_all(mid)?;
    if format == FileFormat::Fastq {
        // the same hits color both lines
        let mut hits = hits;
        write_marked(writer, read.seq, hits.as_deref_mut())?;
        writer.write_all(b"\n+\n")?;
        write_marked(writer, &read.qual[..read.seq.len()], hits)?;
    } else {
        write_marked(writer, read.seq, hits)?;
    }
    writer.write_all(b"\n")
}

/// Where the hits of each mate are, to color them.
pub type PairHits<'a> = Option<[&'a mut Vec<Interval>; 2]>;

/// A thread-local buffer of records waiting to be written by a [`SplitWriter`].
#[derive(Clone, Default)]
pub struct Batch {
    /// Interleaved records, or singlets
    mixed: Vec<u8>,
    /// R1 / R2 when writing to a pair of files
    left: Vec<u8>,
    right: Vec<u8>,
    split: bool,
}
impl Batch {
    pub fn new(writer: &SplitWriter) -> Self {
        Self {
            split: writer.is_split(),
            ..Self::default()
        }
    }

    pub fn clear(&mut self) {
        self.mixed.clear();
        self.left.clear();
        self.right.clear();
    }

    /// Appends the mate(s) selected by `mate`, coloring `hits` if given.
    pub fn push_pair(
        &mut self,
        mate: Mate,
        r1: SeqRead,
        r2: SeqRead,
        hits: PairHits,
        format: FileFormat,
    ) -> Result<()> {
        let [h1, h2] = match hits {
            Some([a, b]) => [Some(a), Some(b)],
            None => [None, None],
        };
        match mate {
            Mate::Both => {
                let first = if self.split {
                    &mut self.left
                } else {
                    &mut self.mixed
                };
                write_record(first, r1, h1, format)?;
                if !r2.seq.is_empty() {
                    let second = if self.split {
                        &mut self.right
                    } else {
                        &mut self.mixed
                    };
                    write_record(second, r2, h2, format)?;
                }
            }
            Mate::One => write_record(&mut self.mixed, r1, h1, format)?,
            Mate::Two => write_record(&mut self.mixed, r2, h2, format)?,
        }
        Ok(())
    }
}
