use std::io::Write;

use anyhow::Result;

use super::Writer;
use crate::cli::{FileFormat, Mate};

pub fn write_fastq_parts<W: Write>(
    writer: &mut W,
    index: &[u8],
    sequence: &[u8],
    quality: &[u8],
) -> std::io::Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(index)?;
    writer.write_all(b"\n")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(quality)?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn write_fasta_parts<W: Write>(
    writer: &mut W,
    index: &[u8],
    sequence: &[u8],
) -> std::io::Result<()> {
    writer.write_all(b">")?;
    writer.write_all(index)?;
    writer.write_all(b"\n")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n")?;
    Ok(())
}

pub fn write_tsv_parts<W: Write>(
    writer: &mut W,
    index: &[u8],
    sequence: &[u8],
) -> std::io::Result<()> {
    writer.write_all(index)?;
    writer.write_all(b"\t")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n")?;
    Ok(())
}

/// Placeholder quality (`?`, Phred 30) of length `len`, growing `buf` as needed.
pub fn fill_qual(buf: &mut Vec<u8>, len: usize) -> &[u8] {
    if buf.len() < len {
        buf.resize(len, b'?');
    }
    &buf[..len]
}

pub enum SplitWriter {
    Interleaved { inner: Writer },
    Split { left: Writer, right: Writer },
}
impl SplitWriter {
    pub fn is_split(&self) -> bool {
        match self {
            Self::Interleaved { .. } => false,
            Self::Split { .. } => true,
        }
    }

    /// Write one batch: `left`/`right` when split, `mixed` when interleaved.
    pub fn write_batch(&mut self, left: &[u8], right: &[u8], mixed: &[u8]) -> std::io::Result<()> {
        match self {
            Self::Interleaved { inner } => inner.write_all(mixed),
            Self::Split { left: l, right: r } => {
                l.write_all(left)?;
                r.write_all(right)
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

pub fn write_record<W: Write>(
    writer: &mut W,
    header: &[u8],
    sequence: &[u8],
    quality: &[u8],
    format: FileFormat,
) -> Result<(), std::io::Error> {
    let qual_buf = &quality[..sequence.len()];
    match format {
        FileFormat::Fasta => write_fasta_parts(writer, header, sequence),
        FileFormat::Fastq => write_fastq_parts(writer, header, sequence, qual_buf),
        FileFormat::Tsv => write_tsv_parts(writer, header, sequence),
        FileFormat::Bam => unimplemented!("Cannot write BAM record from here"),
    }
}

#[allow(clippy::too_many_arguments)]
pub fn write_record_pair<W: Write>(
    left: &mut W,
    right: &mut W,
    mixed: &mut W,
    mate: Option<Mate>,
    split: bool,
    sbuf: &[u8],
    squal: &[u8],
    sheader: &[u8],
    xbuf: &[u8],
    xqual: &[u8],
    xheader: &[u8],
    format: FileFormat,
) -> Result<()> {
    match mate {
        Some(Mate::Both) => {
            if split {
                write_record(left, sheader, sbuf, squal, format)?;
                if !xbuf.is_empty() {
                    write_record(right, xheader, xbuf, xqual, format)?;
                }
            } else {
                write_record(mixed, sheader, sbuf, squal, format)?;
                if !xbuf.is_empty() {
                    write_record(mixed, xheader, xbuf, xqual, format)?;
                }
            }
        }
        Some(Mate::One) | None => {
            write_record(mixed, sheader, sbuf, squal, format)?;
        }
        Some(Mate::Two) => {
            write_record(mixed, xheader, xbuf, xqual, format)?;
        }
    }

    Ok(())
}
