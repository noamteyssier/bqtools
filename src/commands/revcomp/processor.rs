use std::{
    borrow::Cow,
    io::Write,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
};

use binseq::{BinseqRecord, BinseqWriter, ParallelProcessor, SequencingRecordBuilder};

use crate::cli::Mate;

/// Reverse complements a nucleotide sequence buffer in place.
///
/// Any byte outside `ACGTacgt` (e.g. `N`) is left untouched, matching the
/// behavior of 4-bit decoding, which collapses all ambiguity codes to `N`.
fn reverse_complement(buf: &mut [u8]) {
    buf.reverse();
    for base in buf.iter_mut() {
        let comp = match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => continue,
        };
        // keep the original case bit
        *base = comp | (*base & 0x20);
    }
}

/// Returns the sequence and quality, reverse complemented / reversed when `rc` is set.
fn orient<'a>(
    seq: &'a [u8],
    qual: Option<&'a [u8]>,
    rc: bool,
) -> (Cow<'a, [u8]>, Option<Cow<'a, [u8]>>) {
    if !rc {
        return (seq.into(), qual.map(Into::into));
    }
    let mut seq = seq.to_vec();
    reverse_complement(&mut seq);
    let qual = qual.map(|q| q.iter().rev().copied().collect());
    (seq.into(), qual)
}

pub struct RevCompProcessor<W: Write + Send> {
    /// Which mate(s) to reverse complement
    mate: Mate,

    /// Thread-local writer for the processor
    t_writer: BinseqWriter<Vec<u8>>,
    t_count: usize,

    /// Global writer for the processor
    writer: Arc<Mutex<BinseqWriter<W>>>,
    /// Global record count
    count: Arc<AtomicUsize>,
}
impl<W: Write + Send> Clone for RevCompProcessor<W> {
    fn clone(&self) -> Self {
        Self {
            mate: self.mate,
            t_writer: self.t_writer.clone(),
            t_count: 0,
            writer: self.writer.clone(),
            count: self.count.clone(),
        }
    }
}
impl<W: Write + Send> RevCompProcessor<W> {
    pub fn new(writer: BinseqWriter<W>, mate: Mate) -> binseq::Result<Self> {
        let t_writer = writer.new_headless_buffer()?;
        Ok(Self {
            mate,
            t_writer,
            t_count: 0,
            writer: Arc::new(Mutex::new(writer)),
            count: Arc::new(AtomicUsize::new(0)),
        })
    }

    pub fn finish(&mut self) -> binseq::Result<()> {
        self.writer.lock().unwrap().finish()
    }

    pub fn get_global_record_count(&self) -> usize {
        self.count.load(Ordering::Relaxed)
    }
}

impl<W: Write + Send> ParallelProcessor for RevCompProcessor<W> {
    fn process_record<B: BinseqRecord>(&mut self, record: B) -> binseq::Result<()> {
        let is_paired = record.is_paired();
        let has_quality = record.has_quality();
        let rc_primary = matches!(self.mate, Mate::One | Mate::Both);
        let rc_extended = is_paired && matches!(self.mate, Mate::Two | Mate::Both);

        let (s_seq, s_qual) = orient(
            record.sseq(),
            has_quality.then(|| record.squal()),
            rc_primary,
        );
        // an empty `x_seq` would mark the record as paired, so only set mate fields when paired
        let (x_seq, x_qual) = orient(
            record.xseq(),
            (is_paired && has_quality).then(|| record.xqual()),
            rc_extended,
        );

        let rec = SequencingRecordBuilder::default()
            .s_seq(&s_seq)
            .opt_s_qual(s_qual.as_deref())
            .s_header(record.sheader())
            .opt_x_seq(is_paired.then_some(&*x_seq))
            .opt_x_qual(x_qual.as_deref())
            .opt_x_header(is_paired.then(|| record.xheader()))
            .build()?;

        if self.t_writer.push(rec)? {
            self.t_count += 1;
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        self.count.fetch_add(self.t_count, Ordering::Relaxed);
        self.t_count = 0;
        self.writer
            .lock()
            .unwrap()
            .ingest_completed(&mut self.t_writer)
    }

    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.writer.lock().unwrap().ingest(&mut self.t_writer)
    }
}

#[cfg(test)]
mod tests {
    use super::reverse_complement;

    #[test]
    fn test_reverse_complement_basic() {
        let mut seq = b"ACGTACGT".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(seq, b"ACGTACGT");

        let mut seq = b"GATTACA".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(seq, b"TGTAATC");
    }

    #[test]
    fn test_reverse_complement_preserves_n() {
        let mut seq = b"ACGTN".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(seq, b"NACGT");
    }
}
