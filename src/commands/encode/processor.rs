use std::{
    io::Write,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc,
    },
};

use binseq::{BinseqWriter, SequencingRecord, SequencingRecordBuilder};
use paraseq::{
    prelude::{PairedParallelProcessor, ParallelProcessor},
    IntoParaseqError,
};
use std::sync::Mutex;

pub struct Encoder<W: Write + Send> {
    /// Thread-local writer for the encoder.
    t_writer: BinseqWriter<Vec<u8>>,
    /// Thread-local record count for the encoder.
    t_count: usize,
    /// Thread-local skip count for the encoder.
    t_skip: usize,

    /// Global writer for the encoder.
    writer: Arc<Mutex<BinseqWriter<W>>>,
    /// Global record count for the encoder.
    count: Arc<AtomicUsize>,
    /// Global skip count for the encoder.
    skip: Arc<AtomicUsize>,
}
impl<W: Write + Send> Clone for Encoder<W> {
    fn clone(&self) -> Self {
        Self {
            t_writer: self.t_writer.clone(),
            t_count: self.t_count,
            t_skip: self.t_skip,
            writer: self.writer.clone(),
            count: self.count.clone(),
            skip: self.skip.clone(),
        }
    }
}
impl<W: Write + Send> Encoder<W> {
    pub fn new(writer: BinseqWriter<W>) -> binseq::Result<Self> {
        let t_writer = writer.new_headless_buffer()?;
        Ok(Self {
            writer: Arc::new(Mutex::new(writer)),
            t_writer,
            t_count: 0,
            t_skip: 0,
            count: Arc::new(AtomicUsize::new(0)),
            skip: Arc::new(AtomicUsize::new(0)),
        })
    }

    /// Push a record into the thread-local writer, counting it as written or skipped.
    fn push(&mut self, rec: SequencingRecord<'_>) -> binseq::Result<()> {
        if self.t_writer.push(rec)? {
            self.t_count += 1;
        } else {
            self.t_skip += 1;
        }
        Ok(())
    }

    fn batch_complete(&mut self) -> binseq::Result<()> {
        self.count.fetch_add(self.t_count, Ordering::Relaxed);
        self.skip.fetch_add(self.t_skip, Ordering::Relaxed);
        self.t_count = 0;
        self.t_skip = 0;
        self.writer
            .lock()
            .unwrap()
            .ingest_completed(&mut self.t_writer)
    }

    fn thread_complete(&mut self) -> binseq::Result<()> {
        self.writer.lock().unwrap().ingest(&mut self.t_writer)
    }

    pub fn finish(&mut self) -> binseq::Result<()> {
        self.writer.lock().unwrap().finish()
    }

    /// Global `(written, skipped)` record counts.
    pub fn counts(&self) -> (usize, usize) {
        (
            self.count.load(Ordering::Relaxed),
            self.skip.load(Ordering::Relaxed),
        )
    }
}

impl<W: Write + Send, Rf: paraseq::Record> ParallelProcessor<Rf> for Encoder<W> {
    fn process_record(&mut self, record: Rf) -> paraseq::Result<()> {
        let seq = record.seq();
        let rec = SequencingRecordBuilder::default()
            .s_seq(&seq)
            .opt_s_qual(record.qual())
            .s_header(record.id())
            .build()
            .map_err(IntoParaseqError::into_paraseq_error)?;
        self.push(rec).map_err(IntoParaseqError::into_paraseq_error)
    }
    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.batch_complete()
            .map_err(IntoParaseqError::into_paraseq_error)
    }
    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        self.thread_complete()
            .map_err(IntoParaseqError::into_paraseq_error)
    }
}

impl<W: Write + Send, Rf: paraseq::Record> PairedParallelProcessor<Rf> for Encoder<W> {
    fn process_record_pair(&mut self, record1: Rf, record2: Rf) -> paraseq::Result<()> {
        let s_seq = record1.seq();
        let x_seq = record2.seq();
        let rec = SequencingRecordBuilder::default()
            .s_seq(&s_seq)
            .opt_s_qual(record1.qual())
            .s_header(record1.id())
            .x_seq(&x_seq)
            .opt_x_qual(record2.qual())
            .x_header(record2.id())
            .build()
            .map_err(IntoParaseqError::into_paraseq_error)?;
        self.push(rec).map_err(IntoParaseqError::into_paraseq_error)
    }
    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.batch_complete()
            .map_err(IntoParaseqError::into_paraseq_error)
    }
    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        self.thread_complete()
            .map_err(IntoParaseqError::into_paraseq_error)
    }
}
impl<W: Write + Send> binseq::ParallelProcessor for Encoder<W> {
    fn process_record<R: binseq::BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        let rec = if self.t_writer.is_paired() {
            SequencingRecordBuilder::default()
                .s_seq(record.sseq())
                .opt_s_qual(record.has_quality().then(|| record.squal()))
                .s_header(record.sheader())
                .x_seq(record.xseq())
                .opt_x_qual(record.has_quality().then(|| record.xqual()))
                .x_header(record.xheader())
                .build()?
        } else {
            SequencingRecordBuilder::default()
                .s_seq(record.sseq())
                .opt_s_qual(record.has_quality().then(|| record.squal()))
                .s_header(record.sheader())
                .build()?
        };
        self.push(rec)
    }
    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        self.batch_complete()
    }
    fn on_thread_complete(&mut self) -> binseq::Result<()> {
        self.thread_complete()
    }
}
