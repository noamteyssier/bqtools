use std::{
    fs::File,
    io::{BufWriter, Write},
    sync::{Arc, Mutex},
};

use anyhow::Result;
use binseq::ParallelProcessor;
use log::trace;

use super::RecordPair;
use crate::{
    cli::FileFormat,
    commands::{
        decode::{fill_qual, write_record, SeqRead},
        pipe::utils::name_fifo,
    },
};

// The `Arc` is required by the `Clone` bound on `ParallelProcessor`.
#[derive(Clone)]
pub struct PipeProcessor {
    writer: Arc<Mutex<BufWriter<File>>>,
    local: Vec<u8>,
    format: FileFormat,
    pair: RecordPair,
    fallback: Vec<u8>,
}
impl PipeProcessor {
    pub fn new(basename: &str, pid: usize, format: FileFormat, pair: RecordPair) -> Result<Self> {
        let path = name_fifo(basename, pid, pair, format);
        let file = File::options().write(true).open(&path)?;
        trace!("Opened writer at FIFO path: {path}");
        Ok(Self {
            writer: Arc::new(Mutex::new(BufWriter::new(file))),
            local: Vec::new(),
            format,
            pair,
            fallback: Vec::new(),
        })
    }
}
impl ParallelProcessor for PipeProcessor {
    fn process_record<R: binseq::BinseqRecord>(&mut self, record: R) -> binseq::Result<()> {
        let (header, seq, qual, len) = match self.pair {
            RecordPair::Unpaired | RecordPair::R1 => (
                record.sheader(),
                record.sseq(),
                record.has_quality().then(|| record.squal()),
                record.slen() as usize,
            ),
            RecordPair::R2 => (
                record.xheader(),
                record.xseq(),
                record.has_quality().then(|| record.xqual()),
                record.xlen() as usize,
            ),
        };
        // handle missing quality if record has no quality
        let qual = qual.unwrap_or_else(|| fill_qual(&mut self.fallback, len));
        write_record(
            &mut self.local,
            SeqRead { header, seq, qual },
            None,
            self.format,
        )?;
        Ok(())
    }
    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        {
            let mut lock = self.writer.lock().unwrap();
            lock.write_all(&self.local)?;
            lock.flush()?;
        }
        self.local.clear();
        Ok(())
    }
}
