use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use binseq::prelude::*;
use binseq::Result;
use rand::{RngExt, SeedableRng};
use std::sync::Mutex;

use super::{fill_qual, write_record_pair, SplitWriter};
use crate::cli::{FileFormat, Mate};

/// A struct for decoding BINSEQ data back to FASTQ format.
#[derive(Clone)]
pub struct Decoder {
    /// Local write buffers
    mixed: Vec<u8>, // General purpose, interleaved or singlets
    left: Vec<u8>, // Used when writing pairs of files (R1/R2)
    right: Vec<u8>,

    /// Local count of records
    local_count: usize,
    /// Quality buffer (primary)
    squal: Vec<u8>,
    /// Quality buffer (extended)
    xqual: Vec<u8>,

    /// Options
    format: FileFormat,
    mate: Mate,
    is_split: bool,
    /// Optional `(fraction, seed)` keep-filter
    sample: Option<(f64, u64)>,

    /// Global values
    global_writer: Arc<Mutex<SplitWriter>>,
    num_records: Arc<AtomicUsize>,
}

impl Decoder {
    pub fn new(
        writer: SplitWriter,
        format: FileFormat,
        mate: Mate,
        sample: Option<(f64, u64)>,
    ) -> Self {
        Decoder {
            mixed: Vec::new(),
            left: Vec::new(),
            right: Vec::new(),
            local_count: 0,
            squal: Vec::new(),
            xqual: Vec::new(),
            format,
            mate,
            is_split: writer.is_split(),
            sample,
            global_writer: Arc::new(Mutex::new(writer)),
            num_records: Arc::new(AtomicUsize::new(0)),
        }
    }

    pub fn num_records(&self) -> usize {
        self.num_records.load(Ordering::Relaxed)
    }
}

impl ParallelProcessor for Decoder {
    fn process_record<B: BinseqRecord>(&mut self, record: B) -> Result<()> {
        // Keep/drop is a pure function of `(seed, record index)`, so the sample is
        // reproducible regardless of thread count or batch boundaries.
        if let Some((fraction, seed)) = self.sample {
            let index = record.index();
            if !rand::rngs::SmallRng::seed_from_u64(seed.wrapping_add(index)).random_bool(fraction)
            {
                return Ok(());
            }
        }

        let sbuf = record.sseq();
        let xbuf = record.xseq();

        // decode sequences
        let squal = if record.has_quality() {
            record.squal()
        } else {
            fill_qual(&mut self.squal, sbuf.len())
        };

        let xqual = if record.is_paired() && record.has_quality() {
            record.xqual()
        } else {
            fill_qual(&mut self.xqual, xbuf.len())
        };

        write_record_pair(
            &mut self.left,
            &mut self.right,
            &mut self.mixed,
            self.mate,
            self.is_split,
            sbuf,
            squal,
            record.sheader(),
            xbuf,
            xqual,
            record.xheader(),
            self.format,
        )?;

        self.local_count += 1;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        // Lock the mutex to write to the global buffer
        {
            let mut writer = self.global_writer.lock().unwrap();
            writer.write_batch(&self.left, &self.right, &self.mixed)?;
        }
        self.num_records
            .fetch_add(self.local_count, Ordering::Relaxed);

        // Clear the local buffer and reset the local record count
        self.mixed.clear();
        self.left.clear();
        self.right.clear();
        self.local_count = 0;
        Ok(())
    }
}
