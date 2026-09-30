use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use binseq::prelude::*;
use binseq::Result;
use rand::{RngExt, SeedableRng};
use std::sync::Mutex;

use super::{fill_qual, Batch, SeqRead, SplitWriter};
use crate::cli::{FileFormat, Mate};

/// A struct for decoding BINSEQ data back to FASTQ format.
#[derive(Clone)]
pub struct Decoder {
    /// Local write buffers
    batch: Batch,

    /// Local count of records
    local_count: usize,
    /// Quality buffer (primary)
    squal: Vec<u8>,
    /// Quality buffer (extended)
    xqual: Vec<u8>,

    /// Options
    format: FileFormat,
    mate: Mate,
    /// Optional keep-filter
    sample: Option<Sample>,

    /// Global values
    global_writer: Arc<Mutex<SplitWriter>>,
    num_records: Arc<AtomicUsize>,
}

impl Decoder {
    pub fn new(
        writer: SplitWriter,
        format: FileFormat,
        mate: Mate,
        sample: Option<Sample>,
    ) -> Self {
        Decoder {
            batch: Batch::new(&writer),
            local_count: 0,
            squal: Vec::new(),
            xqual: Vec::new(),
            format,
            mate,
            sample,
            global_writer: Arc::new(Mutex::new(writer)),
            num_records: Arc::new(AtomicUsize::new(0)),
        }
    }

    pub fn num_records(&self) -> usize {
        self.num_records.load(Ordering::Relaxed)
    }
}

/// Keep/drop decision as a pure function of `(seed, record index)`.
fn keep(index: u64, fraction: f64, seed: u64) -> bool {
    rand::rngs::SmallRng::seed_from_u64(seed.wrapping_add(index)).random_bool(fraction)
}

/// Record sampling strategy.
#[derive(Clone)]
pub enum Sample {
    /// Keep each record independently with probability `fraction`.
    Fraction { fraction: f64, seed: u64 },
    /// Keep exactly these record indices.
    Exact(Arc<hashbrown::HashSet<u64>>),
}
impl Sample {
    /// Choose exactly `min(n, range.len())` indices from `range`.
    pub fn exact(n: usize, range: std::ops::Range<usize>, seed: u64) -> Self {
        let mut rng = rand::rngs::SmallRng::seed_from_u64(seed);
        let n = n.min(range.len());
        let set = rand::seq::index::sample(&mut rng, range.len(), n)
            .into_iter()
            .map(|i| (range.start + i) as u64)
            .collect();
        Self::Exact(Arc::new(set))
    }

    pub fn keep(&self, index: u64) -> bool {
        match self {
            Self::Fraction { fraction, seed } => keep(index, *fraction, *seed),
            Self::Exact(set) => set.contains(&index),
        }
    }
}

impl ParallelProcessor for Decoder {
    fn process_record<B: BinseqRecord>(&mut self, record: B) -> Result<()> {
        // Keep/drop is a pure function of `(seed, record index)`, so the sample is
        // reproducible regardless of thread count or batch boundaries.
        if let Some(sample) = &self.sample {
            if !sample.keep(record.index()) {
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

        self.batch.push_pair(
            self.mate,
            SeqRead {
                header: record.sheader(),
                seq: sbuf,
                qual: squal,
            },
            SeqRead {
                header: record.xheader(),
                seq: xbuf,
                qual: xqual,
            },
            None,
            self.format,
        )?;

        self.local_count += 1;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> Result<()> {
        // Lock the mutex to write to the global buffer
        {
            let mut writer = self.global_writer.lock().unwrap();
            writer.write_batch(&self.batch)?;
        }
        self.num_records
            .fetch_add(self.local_count, Ordering::Relaxed);

        // Clear the local buffer and reset the local record count
        self.batch.clear();
        self.local_count = 0;
        Ok(())
    }
}
