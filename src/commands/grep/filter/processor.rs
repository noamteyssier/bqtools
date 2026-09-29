use crate::{
    cli::{FileFormat, Mate},
    commands::{
        decode::{fill_qual, Batch, SeqRead, SplitWriter},
        grep::{Engine, SimpleRange, Spans},
    },
};
use binseq::prelude::*;
use fixedbitset::FixedBitSet;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, Mutex,
};

#[derive(Clone)]
#[allow(clippy::struct_excessive_bools)]
pub struct FilterProcessor {
    engine: Engine,

    /// Match logic (true = AND, false = OR)
    and_logic: bool,

    /// Invert the pattern selection
    invert: bool,

    /// Only count the number of matches
    count: bool,

    /// Show count as fraction of total records
    frac: bool,

    /// Match within range
    range: SimpleRange,

    /// Match against the sequence header instead of the sequence
    header: bool,

    /// Local count
    local_count: usize,

    /// Local total records processed
    local_total: usize,

    /// Local pattern hits, and where they are (only filled for colored output)
    bits: FixedBitSet,
    spans: Spans,
    collect_spans: bool,

    /// Local write buffers
    batch: Batch,

    /// Quality buffers
    squal: Vec<u8>,
    xqual: Vec<u8>,

    /// Write Options
    format: FileFormat,
    mate: Option<Mate>,

    /// Global values
    global_writer: Arc<Mutex<SplitWriter>>,
    global_count: Arc<AtomicUsize>,
    global_total: Arc<AtomicUsize>,
}
impl FilterProcessor {
    #[allow(clippy::fn_params_excessive_bools)]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        engine: Engine,
        and_logic: bool,
        invert: bool,
        count: bool,
        frac: bool,
        range: SimpleRange,
        header: bool,
        writer: SplitWriter,
        format: FileFormat,
        mate: Option<Mate>,
        color: bool,
    ) -> Self {
        Self {
            batch: Batch::new(&writer),
            squal: Vec::new(),
            xqual: Vec::new(),
            bits: engine.bitset(),
            spans: Spans::default(),
            // hit positions are only drawn for colored, non-inverted, written matches
            collect_spans: color && !invert && !count,
            engine,
            and_logic,
            invert,
            count,
            frac,
            range,
            header,
            format,
            mate,
            global_writer: Arc::new(Mutex::new(writer)),
            local_count: 0,
            local_total: 0,
            global_count: Arc::new(AtomicUsize::new(0)),
            global_total: Arc::new(AtomicUsize::new(0)),
        }
    }
    fn pattern_match(&mut self, sbuf: &[u8], xbuf: &[u8]) -> bool {
        let (primary, extended) = (self.range.slice(sbuf), self.range.slice(xbuf));

        let found = match (self.and_logic, self.collect_spans) {
            // AND gives up at the first missing pattern; positions are only
            // worth collecting once everything is known to hit
            (true, _) => self.engine.all(primary, extended, &mut self.bits),
            // OR without positions only needs to know whether anything hits
            (false, false) => self.engine.any(primary, extended),
            // a separate yes/no pass would scan matching records twice
            (false, true) => {
                self.collect(primary, extended);
                !self.bits.is_clear()
            }
        };
        if found && self.and_logic && self.collect_spans {
            self.collect(primary, extended);
        }
        found != self.invert
    }

    /// Finds every pattern and where it hits, for colored output.
    fn collect(&mut self, primary: &[u8], extended: &[u8]) {
        self.bits.clear();
        self.spans.clear();
        self.engine
            .hit(primary, extended, &mut self.bits, Some(&mut self.spans));
        self.spans.shift(self.range.offset());
    }
    pub fn pprint_counts(&self) {
        let count = self.global_count.load(Ordering::Relaxed);
        if self.frac {
            let total = self.global_total.load(Ordering::Relaxed);
            let frac = if total > 0 {
                count as f64 / total as f64
            } else {
                0.0
            };
            println!("count\ttotal\tfrac");
            println!("{count}\t{total}\t{frac:.4}");
        } else {
            println!("{count}");
        }
    }
}

impl ParallelProcessor for FilterProcessor {
    fn process_record<B: BinseqRecord>(&mut self, record: B) -> binseq::Result<()> {
        self.local_total += 1;

        let sbuf = record.sseq();
        let xbuf = record.xseq();
        let matched = if self.header {
            self.pattern_match(record.sheader(), record.xheader())
        } else {
            self.pattern_match(sbuf, xbuf)
        };
        if matched {
            self.local_count += 1;
            if self.count {
                // No further processing needed
                return Ok(());
            }

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
                self.mate.unwrap_or(Mate::One),
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
                self.collect_spans
                    .then_some([&mut self.spans.primary, &mut self.spans.secondary]),
                self.format,
            )?;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        // Lock the mutex to write to the global buffer
        if !self.count {
            let mut writer = self.global_writer.lock().unwrap();
            writer.write_batch(&self.batch)?;
        }

        // Clear the local buffer and reset the local record count
        self.batch.clear();

        // Increment the global count and reset local
        self.global_count
            .fetch_add(std::mem::take(&mut self.local_count), Ordering::Relaxed);

        // Increment the global total and reset local
        self.global_total
            .fetch_add(std::mem::take(&mut self.local_total), Ordering::Relaxed);

        Ok(())
    }
}
