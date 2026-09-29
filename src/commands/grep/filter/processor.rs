use crate::{
    cli::{FileFormat, Mate},
    commands::{
        decode::{fill_qual, write_record_pair, SplitWriter},
        grep::{color::write_colored_record_pair, Engine, SimpleRange, Spans},
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
    mixed: Vec<u8>, // General purpose, interleaved or singlets
    left: Vec<u8>, // Used when writing pairs of files (R1/R2)
    right: Vec<u8>,

    /// Quality buffers
    squal: Vec<u8>,
    xqual: Vec<u8>,

    /// Write Options
    format: FileFormat,
    mate: Option<Mate>,
    is_split: bool,
    color: bool,

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
            mixed: Vec::new(),
            left: Vec::new(),
            right: Vec::new(),
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
            color,
            is_split: writer.is_split(),
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

            if self.color {
                write_colored_record_pair(
                    &mut self.mixed,
                    self.mate,
                    sbuf,
                    squal,
                    record.sheader(),
                    xbuf,
                    xqual,
                    record.xheader(),
                    &mut self.spans.primary,
                    &mut self.spans.secondary,
                    self.format,
                )
            } else {
                write_record_pair(
                    &mut self.left,
                    &mut self.right,
                    &mut self.mixed,
                    self.mate.unwrap_or(Mate::One),
                    self.is_split,
                    sbuf,
                    squal,
                    record.sheader(),
                    xbuf,
                    xqual,
                    record.xheader(),
                    self.format,
                )
            }?;
        }

        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        // Lock the mutex to write to the global buffer
        if !self.count {
            let mut writer = self.global_writer.lock().unwrap();
            writer.write_batch(&self.left, &self.right, &self.mixed)?;
        }

        // Clear the local buffer and reset the local record count
        self.mixed.clear();
        self.left.clear();
        self.right.clear();

        // Increment the global count and reset local
        self.global_count
            .fetch_add(std::mem::take(&mut self.local_count), Ordering::Relaxed);

        // Increment the global total and reset local
        self.global_total
            .fetch_add(std::mem::take(&mut self.local_total), Ordering::Relaxed);

        Ok(())
    }
}
