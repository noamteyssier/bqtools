use crate::{
    cli::{FileFormat, Mate},
    commands::{
        decode::{fill_qual, write_record_pair, SplitWriter},
        grep::{color::write_colored_record_pair, SimpleRange},
    },
};
use binseq::prelude::*;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, Mutex,
};

use super::{MatchRanges, PatternMatch};

#[derive(Clone)]
#[allow(clippy::struct_excessive_bools)]
pub struct FilterProcessor<Pm: PatternMatch> {
    matcher: Pm,

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

    /// Local primary/extended sequence match indices
    smatches: MatchRanges,
    xmatches: MatchRanges,

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
impl<Pm: PatternMatch> FilterProcessor<Pm> {
    #[allow(clippy::fn_params_excessive_bools)]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        matcher: Pm,
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
            smatches: MatchRanges::default(),
            xmatches: MatchRanges::default(),
            matcher,
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
    pub fn clear_matches(&mut self) {
        self.smatches.clear();
        self.xmatches.clear();
    }

    pub fn pattern_match(&mut self, sbuf: &[u8], xbuf: &[u8]) -> bool {
        let (primary, extended) = (self.range.slice(sbuf), self.range.slice(xbuf));

        let found_either = self.matcher.match_either(
            primary,
            extended,
            &mut self.smatches,
            &mut self.xmatches,
            self.and_logic,
        );
        let found_primary = self
            .matcher
            .match_primary(primary, &mut self.smatches, self.and_logic);
        let found_secondary =
            self.matcher
                .match_secondary(extended, &mut self.xmatches, self.and_logic);

        let pred = if self.and_logic {
            found_either && found_primary && found_secondary
        } else {
            !self.smatches.is_empty() || !self.xmatches.is_empty()
        };

        if self.invert {
            self.clear_matches(); // ensure no partial matches are highlighted
            !pred
        } else {
            pred
        }
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

impl<Pm: PatternMatch> ParallelProcessor for FilterProcessor<Pm> {
    fn process_record<B: BinseqRecord>(&mut self, record: B) -> binseq::Result<()> {
        self.clear_matches();
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
                    &mut self.smatches,
                    &mut self.xmatches,
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
