use std::{
    io::stdout,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
};

use anyhow::Result;
use binseq::{BinseqRecord, ParallelProcessor};
use serde::Serialize;

use crate::commands::grep::SimpleRange;

use super::PatternCounter;

#[derive(Serialize)]
pub struct PatternCountResult<'a> {
    /// The name of the pattern OR the pattern if anonymous
    name: &'a str,
    /// Number of sequences containing the pattern
    count: usize,
    /// Fraction of total sequences containing the pattern
    frac_total: f64,
}
impl<'a> PatternCountResult<'a> {
    pub fn new(name: &'a str, count: usize, total: usize) -> Self {
        Self {
            name,
            count,
            frac_total: if total > 0 {
                count as f64 / total as f64
            } else {
                0.0
            },
        }
    }
}

#[derive(Clone)]
pub struct PatternCountProcessor {
    counter: PatternCounter,
    range: SimpleRange,
    header: bool,

    local_pattern_count: Vec<usize>,
    local_total: usize, // total number of reads processed (not just matches)

    /// Global values
    global_pattern_count: Arc<Mutex<Vec<usize>>>,
    global_total: Arc<AtomicUsize>, // total number of reads processed
}
impl PatternCountProcessor {
    pub fn new(counter: PatternCounter, range: SimpleRange, header: bool) -> Self {
        let num_patterns = counter.num_patterns();
        Self {
            counter,
            range,
            header,
            local_pattern_count: vec![0; num_patterns],
            local_total: 0,
            global_pattern_count: Arc::new(Mutex::new(vec![0; num_patterns])),
            global_total: Arc::new(AtomicUsize::new(0)),
        }
    }
    pub fn pprint_pattern_counts(&self) -> Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_writer(stdout());

        let total_records = self.global_total.load(Ordering::Relaxed);
        let counts = self.global_pattern_count.lock().unwrap();
        for (name, count) in self.counter.pattern_names().iter().zip(counts.iter()) {
            writer.serialize(PatternCountResult::new(name, *count, total_records))?;
        }
        writer.flush()?;
        Ok(())
    }
}
impl ParallelProcessor for PatternCountProcessor {
    fn process_record<B: BinseqRecord>(&mut self, record: B) -> binseq::Result<()> {
        let (primary, extended) = if self.header {
            (record.sheader(), record.xheader())
        } else {
            (
                self.range.slice(record.sseq()),
                self.range.slice(record.xseq()),
            )
        };

        self.counter
            .count_patterns(primary, extended, &mut self.local_pattern_count);
        self.local_total += 1;
        Ok(())
    }

    fn on_batch_complete(&mut self) -> binseq::Result<()> {
        // update the local and global pattern counts
        let mut global = self.global_pattern_count.lock().unwrap();
        for (local, global) in self.local_pattern_count.iter_mut().zip(global.iter_mut()) {
            *global += std::mem::take(local);
        }

        // update the local and global total records processed
        self.global_total
            .fetch_add(std::mem::take(&mut self.local_total), Ordering::Relaxed);

        Ok(())
    }
}
