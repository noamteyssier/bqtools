use anyhow::Result;
use binseq::BinseqRecord;
use serde::Serialize;
use std::{io::Write, path::Path};

use super::report::{add_assign, stats, table, write_tsv, Hist, Pair};

/// Percentage bins: 0..=100
const NUM_GC_BINS: usize = 101;

fn is_gc(base: u8) -> bool {
    matches!(base, b'G' | b'g' | b'C' | b'c')
}

#[derive(Serialize)]
struct GcContentRecord {
    pct_gc: usize,
    count: usize,
}

#[derive(Clone)]
struct GcHistogram {
    inner: [usize; NUM_GC_BINS],
}
impl Default for GcHistogram {
    fn default() -> Self {
        Self {
            inner: [0; NUM_GC_BINS],
        }
    }
}
impl GcHistogram {
    /// Bin a whole read by the percentage of G/C bases it contains.
    #[allow(clippy::cast_sign_loss)]
    fn push(&mut self, seq: &[u8]) {
        if seq.is_empty() {
            return;
        }
        let gc = seq.iter().copied().filter(|&b| is_gc(b)).count();
        let pct_gc = ((gc as f64 / seq.len() as f64) * 100.0).round() as usize;
        self.inner[pct_gc.min(self.inner.len() - 1)] += 1;
    }
}
impl Hist for GcHistogram {
    fn is_empty(&self) -> bool {
        self.inner.iter().all(|&c| c == 0)
    }

    fn ingest(&mut self, other: &mut Self) {
        add_assign(&mut self.inner, &mut other.inner);
    }

    fn serialize_to<W: Write>(&self, wtr: &mut W) -> Result<()> {
        if self.is_empty() {
            return Ok(());
        }

        write_tsv(
            wtr,
            self.inner
                .iter()
                .enumerate()
                .map(|(pct_gc, &count)| GcContentRecord { pct_gc, count }),
        )
    }

    fn summary_table(&self) -> Option<String> {
        if self.is_empty() {
            return None;
        }
        let (total, mean, median, mode) = stats(&self.inner);
        Some(table(
            &["Metric", "Value"],
            &[
                vec!["Reads".into(), total.to_string()],
                vec!["Mean GC%".into(), format!("{mean:.2}%")],
                vec!["Median GC%".into(), format!("{median}%")],
                vec!["Mode GC%".into(), format!("{mode}%")],
            ],
        ))
    }
}

#[derive(Default, Clone)]
pub struct PerSequenceGcContent(Pair<GcHistogram>);
impl PerSequenceGcContent {
    pub fn push<R: BinseqRecord>(&mut self, record: &R) {
        self.0.t[0].push(record.sseq());
        self.0.t[1].push(record.xseq());
    }

    pub fn sync_final(&mut self) {
        self.0.sync_final();
    }

    pub fn finish(&mut self, outdir: &Path) -> Result<()> {
        self.0.write(outdir, "gc_content")
    }

    pub fn summarize(&self) -> String {
        self.0.summarize("Per-Sequence GC Content")
    }
}

#[cfg(test)]
// Expected values below are exact (small-integer division that lands on a
// representable value), so strict float equality is correct here.
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    #[test]
    fn starts_empty() {
        assert!(GcHistogram::default().is_empty());
    }

    #[test]
    fn is_gc_matches_c_and_g_case_insensitively() {
        assert!(is_gc(b'G'));
        assert!(is_gc(b'g'));
        assert!(is_gc(b'C'));
        assert!(is_gc(b'c'));
        assert!(!is_gc(b'A'));
        assert!(!is_gc(b'T'));
        assert!(!is_gc(b'N'));
    }

    #[test]
    fn push_ignores_empty_sequence() {
        let mut hist = GcHistogram::default();
        hist.push(b"");
        assert!(hist.is_empty());
    }

    #[test]
    fn push_bins_by_rounded_gc_percentage() {
        let mut hist = GcHistogram::default();
        hist.push(b"GCAT"); // 2/4 = 50% GC
        assert!(!hist.is_empty());
        assert_eq!(stats(&hist.inner).0, 1);
        assert_eq!(stats(&hist.inner).1, 50.0);
    }

    #[test]
    fn mean_median_mode_over_multiple_reads() {
        let mut hist = GcHistogram::default();
        hist.push(b"AAAA"); // 0% GC
        hist.push(b"AAAA"); // 0% GC
        hist.push(b"GGGG"); // 100% GC
        assert_eq!(stats(&hist.inner).0, 3);
        assert!((stats(&hist.inner).1 - 33.333_333_333_333_336).abs() < 1e-9);
        assert_eq!(stats(&hist.inner).2, 0);
        assert_eq!(stats(&hist.inner).3, 0);
    }

    #[test]
    fn summary_table_none_when_empty() {
        assert!(GcHistogram::default().summary_table().is_none());
    }

    #[test]
    fn summary_table_reports_headline_stats() {
        let mut hist = GcHistogram::default();
        hist.push(b"GCGC"); // 100% GC
        hist.push(b"ATAT"); // 0% GC
        let summary = hist.summary_table().expect("non-empty histogram");
        assert!(summary.contains("| Reads | 2 |"));
        assert!(summary.contains("| Mean GC% | 50.00% |"));
    }

    #[test]
    fn ingest_merges_counts_and_zeroes_source() {
        let mut a = GcHistogram::default();
        let mut b = GcHistogram::default();
        a.push(b"AAAA"); // 0%
        b.push(b"GGGG"); // 100%

        a.ingest(&mut b);

        assert_eq!(stats(&a.inner).0, 2);
        assert_eq!(stats(&a.inner).1, 50.0);
        assert!(b.is_empty());
    }
}
