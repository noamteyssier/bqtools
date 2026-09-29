use std::{io::Write, path::Path};

use anyhow::Result;
use binseq::BinseqRecord;
use serde::Serialize;

use super::report::{add_assign, stats, table, write_tsv, Hist, Pair};

#[derive(Serialize)]
struct SeqLenRecord {
    len: usize,
    count: usize,
}

#[derive(Debug, Clone, Default)]
struct SeqLenHistogram {
    /// Indexed directly by sequence length
    inner: Vec<usize>,
}
impl SeqLenHistogram {
    fn len(&self) -> usize {
        self.inner.len()
    }
    /// Track a single read's length
    fn push(&mut self, len: usize) {
        if len == 0 {
            return;
        }
        if self.inner.len() <= len {
            self.inner.resize(len + 1, 0);
        }
        self.inner[len] += 1;
    }

    fn min_len(&self) -> Option<usize> {
        self.inner.iter().position(|&c| c > 0)
    }

    fn max_len(&self) -> Option<usize> {
        self.inner.iter().rposition(|&c| c > 0)
    }
}
impl Hist for SeqLenHistogram {
    fn is_empty(&self) -> bool {
        self.inner.iter().all(|&c| c == 0)
    }

    fn ingest(&mut self, other: &mut Self) {
        if self.len() < other.len() {
            self.inner.resize(other.len(), 0);
        }
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
                .filter(|(_, &count)| count > 0)
                .map(|(len, &count)| SeqLenRecord { len, count }),
        )
    }

    fn summary_table(&self) -> Option<String> {
        if self.is_empty() {
            return None;
        }
        let (total, mean, _, mode) = stats(&self.inner);
        Some(table(
            &["Metric", "Value"],
            &[
                vec!["Reads".into(), total.to_string()],
                vec!["Min Length".into(), self.min_len().unwrap_or(0).to_string()],
                vec!["Max Length".into(), self.max_len().unwrap_or(0).to_string()],
                vec!["Mean Length".into(), format!("{mean:.2}")],
                vec!["Mode Length".into(), mode.to_string()],
            ],
        ))
    }
}

#[derive(Clone, Default)]
pub struct SequenceLengthDistribution(Pair<SeqLenHistogram>);
impl SequenceLengthDistribution {
    pub fn push<R: BinseqRecord>(&mut self, record: &R) {
        self.0.t[0].push(record.slen() as usize);
        self.0.t[1].push(record.xlen() as usize);
    }

    pub fn sync_final(&mut self) {
        self.0.sync_final();
    }

    pub fn finish(&mut self, outdir: &Path) -> Result<()> {
        self.0.write(outdir, "seq_length")
    }

    pub fn summarize(&self) -> String {
        self.0.summarize("Sequence Length Distribution")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn starts_empty() {
        assert!(SeqLenHistogram::default().is_empty());
    }

    #[test]
    fn push_ignores_zero_length() {
        let mut hist = SeqLenHistogram::default();
        hist.push(0);
        assert!(hist.is_empty());
    }

    #[test]
    fn push_tracks_length_counts() {
        let mut hist = SeqLenHistogram::default();
        hist.push(100);
        hist.push(100);
        hist.push(150);
        assert!(!hist.is_empty());
        assert_eq!(stats(&hist.inner).0, 3);
        assert_eq!(hist.min_len(), Some(100));
        assert_eq!(hist.max_len(), Some(150));
        assert_eq!(stats(&hist.inner).3, 100);
    }

    #[test]
    fn mean_is_weighted_by_count() {
        let mut hist = SeqLenHistogram::default();
        hist.push(100);
        hist.push(100);
        hist.push(200);
        assert!((stats(&hist.inner).1 - 133.333_333_333_333_33).abs() < 1e-6);
    }

    #[test]
    fn summary_table_none_when_empty() {
        assert!(SeqLenHistogram::default().summary_table().is_none());
    }

    #[test]
    fn summary_table_reports_headline_stats() {
        let mut hist = SeqLenHistogram::default();
        hist.push(50);
        hist.push(50);
        let summary = hist.summary_table().expect("non-empty histogram");
        assert!(summary.contains("| Reads | 2 |"));
        assert!(summary.contains("| Min Length | 50 |"));
        assert!(summary.contains("| Max Length | 50 |"));
        assert!(summary.contains("| Mean Length | 50.00 |"));
        assert!(summary.contains("| Mode Length | 50 |"));
    }

    #[test]
    fn ingest_merges_counts_and_zeroes_source() {
        let mut a = SeqLenHistogram::default();
        let mut b = SeqLenHistogram::default();
        a.push(100);
        b.push(200);

        a.ingest(&mut b);

        assert_eq!(stats(&a.inner).0, 2);
        assert_eq!(a.min_len(), Some(100));
        assert_eq!(a.max_len(), Some(200));
        assert!(b.is_empty());
    }
}
