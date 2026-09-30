use anyhow::Result;
use binseq::BinseqRecord;
use serde::Serialize;
use serde_json::{json, Value};
use std::{io::Write, path::Path};

use super::{
    report::{add_assign, stats, table, write_tsv, Hist, Pair},
    QualAbundance, DEFAULT_QUAL_ABUNDANCE, PHRED_OFFSET,
};

#[derive(Serialize)]
struct SeqQualityRecord {
    qual: usize,
    count: usize,
}

#[derive(Clone)]
struct QualHistogram {
    inner: QualAbundance,
}
impl Default for QualHistogram {
    fn default() -> Self {
        Self {
            inner: DEFAULT_QUAL_ABUNDANCE,
        }
    }
}
impl QualHistogram {
    #[allow(clippy::cast_sign_loss)]
    fn push(&mut self, qual: &[u8]) {
        if qual.is_empty() {
            return;
        }
        let total: usize = qual
            .iter()
            .map(|x| x.saturating_sub(PHRED_OFFSET) as usize)
            .sum();
        let binned_mean = (total as f64 / qual.len() as f64).round() as usize;
        self.inner[binned_mean.min(self.inner.len() - 1)] += 1;
    }
}
impl Hist for QualHistogram {
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
                .map(|(qual, &count)| SeqQualityRecord { qual, count }),
        )
    }

    fn summary_table(&self) -> Option<String> {
        if self.is_empty() {
            return None;
        }
        let (total, mean, median, _) = stats(&self.inner);
        Some(table(
            &["Metric", "Value"],
            &[
                vec!["Reads".into(), total.to_string()],
                vec!["Mean Quality".into(), format!("{mean:.2}")],
                vec!["Median Quality".into(), median.to_string()],
            ],
        ))
    }

    fn json(&self) -> Option<Value> {
        if self.is_empty() {
            return None;
        }
        let (reads, mean, median, _) = stats(&self.inner);
        Some(json!({"reads": reads, "mean_quality": mean, "median_quality": median}))
    }
}

#[derive(Default, Clone)]
pub struct PerSequenceQuality(Pair<QualHistogram>);
impl PerSequenceQuality {
    pub fn push<R: BinseqRecord>(&mut self, record: &R) {
        self.0.t[0].push(record.squal());
        self.0.t[1].push(record.xqual());
    }

    pub fn sync_final(&mut self) {
        self.0.sync_final();
    }

    pub fn finish(&mut self, outdir: &Path) -> Result<()> {
        self.0.write(outdir, "seq_quality")
    }

    pub fn summarize(&self) -> String {
        self.0.summarize("Per-Sequence Quality")
    }

    pub fn json(&self) -> Value {
        self.0.json()
    }
}

#[cfg(test)]
// Expected values below are exact (small-integer division that lands on a
// representable value), so strict float equality is correct here.
#[allow(clippy::float_cmp)]
mod tests {
    use super::*;

    fn phred(scores: &[u8]) -> Vec<u8> {
        scores.iter().map(|&s| s + PHRED_OFFSET).collect()
    }

    #[test]
    fn starts_empty() {
        assert!(QualHistogram::default().is_empty());
    }

    #[test]
    fn push_ignores_empty_quality() {
        let mut hist = QualHistogram::default();
        hist.push(&[]);
        assert!(hist.is_empty());
    }

    #[test]
    fn push_bins_by_rounded_mean_quality() {
        let mut hist = QualHistogram::default();
        hist.push(&phred(&[10, 20])); // mean 15
        assert!(!hist.is_empty());
        assert_eq!(stats(&hist.inner).0, 1);
        assert_eq!(stats(&hist.inner).1, 15.0);
    }

    #[test]
    fn mean_and_median_over_multiple_reads() {
        let mut hist = QualHistogram::default();
        hist.push(&phred(&[10]));
        hist.push(&phred(&[20]));
        hist.push(&phred(&[30]));
        assert_eq!(stats(&hist.inner).0, 3);
        assert_eq!(stats(&hist.inner).1, 20.0);
        assert_eq!(stats(&hist.inner).2, 20);
    }

    #[test]
    fn summary_table_none_when_empty() {
        assert!(QualHistogram::default().summary_table().is_none());
    }

    #[test]
    fn summary_table_reports_headline_stats() {
        let mut hist = QualHistogram::default();
        hist.push(&phred(&[10]));
        hist.push(&phred(&[30]));
        let summary = hist.summary_table().expect("non-empty histogram");
        assert!(summary.contains("| Reads | 2 |"));
        assert!(summary.contains("| Mean Quality | 20.00 |"));
        assert!(summary.contains("| Median Quality | 30 |"));
    }

    #[test]
    fn ingest_merges_counts_and_zeroes_source() {
        let mut a = QualHistogram::default();
        let mut b = QualHistogram::default();
        a.push(&phred(&[10]));
        b.push(&phred(&[30]));

        a.ingest(&mut b);

        assert_eq!(stats(&a.inner).0, 2);
        assert_eq!(stats(&a.inner).1, 20.0);
        assert!(b.is_empty());
    }
}
