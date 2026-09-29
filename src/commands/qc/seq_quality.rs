use anyhow::Result;
use binseq::BinseqRecord;
use serde::Serialize;
use std::sync::Mutex;
use std::{io::Write, ops::Div, path::Path, sync::Arc};

use super::{
    report::{add_assign, stats, table, write_tsv},
    QualAbundance, DEFAULT_QUAL_ABUNDANCE, PHRED_OFFSET,
};
use crate::commands::{match_output, qc::modules::QcModule};

const SEQ_QUALITY_PRIMARY_PATH: &str = "seq_quality_R1.tsv";
const SEQ_QUALITY_EXTENDED_PATH: &str = "seq_quality_R2.tsv";

#[derive(Serialize)]
struct SeqQualityRecord {
    qual: usize,
    count: usize,
}

#[derive(Clone)]
pub struct QualHistogram {
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
    fn is_empty(&self) -> bool {
        self.inner.iter().copied().sum::<usize>() == 0
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

    #[allow(clippy::cast_sign_loss)]
    fn push(&mut self, qual: &[u8]) {
        if qual.is_empty() {
            return;
        }
        let total: usize = qual
            .iter()
            .map(|x| x.saturating_sub(PHRED_OFFSET) as usize)
            .sum();
        let binned_mean = (total as f64).div(&(qual.len() as f64)).round() as usize;
        self.inner[binned_mean.min(self.inner.len() - 1)] += 1;
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
}

#[derive(Default, Clone)]
pub struct PerSequenceQuality {
    t_seq_squal: QualHistogram,
    t_seq_xqual: QualHistogram,

    seq_squal: Arc<Mutex<QualHistogram>>,
    seq_xqual: Arc<Mutex<QualHistogram>>,
}
impl QcModule for PerSequenceQuality {
    fn push<R: BinseqRecord>(&mut self, record: &R) {
        self.t_seq_squal.push(record.squal());
        self.t_seq_xqual.push(record.xqual());
    }

    fn sync_final(&mut self) {
        self.seq_squal.lock().unwrap().ingest(&mut self.t_seq_squal);
        self.seq_xqual.lock().unwrap().ingest(&mut self.t_seq_xqual);
    }

    fn finish<P: AsRef<Path>>(&mut self, outdir: P) -> Result<()> {
        if !outdir.as_ref().exists() {
            std::fs::create_dir_all(outdir.as_ref())?;
        }

        let write_to = |seq_qual: &QualHistogram, primary: bool| -> Result<()> {
            if seq_qual.is_empty() {
                return Ok(());
            }
            let mut handle = if primary {
                match_output(Some(outdir.as_ref().join(SEQ_QUALITY_PRIMARY_PATH)))
            } else {
                match_output(Some(outdir.as_ref().join(SEQ_QUALITY_EXTENDED_PATH)))
            }?;
            seq_qual.serialize_to(&mut handle)
        };

        write_to(&self.seq_squal.lock().unwrap(), true)?;
        write_to(&self.seq_xqual.lock().unwrap(), false)?;

        Ok(())
    }

    fn summarize(&self) -> String {
        let primary = self.seq_squal.lock().unwrap().summary_table();
        let extended = self.seq_xqual.lock().unwrap().summary_table();
        super::report::dual_section("Per-Sequence Quality", primary, extended)
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
