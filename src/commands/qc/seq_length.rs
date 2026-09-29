use std::{io::Write, path::Path, sync::Arc};

use anyhow::Result;
use binseq::BinseqRecord;
use serde::Serialize;
use std::sync::Mutex;

use super::report::{add_assign, stats, table, write_tsv};
use crate::commands::{match_output, qc::modules::QcModule};

const SEQ_LENGTH_PRIMARY_PATH: &str = "seq_length_R1.tsv";
const SEQ_LENGTH_EXTENDED_PATH: &str = "seq_length_R2.tsv";

#[derive(Serialize)]
pub struct SeqLenRecord {
    len: usize,
    count: usize,
}

#[derive(Debug, Clone, Default)]
pub struct SeqLenHistogram {
    /// Indexed directly by sequence length
    inner: Vec<usize>,
}
impl SeqLenHistogram {
    fn is_empty(&self) -> bool {
        self.inner.iter().copied().sum::<usize>() == 0
    }
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

    fn min_len(&self) -> Option<usize> {
        self.inner.iter().position(|&c| c > 0)
    }

    fn max_len(&self) -> Option<usize> {
        self.inner.iter().rposition(|&c| c > 0)
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
pub struct SequenceLengthDistribution {
    /// thread - sequence length distribution (primary)
    t_slen: SeqLenHistogram,
    /// thread - sequence length distribution (extended)
    t_xlen: SeqLenHistogram,

    /// global - sequence length distribution (primary)
    slen: Arc<Mutex<SeqLenHistogram>>,
    /// global - sequence length distribution (extended)
    xlen: Arc<Mutex<SeqLenHistogram>>,
}
impl QcModule for SequenceLengthDistribution {
    fn push<R: BinseqRecord>(&mut self, record: &R) {
        self.t_slen.push(record.slen() as usize);
        self.t_xlen.push(record.xlen() as usize);
    }

    fn sync_final(&mut self) {
        self.slen.lock().unwrap().ingest(&mut self.t_slen);
        self.xlen.lock().unwrap().ingest(&mut self.t_xlen);
    }

    fn finish<P: AsRef<Path>>(&mut self, outdir: P) -> Result<()> {
        if !outdir.as_ref().exists() {
            std::fs::create_dir_all(outdir.as_ref())?;
        }

        let write_to = |hist: &SeqLenHistogram, primary: bool| -> Result<()> {
            if hist.is_empty() {
                return Ok(());
            }
            let mut handle = if primary {
                match_output(Some(outdir.as_ref().join(SEQ_LENGTH_PRIMARY_PATH)))
            } else {
                match_output(Some(outdir.as_ref().join(SEQ_LENGTH_EXTENDED_PATH)))
            }?;
            hist.serialize_to(&mut handle)
        };

        write_to(&self.slen.lock().unwrap(), true)?;
        write_to(&self.xlen.lock().unwrap(), false)?;

        Ok(())
    }

    fn summarize(&self) -> String {
        let primary = self.slen.lock().unwrap().summary_table();
        let extended = self.xlen.lock().unwrap().summary_table();
        super::report::dual_section("Sequence Length Distribution", primary, extended)
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
