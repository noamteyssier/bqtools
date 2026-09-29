use std::{io::Write, path::Path};

use anyhow::Result;
use binseq::BinseqRecord;
use serde::Serialize;

use super::{
    report::{add_assign, stats, table, write_tsv, Hist, Pair},
    QualAbundance, DEFAULT_QUAL_ABUNDANCE, PHRED_OFFSET,
};

#[derive(Serialize)]
struct BaseQualityRecord {
    pos: usize,
    qual: usize,
    count: usize,
}

#[derive(Debug, Clone, Default)]
struct BaseHistogram {
    /// Outer: position
    /// Inner: quality
    inner: Vec<QualAbundance>,
}
impl BaseHistogram {
    /// Number of positions tracked
    fn len(&self) -> usize {
        self.inner.len()
    }
    /// Track quality score histogram over positions
    fn push(&mut self, qual: &[u8]) {
        if qual.is_empty() {
            return;
        }
        if self.inner.len() <= qual.len() {
            self.inner.resize(qual.len(), DEFAULT_QUAL_ABUNDANCE);
        }
        qual.iter()
            .map(|q| q.saturating_sub(PHRED_OFFSET) as usize)
            .zip(self.inner.iter_mut())
            .for_each(|(q, pos_vec)| {
                pos_vec[q.min(pos_vec.len() - 1)] += 1;
            });
    }

    /// Mean quality score at each position.
    fn position_means(&self) -> Vec<f64> {
        self.inner.iter().map(|counts| stats(counts).1).collect()
    }
}
impl Hist for BaseHistogram {
    /// Checks if empty
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn ingest(&mut self, other: &mut Self) {
        if self.len() < other.len() {
            self.inner.resize(other.len(), DEFAULT_QUAL_ABUNDANCE);
        }
        for (dst, src) in self.inner.iter_mut().zip(&mut other.inner) {
            add_assign(dst, src);
        }
    }

    fn serialize_to<W: Write>(&self, wtr: &mut W) -> Result<()> {
        if self.is_empty() {
            return Ok(());
        }

        write_tsv(
            wtr,
            self.inner.iter().enumerate().flat_map(|(pos, inner)| {
                inner
                    .iter()
                    .enumerate()
                    .filter(|(_, &count)| count > 0)
                    .map(move |(qual, &count)| BaseQualityRecord { pos, qual, count })
            }),
        )
    }

    fn summary_table(&self) -> Option<String> {
        if self.is_empty() {
            return None;
        }

        let mut num = 0usize;
        let mut den = 0usize;
        for counts in &self.inner {
            for (q, &c) in counts.iter().enumerate() {
                num += q * c;
                den += c;
            }
        }
        let overall_mean = if den == 0 {
            0.0
        } else {
            num as f64 / den as f64
        };

        let means = self.position_means();
        let (min_pos, min_mean) = means
            .iter()
            .copied()
            .enumerate()
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap_or((0, 0.0));
        let (max_pos, max_mean) = means
            .iter()
            .copied()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap_or((0, 0.0));

        Some(table(
            &["Metric", "Value"],
            &[
                vec!["Positions".into(), means.len().to_string()],
                vec!["Mean Quality".into(), format!("{overall_mean:.2}")],
                vec![
                    "Lowest Mean Quality".into(),
                    format!("{min_mean:.2} (pos {min_pos})"),
                ],
                vec![
                    "Highest Mean Quality".into(),
                    format!("{max_mean:.2} (pos {max_pos})"),
                ],
            ],
        ))
    }
}

#[derive(Clone, Default)]
pub struct PerBaseSequenceQuality(Pair<BaseHistogram>);
impl PerBaseSequenceQuality {
    pub fn push<R: BinseqRecord>(&mut self, record: &R) {
        self.0.t[0].push(record.squal());
        self.0.t[1].push(record.xqual());
    }

    pub fn sync_final(&mut self) {
        self.0.sync_final();
    }

    pub fn finish(&mut self, outdir: &Path) -> Result<()> {
        self.0.write(outdir, "base_quality")
    }

    pub fn summarize(&self) -> String {
        self.0.summarize("Per-Base Sequence Quality")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn phred(scores: &[u8]) -> Vec<u8> {
        scores.iter().map(|&s| s + PHRED_OFFSET).collect()
    }

    #[test]
    fn starts_empty() {
        assert!(BaseHistogram::default().is_empty());
    }

    #[test]
    fn push_ignores_empty_quality() {
        let mut hist = BaseHistogram::default();
        hist.push(&[]);
        assert!(hist.is_empty());
    }

    #[test]
    fn push_tracks_per_position_quality() {
        let mut hist = BaseHistogram::default();
        hist.push(&phred(&[10, 20, 30]));
        assert!(!hist.is_empty());
        assert_eq!(hist.len(), 3);
        assert_eq!(hist.position_means(), vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn push_accumulates_across_reads() {
        let mut hist = BaseHistogram::default();
        hist.push(&phred(&[10, 10]));
        hist.push(&phred(&[30, 30]));
        assert_eq!(hist.position_means(), vec![20.0, 20.0]);
    }

    #[test]
    fn summary_table_none_when_empty() {
        assert!(BaseHistogram::default().summary_table().is_none());
    }

    #[test]
    fn summary_table_reports_overall_and_extreme_positions() {
        let mut hist = BaseHistogram::default();
        hist.push(&phred(&[10, 40]));
        let summary = hist.summary_table().expect("non-empty histogram");
        assert!(summary.contains("| Positions | 2 |"));
        assert!(summary.contains("| Mean Quality | 25.00 |"));
        assert!(summary.contains("| Lowest Mean Quality | 10.00 (pos 0) |"));
        assert!(summary.contains("| Highest Mean Quality | 40.00 (pos 1) |"));
    }

    #[test]
    fn ingest_merges_counts_and_zeroes_source() {
        let mut a = BaseHistogram::default();
        let mut b = BaseHistogram::default();
        a.push(&phred(&[10, 10]));
        b.push(&phred(&[30, 30]));

        a.ingest(&mut b);

        assert_eq!(a.position_means(), vec![20.0, 20.0]);
        assert_eq!(b.position_means(), vec![0.0, 0.0]);
    }
}
