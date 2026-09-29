use std::{io::Write, path::Path};

use anyhow::Result;
use binseq::BinseqRecord;
use serde::Serialize;

use super::report::{add_assign, pct, table, write_tsv, Hist, Pair};

const NUM_BASES: usize = 5;
const IDX_A: usize = 0;
const IDX_C: usize = 1;
const IDX_G: usize = 2;
const IDX_T: usize = 3;
const IDX_N: usize = 4;

pub type BaseAbundance = [usize; NUM_BASES];
pub const DEFAULT_BASE_ABUNDANCE: BaseAbundance = [0; NUM_BASES];

/// Byte -> histogram index, built once at compile time.
///
/// A `match` on `A/C/G/T` (case-insensitive) compiles to a chain of
/// unpredictable branches, since those bytes aren't a contiguous range. A
/// 256-entry table turns the per-base lookup into one branchless load.
const BASE_LUT: [u8; 256] = {
    let mut lut = [IDX_N as u8; 256];
    lut[b'A' as usize] = IDX_A as u8;
    lut[b'a' as usize] = IDX_A as u8;
    lut[b'C' as usize] = IDX_C as u8;
    lut[b'c' as usize] = IDX_C as u8;
    lut[b'G' as usize] = IDX_G as u8;
    lut[b'g' as usize] = IDX_G as u8;
    lut[b'T' as usize] = IDX_T as u8;
    lut[b't' as usize] = IDX_T as u8;
    lut
};

/// Buckets a decoded base into its histogram index.
///
/// Anything outside `ACGT` (case-insensitive) - ambiguity codes included -
/// is folded into the `N` bucket.
#[inline]
fn base_index(base: u8) -> usize {
    BASE_LUT[base as usize] as usize
}

#[derive(Serialize)]
pub struct BaseContentRecord {
    pos: usize,
    a: usize,
    c: usize,
    g: usize,
    t: usize,
    n: usize,
    pct_a: f64,
    pct_c: f64,
    pct_g: f64,
    pct_t: f64,
    pct_n: f64,
}

#[derive(Debug, Clone, Default)]
pub struct BaseContentHistogram {
    /// Outer: position
    /// Inner: base abundance (A, C, G, T, N)
    inner: Vec<BaseAbundance>,
}
impl BaseContentHistogram {
    /// Number of positions tracked
    fn len(&self) -> usize {
        self.inner.len()
    }
    /// Track base identity histogram over positions
    fn push(&mut self, seq: &[u8]) {
        if seq.is_empty() {
            return;
        }
        if self.inner.len() < seq.len() {
            self.inner.resize(seq.len(), DEFAULT_BASE_ABUNDANCE);
        }
        seq.iter()
            .zip(self.inner.iter_mut())
            .for_each(|(&base, pos_vec)| {
                pos_vec[base_index(base)] += 1;
            });
    }

    /// Aggregate base composition across all positions.
    fn totals(&self) -> BaseAbundance {
        let mut totals = DEFAULT_BASE_ABUNDANCE;
        for counts in &self.inner {
            for (t, &c) in totals.iter_mut().zip(counts.iter()) {
                *t += c;
            }
        }
        totals
    }
}
impl Hist for BaseContentHistogram {
    /// Checks if empty
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    fn ingest(&mut self, other: &mut Self) {
        if self.len() < other.len() {
            self.inner.resize(other.len(), DEFAULT_BASE_ABUNDANCE);
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
            self.inner.iter().enumerate().map(|(pos, counts)| {
                let total: usize = counts.iter().sum();
                let pct = |c: usize| pct(c, total);
                BaseContentRecord {
                    pos,
                    a: counts[IDX_A],
                    c: counts[IDX_C],
                    g: counts[IDX_G],
                    t: counts[IDX_T],
                    n: counts[IDX_N],
                    pct_a: pct(counts[IDX_A]),
                    pct_c: pct(counts[IDX_C]),
                    pct_g: pct(counts[IDX_G]),
                    pct_t: pct(counts[IDX_T]),
                    pct_n: pct(counts[IDX_N]),
                }
            }),
        )
    }

    fn summary_table(&self) -> Option<String> {
        if self.is_empty() {
            return None;
        }

        let totals = self.totals();
        let total: usize = totals.iter().sum();

        Some(table(
            &["Base", "Count", "Pct"],
            &[
                vec![
                    "A".into(),
                    totals[IDX_A].to_string(),
                    format!("{:.2}%", pct(totals[IDX_A], total)),
                ],
                vec![
                    "C".into(),
                    totals[IDX_C].to_string(),
                    format!("{:.2}%", pct(totals[IDX_C], total)),
                ],
                vec![
                    "G".into(),
                    totals[IDX_G].to_string(),
                    format!("{:.2}%", pct(totals[IDX_G], total)),
                ],
                vec![
                    "T".into(),
                    totals[IDX_T].to_string(),
                    format!("{:.2}%", pct(totals[IDX_T], total)),
                ],
                vec![
                    "N".into(),
                    totals[IDX_N].to_string(),
                    format!("{:.2}%", pct(totals[IDX_N], total)),
                ],
            ],
        ))
    }
}

#[derive(Clone, Default)]
pub struct PerBaseSequenceContent(Pair<BaseContentHistogram>);
impl PerBaseSequenceContent {
    pub fn push<R: BinseqRecord>(&mut self, record: &R) {
        self.0.t[0].push(record.sseq());
        self.0.t[1].push(record.xseq());
    }

    pub fn sync_final(&mut self) {
        self.0.sync_final();
    }

    pub fn finish(&mut self, outdir: &Path) -> Result<()> {
        if !outdir.exists() {
            std::fs::create_dir_all(outdir)?;
        }
        self.0.write(outdir, "base_content")
    }

    pub fn summarize(&self) -> String {
        self.0.summarize("Per-Base Sequence Content")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn starts_empty() {
        assert!(BaseContentHistogram::default().is_empty());
    }

    #[test]
    fn push_ignores_empty_sequence() {
        let mut hist = BaseContentHistogram::default();
        hist.push(b"");
        assert!(hist.is_empty());
    }

    #[test]
    fn push_tracks_base_identity_per_position() {
        let mut hist = BaseContentHistogram::default();
        hist.push(b"ACGTN");
        assert!(!hist.is_empty());
        assert_eq!(hist.len(), 5);
    }

    #[test]
    fn push_folds_ambiguity_codes_and_lowercase_into_expected_buckets() {
        let mut hist = BaseContentHistogram::default();
        hist.push(b"ACGTRacgt");
        let totals = hist.totals();
        assert_eq!(totals[IDX_A], 2); // 'A' and 'a'
        assert_eq!(totals[IDX_C], 2); // 'C' and 'c'
        assert_eq!(totals[IDX_G], 2); // 'G' and 'g'
        assert_eq!(totals[IDX_T], 2); // 'T' and 't'
        assert_eq!(totals[IDX_N], 1); // 'R' (ambiguity code)
    }

    #[test]
    fn summary_table_none_when_empty() {
        assert!(BaseContentHistogram::default().summary_table().is_none());
    }

    #[test]
    fn summary_table_reports_base_composition() {
        let mut hist = BaseContentHistogram::default();
        hist.push(b"AACG");
        let summary = hist.summary_table().expect("non-empty histogram");
        assert!(summary.contains("| A | 2 | 50.00% |"));
        assert!(summary.contains("| C | 1 | 25.00% |"));
        assert!(summary.contains("| G | 1 | 25.00% |"));
        assert!(summary.contains("| T | 0 | 0.00% |"));
        assert!(summary.contains("| N | 0 | 0.00% |"));
    }

    #[test]
    fn ingest_merges_counts() {
        let mut a = BaseContentHistogram::default();
        let mut b = BaseContentHistogram::default();
        a.push(b"AA");
        b.push(b"CC");

        a.ingest(&mut b);

        let totals = a.totals();
        assert_eq!(totals[IDX_A], 2);
        assert_eq!(totals[IDX_C], 2);
    }
}
