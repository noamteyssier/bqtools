use std::{
    io::Write,
    path::Path,
    sync::{Arc, Mutex},
};

use anyhow::Result;
use serde::Serialize;

use crate::commands::match_output;

/// `n` as a percentage of `total` (0 when `total` is 0).
pub fn pct(n: usize, total: usize) -> f64 {
    if total == 0 {
        0.0
    } else {
        (n as f64 / total as f64) * 100.0
    }
}

/// Writes `rows` as a tab-separated file with a header row.
pub fn write_tsv<W: Write, S: Serialize>(wtr: W, rows: impl IntoIterator<Item = S>) -> Result<()> {
    let mut ser = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .from_writer(wtr);
    for row in rows {
        ser.serialize(row)?;
    }
    ser.flush()?;
    Ok(())
}

/// Adds `src` into `dst` element-wise and zeroes `src`.
pub fn add_assign(dst: &mut [usize], src: &mut [usize]) {
    for (d, s) in dst.iter_mut().zip(src) {
        *d += std::mem::take(s);
    }
}

/// `(total, mean, median, mode)` of a histogram indexed by value.
pub fn stats(counts: &[usize]) -> (usize, f64, usize, usize) {
    let total: usize = counts.iter().sum();
    if total == 0 {
        return (0, 0.0, 0, counts.len().saturating_sub(1));
    }
    let sum: usize = counts.iter().enumerate().map(|(v, &c)| v * c).sum();
    let half = total / 2;
    let mut cum = 0;
    let median = counts
        .iter()
        .position(|&c| {
            cum += c;
            cum > half
        })
        .unwrap_or(0);
    let mode = counts
        .iter()
        .enumerate()
        .max_by_key(|&(_, &c)| c)
        .map_or(0, |(v, _)| v);
    (total, sum as f64 / total as f64, median, mode)
}

/// A per-side (R1 or R2) accumulator merged across threads.
pub trait Hist: Default {
    fn is_empty(&self) -> bool;
    /// Merges `other` into `self`, zeroing `other`.
    fn ingest(&mut self, other: &mut Self);
    fn serialize_to<W: Write>(&self, wtr: &mut W) -> Result<()>;
    fn summary_table(&self) -> Option<String>;
}

/// Thread-local (`t`) and shared (`g`) accumulators for the primary (R1) and
/// extended (R2) reads. Cloning shares `g` and starts with fresh `t`.
#[derive(Default)]
pub struct Pair<H> {
    pub t: [H; 2],
    g: Arc<Mutex<[H; 2]>>,
}
impl<H: Hist> Clone for Pair<H> {
    fn clone(&self) -> Self {
        Self {
            t: Default::default(),
            g: Arc::clone(&self.g),
        }
    }
}
impl<H: Hist> Pair<H> {
    /// Merges this thread's accumulators into the shared ones.
    pub fn sync_final(&mut self) {
        let mut g = self.g.lock().unwrap();
        for (g, t) in g.iter_mut().zip(&mut self.t) {
            g.ingest(t);
        }
    }

    /// Runs `f` on each non-empty side with its `R1`/`R2` label.
    pub fn each(&self, mut f: impl FnMut(&H, &str) -> Result<()>) -> Result<()> {
        let g = self.g.lock().unwrap();
        for (h, side) in g.iter().zip(["R1", "R2"]) {
            if !h.is_empty() {
                f(h, side)?;
            }
        }
        Ok(())
    }

    /// Writes `{stem}_R1.tsv` and `{stem}_R2.tsv` into `outdir` for each non-empty side.
    pub fn write(&self, outdir: &Path, stem: &str) -> Result<()> {
        self.each(|h, side| {
            let mut handle = match_output(Some(outdir.join(format!("{stem}_{side}.tsv"))))?;
            h.serialize_to(&mut handle)
        })
    }

    /// Maps `f` over the shared primary and extended accumulators.
    pub fn map<T>(&self, f: impl Fn(&H) -> T) -> (T, T) {
        let g = self.g.lock().unwrap();
        (f(&g[0]), f(&g[1]))
    }

    pub fn summarize(&self, title: &str) -> String {
        let (primary, extended) = self.map(H::summary_table);
        dual_section(title, primary, extended)
    }
}

/// Renders a markdown table. Returns an empty string if `rows` is empty, so
/// callers can unconditionally splice the result into a report without an
/// extra emptiness check.
pub fn table(headers: &[&str], rows: &[Vec<String>]) -> String {
    use std::fmt::Write as _;

    if rows.is_empty() {
        return String::new();
    }

    let mut out = format!("| {} |\n", headers.join(" | "));
    let _ = writeln!(
        out,
        "|{}|",
        headers.iter().map(|_| "---").collect::<Vec<_>>().join("|")
    );
    for row in rows {
        let _ = writeln!(out, "| {} |", row.join(" | "));
    }
    out
}

/// Wraps a module's already-rendered primary/extended (R1/R2) summary
/// bodies under one heading. `None` means that side had no data to report
/// (e.g. `extended` is always `None` for single-end input). Returns an empty
/// string if both sides are `None`, so callers can splice the result in
/// unconditionally.
pub fn dual_section(title: &str, primary: Option<String>, extended: Option<String>) -> String {
    if primary.is_none() && extended.is_none() {
        return String::new();
    }

    let split = primary.is_some() && extended.is_some();
    let mut out = format!("## {title}\n\n");
    if let Some(body) = primary {
        if split {
            out.push_str("### R1\n\n");
        }
        out.push_str(&body);
        out.push('\n');
    }
    if let Some(body) = extended {
        if split {
            out.push_str("### R2\n\n");
        }
        out.push_str(&body);
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_with_no_rows_is_empty() {
        assert_eq!(table(&["A", "B"], &[]), "");
    }

    #[test]
    fn table_renders_header_separator_and_rows() {
        let rows = vec![
            vec!["1".to_string(), "2".to_string()],
            vec!["3".to_string(), "4".to_string()],
        ];
        assert_eq!(
            table(&["A", "B"], &rows),
            "| A | B |\n|---|---|\n| 1 | 2 |\n| 3 | 4 |\n"
        );
    }

    #[test]
    fn dual_section_with_neither_side_is_empty() {
        assert_eq!(dual_section("Title", None, None), "");
    }

    #[test]
    fn dual_section_primary_only_has_no_r1_r2_headings() {
        let out = dual_section("Title", Some("body\n".to_string()), None);
        assert!(out.contains("## Title"));
        assert!(out.contains("body"));
        assert!(!out.contains("### R1"));
        assert!(!out.contains("### R2"));
    }

    #[test]
    fn dual_section_extended_only_has_no_r1_r2_headings() {
        let out = dual_section("Title", None, Some("body\n".to_string()));
        assert!(out.contains("## Title"));
        assert!(out.contains("body"));
        assert!(!out.contains("### R1"));
        assert!(!out.contains("### R2"));
    }

    #[test]
    fn dual_section_both_sides_split_r1_before_r2() {
        let out = dual_section(
            "Title",
            Some("primary body\n".to_string()),
            Some("extended body\n".to_string()),
        );
        assert!(out.contains("primary body"));
        assert!(out.contains("extended body"));
        assert!(out.find("### R1").unwrap() < out.find("### R2").unwrap());
    }
}
