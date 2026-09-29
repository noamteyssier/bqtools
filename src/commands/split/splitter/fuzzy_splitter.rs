use anyhow::Result;
use fixedbitset::FixedBitSet;
use sassy::{profiles::Iupac, EncodedPatterns, Searcher};

use crate::commands::{grep::PatternCollection, utils::build_fuzzy_searcher};

type Profile = Iupac;

/// A searcher with its encoded patterns; `None` when the pattern set is empty.
type Set = Option<(Searcher<Profile>, EncodedPatterns<Profile>)>;

/// Fuzzy (edit-distance) matching via `sassy` over the three pattern sets
/// (primary-only, secondary-only, either).
#[derive(Clone)]
pub struct FuzzySplitter {
    sets: [Set; 3],

    /// Global pattern index of the first pattern in each set
    offsets: [usize; 3],

    /// Maximum edit distance to accept
    k: usize,
    /// Whether to only accept inexact matches
    inexact: bool,
}

impl FuzzySplitter {
    pub fn new(
        pat1: &PatternCollection,
        pat2: &PatternCollection,
        pat: &PatternCollection,
        k: usize,
        inexact: bool,
        max_n_frac: Option<f32>,
    ) -> Result<Self> {
        // validate lengths, resolve max_n_frac, and encode patterns per pattern set
        let build = |p: &PatternCollection| -> Result<Set> {
            let (searcher, enc) = build_fuzzy_searcher(&p.bytes(), k, max_n_frac)?;
            Ok(enc.map(|e| (searcher, e)))
        };
        Ok(Self {
            sets: [build(pat1)?, build(pat2)?, build(pat)?],
            offsets: [0, pat1.len(), pat1.len() + pat2.len()],
            k,
            inexact,
        })
    }

    /// Sets a bit in `bits` for every pattern found in the sequences.
    pub fn hit(&mut self, primary: &[u8], secondary: &[u8], bits: &mut FixedBitSet) {
        let seqs: [&[&[u8]]; 3] = [&[primary], &[secondary], &[primary, secondary]];
        for ((set, offset), seqs) in self.sets.iter_mut().zip(self.offsets).zip(seqs) {
            let Some((searcher, patterns)) = set else {
                continue;
            };
            for seq in seqs.iter().filter(|s| !s.is_empty()) {
                for m in searcher.search_encoded_patterns(patterns, seq, self.k) {
                    if !(self.inexact && m.cost == 0) {
                        bits.insert(offset + m.pattern_idx);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::FuzzySplitter;
    use crate::commands::grep::Pattern;
    use crate::commands::grep::PatternCollection;
    use crate::commands::split::splitter::{Matcher, Splitter};

    fn pc(patterns: &[&[u8]], name: &str) -> PatternCollection {
        PatternCollection(
            patterns
                .iter()
                .map(|p| Pattern {
                    name: Some(name.to_string()),
                    sequence: p.to_vec(),
                })
                .collect(),
        )
    }

    // Mirrors https://github.com/RagnarGrootKoerkamp/sassy/issues/66: the Iupac
    // profile treats `N` as a wildcard, so without an N-fraction filter a needle
    // matches a haystack made entirely of `N`s.
    #[test]
    fn test_fuzzy_splitter_default_max_n_frac_rejects_all_n_match() {
        let pat1 = pc(&[b"ACGTACGTACGT"], "alias");
        let empty = PatternCollection(vec![]);
        let m = FuzzySplitter::new(&pat1, &empty, &empty, 1, false, None).unwrap();
        let mut splitter = Splitter::new(Matcher::Fuzzy(Box::new(m)), &pat1, &empty, &empty);

        let all_n = b"NNNNNNNNNNNNNNNNNN";
        assert_eq!(
            splitter.split_idx(all_n, b""),
            None,
            "default max_n_frac (k/pattern_len) should reject an all-N match"
        );
    }

    #[test]
    fn test_fuzzy_splitter_max_n_frac_override_allows_all_n_match() {
        let pat1 = pc(&[b"ACGTACGTACGT"], "alias");
        let empty = PatternCollection(vec![]);
        let m = FuzzySplitter::new(&pat1, &empty, &empty, 1, false, Some(1.0)).unwrap();
        let mut splitter = Splitter::new(Matcher::Fuzzy(Box::new(m)), &pat1, &empty, &empty);

        let all_n = b"NNNNNNNNNNNNNNNNNN";
        assert_eq!(
            splitter.split_idx(all_n, b""),
            Some(0),
            "max_n_frac=1.0 should disable the N-fraction filter"
        );
    }

    // sassy's `Searcher::encode_patterns` panics (`assert!`) when a pattern set
    // contains mixed lengths; these tests confirm we catch that up front and
    // return an `Err` instead of letting the panic reach the caller.
    #[test]
    fn test_fuzzy_splitter_rejects_mismatched_pattern_lengths_primary() {
        let pat1 = pc(&[b"AAAA", b"AAAAA"], "alias");
        let empty = PatternCollection(vec![]);
        let result = FuzzySplitter::new(&pat1, &empty, &empty, 1, false, None);
        assert!(
            result.is_err(),
            "mismatched primary pattern lengths should error, not panic"
        );
    }

    #[test]
    fn test_fuzzy_splitter_rejects_mismatched_pattern_lengths_secondary() {
        let pat2 = pc(&[b"AAAA", b"AAAAA"], "alias");
        let empty = PatternCollection(vec![]);
        let result = FuzzySplitter::new(&empty, &pat2, &empty, 1, false, None);
        assert!(
            result.is_err(),
            "mismatched secondary pattern lengths should error, not panic"
        );
    }

    #[test]
    fn test_fuzzy_splitter_rejects_mismatched_pattern_lengths_either() {
        let pat = pc(&[b"AAAA", b"AAAAA"], "alias");
        let empty = PatternCollection(vec![]);
        let result = FuzzySplitter::new(&empty, &empty, &pat, 1, false, None);
        assert!(
            result.is_err(),
            "mismatched either-set pattern lengths should error, not panic"
        );
    }

    #[test]
    fn test_fuzzy_splitter_accepts_uniform_pattern_lengths() {
        let pat1 = pc(&[b"AAAA", b"TTTT", b"CCCC"], "alias");
        let empty = PatternCollection(vec![]);
        let result = FuzzySplitter::new(&pat1, &empty, &empty, 1, false, None);
        assert!(result.is_ok(), "uniform pattern lengths should not error");
    }
}
