//! Pattern matching shared by grep, `grep -P` and split.
//!
//! Patterns come in three sets: primary-only, extended-only and either. They
//! are numbered consecutively in that order, and [`Engine::hit`] sets the bit
//! of every pattern found in a record.

use aho_corasick::AhoCorasick;
use anyhow::Result;
use fixedbitset::FixedBitSet;
#[cfg(feature = "fuzzy")]
use sassy::{profiles::Iupac, EncodedPatterns, Searcher};

use super::PatternSets;
#[cfg(feature = "fuzzy")]
use crate::commands::utils::build_fuzzy_searcher;
use crate::commands::utils::corasick_builder;

/// Half-open `(start, end)` of a hit within the sequence that was searched.
pub type Span = (usize, usize);

/// Where patterns hit in the primary and extended sequences. May hold
/// duplicates and overlaps.
#[derive(Clone, Debug, Default)]
pub struct Spans {
    pub primary: Vec<Span>,
    pub secondary: Vec<Span>,
}
impl Spans {
    pub fn clear(&mut self) {
        self.primary.clear();
        self.secondary.clear();
    }

    /// Moves every span right by `offset` (e.g. to undo a sliced range).
    pub fn shift(&mut self, offset: usize) {
        if offset == 0 {
            return;
        }
        for (start, end) in self.primary.iter_mut().chain(&mut self.secondary) {
            *start += offset;
            *end += offset;
        }
    }
}

/// One pattern set compiled for a matching strategy.
#[derive(Clone)]
enum Set {
    /// `None` when the set is empty.
    AhoCorasick(Option<AhoCorasick>),
    Regex(Vec<regex::bytes::Regex>),
    #[cfg(feature = "fuzzy")]
    Fuzzy(Box<FuzzySet>),
}

#[cfg(feature = "fuzzy")]
#[derive(Clone)]
struct FuzzySet {
    /// `None` when the set is empty.
    inner: Option<(Searcher<Iupac>, EncodedPatterns<Iupac>)>,
    /// Maximum edit distance to accept
    k: usize,
    /// Whether to only accept inexact (nonzero cost) matches
    inexact: bool,
}

impl Set {
    /// Sets `bits[offset + i]` for every pattern `i` found in `seq`, and
    /// records where when `spans` is given (regex only scans for positions then).
    fn scan(
        &mut self,
        seq: &[u8],
        offset: usize,
        bits: &mut FixedBitSet,
        mut spans: Option<&mut Vec<Span>>,
    ) {
        match self {
            Set::AhoCorasick(ac) => {
                let Some(ac) = ac else { return };
                for m in ac.find_overlapping_iter(seq) {
                    bits.insert(offset + m.pattern().as_usize());
                    if let Some(spans) = spans.as_deref_mut() {
                        spans.push((m.start(), m.end()));
                    }
                }
            }
            Set::Regex(regexes) => {
                for (idx, re) in regexes.iter().enumerate() {
                    if let Some(spans) = spans.as_deref_mut() {
                        for m in re.find_iter(seq) {
                            bits.insert(offset + idx);
                            spans.push((m.start(), m.end()));
                        }
                    } else if re.is_match(seq) {
                        bits.insert(offset + idx);
                    }
                }
            }
            #[cfg(feature = "fuzzy")]
            Set::Fuzzy(set) => {
                let (k, inexact) = (set.k, set.inexact);
                let Some((searcher, patterns)) = &mut set.inner else {
                    return;
                };
                for m in searcher.search_encoded_patterns(patterns, seq, k) {
                    if inexact && m.cost == 0 {
                        continue;
                    }
                    bits.insert(offset + m.pattern_idx);
                    if let Some(spans) = spans.as_deref_mut() {
                        spans.push((m.text_start, m.text_end));
                    }
                }
            }
        }
    }
}

impl Set {
    /// Whether every pattern of the set (`range` of the global numbering) is
    /// found in at least one of `seqs`. Regex stops at the first missing pattern.
    fn scan_all(
        &mut self,
        seqs: &[&[u8]],
        range: std::ops::Range<usize>,
        bits: &mut FixedBitSet,
    ) -> bool {
        if let Set::Regex(regexes) = self {
            return regexes
                .iter()
                .all(|re| seqs.iter().any(|s| !s.is_empty() && re.is_match(s)));
        }
        for seq in seqs.iter().filter(|s| !s.is_empty()) {
            self.scan(seq, range.start, bits, None);
        }
        range.into_iter().all(|idx| bits.contains(idx))
    }

    /// Whether any pattern is found in `seq`; stops at the first hit.
    fn scan_any(&mut self, seq: &[u8]) -> bool {
        match self {
            Set::AhoCorasick(ac) => ac.as_ref().is_some_and(|ac| ac.is_match(seq)),
            Set::Regex(regexes) => regexes.iter().any(|re| re.is_match(seq)),
            #[cfg(feature = "fuzzy")]
            Set::Fuzzy(set) => {
                let (k, inexact) = (set.k, set.inexact);
                let Some((searcher, patterns)) = &mut set.inner else {
                    return false;
                };
                searcher
                    .search_encoded_patterns(patterns, seq, k)
                    .iter()
                    .any(|m| !(inexact && m.cost == 0))
            }
        }
    }
}

/// Finds which of the three pattern sets' patterns occur in a record.
#[derive(Clone)]
pub struct Engine {
    sets: [Set; 3],
    /// Global pattern index of the first pattern in each set
    offsets: [usize; 3],
    /// Total number of patterns
    len: usize,
}
impl Engine {
    fn build(
        patterns: &PatternSets,
        mut compile: impl FnMut(&super::PatternCollection) -> Result<Set>,
    ) -> Result<Self> {
        let [p1, p2, p] = patterns.each();
        Ok(Self {
            sets: [compile(p1)?, compile(p2)?, compile(p)?],
            offsets: [0, p1.len(), p1.len() + p2.len()],
            len: patterns.len(),
        })
    }

    /// Fixed-string matching.
    pub fn aho_corasick(patterns: &PatternSets, no_dfa: bool) -> Result<Self> {
        Self::build(patterns, |p| {
            let ac = (!p.is_empty())
                .then(|| corasick_builder(&p.bytes(), no_dfa))
                .transpose()?;
            Ok(Set::AhoCorasick(ac))
        })
    }

    /// Regular-expression matching.
    pub fn regex(patterns: &PatternSets) -> Result<Self> {
        Self::build(patterns, |p| Ok(Set::Regex(p.regexes()?)))
    }

    /// Edit-distance matching via `sassy`. Patterns within a set must share a length.
    #[cfg(feature = "fuzzy")]
    pub fn fuzzy(
        patterns: &PatternSets,
        k: usize,
        inexact: bool,
        max_n_frac: Option<f32>,
    ) -> Result<Self> {
        Self::build(patterns, |p| {
            let (searcher, encoded) = build_fuzzy_searcher(&p.bytes(), k, max_n_frac)?;
            Ok(Set::Fuzzy(Box::new(FuzzySet {
                inner: encoded.map(|e| (searcher, e)),
                k,
                inexact,
            })))
        })
    }

    /// Whether any pattern occurs in the record, searched like [`Self::hit`]
    /// but stopping at the first hit.
    pub fn any(&mut self, primary: &[u8], secondary: &[u8]) -> bool {
        let [s1, s2, s] = &mut self.sets;
        (!primary.is_empty() && (s1.scan_any(primary) || s.scan_any(primary)))
            || (!secondary.is_empty() && (s2.scan_any(secondary) || s.scan_any(secondary)))
    }

    /// Whether every pattern occurs in the record (an either-pattern in either
    /// sequence), searched like [`Self::hit`] but giving up at the first
    /// missing pattern where the backend allows. `bits` is scratch space from
    /// [`Self::bitset`].
    pub fn all(&mut self, primary: &[u8], secondary: &[u8], bits: &mut FixedBitSet) -> bool {
        bits.clear();
        let [s1, s2, s] = &mut self.sets;
        let [o1, o2, o] = self.offsets;
        s1.scan_all(&[primary], o1..o2, bits)
            && s2.scan_all(&[secondary], o2..o, bits)
            && s.scan_all(&[primary, secondary], o..self.len, bits)
    }

    /// An empty bitset with one bit per pattern, for [`Self::hit`].
    pub fn bitset(&self) -> FixedBitSet {
        FixedBitSet::with_capacity(self.len)
    }

    /// Sets a bit in `bits` for every pattern found: primary-only patterns are
    /// searched in `primary`, extended-only in `secondary`, and either-patterns
    /// in both. Empty sequences are skipped. When `spans` is given, hit
    /// positions (relative to each sequence) are appended to it.
    pub fn hit(
        &mut self,
        primary: &[u8],
        secondary: &[u8],
        bits: &mut FixedBitSet,
        spans: Option<&mut Spans>,
    ) {
        let (mut sp, mut sx) = match spans {
            Some(s) => (Some(&mut s.primary), Some(&mut s.secondary)),
            None => (None, None),
        };
        let [s1, s2, s] = &mut self.sets;
        let [o1, o2, o] = self.offsets;
        if !primary.is_empty() {
            s1.scan(primary, o1, bits, sp.as_deref_mut());
            s.scan(primary, o, bits, sp);
        }
        if !secondary.is_empty() {
            s2.scan(secondary, o2, bits, sx.as_deref_mut());
            s.scan(secondary, o, bits, sx);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Engine, Spans};
    use crate::commands::grep::{Pattern, PatternCollection, PatternSets};

    fn pc(patterns: &[&str]) -> PatternCollection {
        PatternCollection(
            patterns
                .iter()
                .map(|p| Pattern {
                    name: None,
                    sequence: p.as_bytes().to_vec(),
                })
                .collect(),
        )
    }

    fn sets(pat1: &[&str], pat2: &[&str], pat: &[&str]) -> PatternSets {
        PatternSets {
            pat1: pc(pat1),
            pat2: pc(pat2),
            pat: pc(pat),
        }
    }

    /// Every backend, exact matching. Fuzzy needs uniform lengths per set.
    fn engines(sets: &PatternSets) -> Vec<(&'static str, Engine)> {
        vec![
            ("aho-corasick", Engine::aho_corasick(sets, false).unwrap()),
            (
                "aho-corasick (no dfa)",
                Engine::aho_corasick(sets, true).unwrap(),
            ),
            ("regex", Engine::regex(sets).unwrap()),
            #[cfg(feature = "fuzzy")]
            ("fuzzy", Engine::fuzzy(sets, 0, false, None).unwrap()),
        ]
    }

    /// Indices of the patterns found, and their spans.
    fn hit(engine: &mut Engine, primary: &str, secondary: &str) -> (Vec<usize>, Spans) {
        let mut bits = engine.bitset();
        let mut spans = Spans::default();
        engine.hit(
            primary.as_bytes(),
            secondary.as_bytes(),
            &mut bits,
            Some(&mut spans),
        );
        (bits.ones().collect(), spans)
    }

    #[test]
    fn test_pattern_sets_are_numbered_and_scoped() {
        // 0: primary-only, 1: extended-only, 2: either
        let sets = sets(&["AAAA"], &["CCCC"], &["GGGG"]);
        for (name, mut e) in engines(&sets) {
            assert_eq!(e.bitset().len(), 3, "{name}");
            // a primary-only pattern is not searched in the extended sequence, and vice versa
            assert_eq!(hit(&mut e, "TTTT", "AAAA").0, [] as [usize; 0], "{name}");
            assert_eq!(hit(&mut e, "CCCC", "TTTT").0, [] as [usize; 0], "{name}");
            assert_eq!(hit(&mut e, "AAAA", "CCCC").0, [0, 1], "{name}");
            // an either-pattern hits in whichever sequence has it
            assert_eq!(hit(&mut e, "GGGG", "TTTT").0, [2], "{name}");
            assert_eq!(hit(&mut e, "TTTT", "GGGG").0, [2], "{name}");
            assert_eq!(hit(&mut e, "AAAAGGGG", "CCCC").0, [0, 1, 2], "{name}");
        }
    }

    #[test]
    fn test_any_is_scoped_like_hit() {
        let sets = sets(&["AAAA"], &["CCCC"], &["GGGG"]);
        for (name, mut e) in engines(&sets) {
            assert!(!e.any(b"", b""), "{name}");
            assert!(!e.any(b"TTTT", b"TTTT"), "{name}");
            // a primary-only pattern is not searched in the extended sequence, and vice versa
            assert!(!e.any(b"TTTT", b"AAAA"), "{name}");
            assert!(!e.any(b"CCCC", b"TTTT"), "{name}");
            assert!(e.any(b"AAAA", b""), "{name}");
            assert!(e.any(b"", b"CCCC"), "{name}");
            assert!(e.any(b"TTTT", b"GGGG"), "{name}");
            assert!(e.any(b"GGGG", b""), "{name}");
        }
    }

    #[test]
    fn test_all_requires_every_pattern() {
        // either-patterns may be found in different mates
        let sets = sets(&["AAAA"], &["CCCC"], &["GGGG", "TTTT"]);
        for (name, mut e) in engines(&sets) {
            let mut bits = e.bitset();
            assert!(e.all(b"AAAAGGGG", b"CCCCTTTT", &mut bits), "{name}");
            assert!(e.all(b"AAAAGGGGTTTT", b"CCCC", &mut bits), "{name}");
            assert!(
                !e.all(b"AAAAGGGG", b"CCCC", &mut bits),
                "{name}: TTTT is missing"
            );
            assert!(
                !e.all(b"AAAAGGGGTTTT", b"", &mut bits),
                "{name}: no extended sequence"
            );
            assert!(
                !e.all(b"", b"CCCCGGGGTTTT", &mut bits),
                "{name}: AAAA is primary-only"
            );
            assert!(
                !e.all(b"CCCC", b"AAAA", &mut bits),
                "{name}: sets are scoped"
            );
        }
    }

    #[test]
    fn test_all_of_no_patterns_holds() {
        for (name, mut e) in engines(&sets(&[], &[], &[])) {
            let mut bits = e.bitset();
            assert!(e.all(b"AAAA", b"", &mut bits), "{name}");
        }
    }

    #[test]
    fn test_empty_sequences_are_skipped() {
        let sets = sets(&["AAAA"], &["CCCC"], &["GGGG"]);
        for (name, mut e) in engines(&sets) {
            assert_eq!(hit(&mut e, "", "").0, [] as [usize; 0], "{name}");
            assert_eq!(hit(&mut e, "AAAAGGGG", "").0, [0, 2], "{name}");
            assert_eq!(hit(&mut e, "", "CCCCGGGG").0, [1, 2], "{name}");
        }
    }

    #[test]
    fn test_empty_pattern_sets() {
        let sets = sets(&[], &[], &[]);
        for (name, mut e) in engines(&sets) {
            assert_eq!(e.bitset().len(), 0, "{name}");
            assert_eq!(hit(&mut e, "AAAA", "AAAA").0, [] as [usize; 0], "{name}");
        }
    }

    #[test]
    fn test_spans_are_relative_to_each_sequence() {
        let sets = sets(&["AAAA"], &["CCCC"], &["GGGG"]);
        for (name, mut e) in engines(&sets) {
            let (_, spans) = hit(&mut e, "TAAAAGGGG", "CCCCTGGGG");
            let (mut primary, mut secondary) = (spans.primary, spans.secondary);
            primary.sort_unstable();
            secondary.sort_unstable();
            assert_eq!(primary, [(1, 5), (5, 9)], "{name}");
            assert_eq!(secondary, [(0, 4), (5, 9)], "{name}");
        }
    }

    #[test]
    fn test_no_spans_requested_still_sets_bits() {
        let sets = sets(&["AAAA"], &[], &["GGGG"]);
        for (name, mut e) in engines(&sets) {
            let mut bits = e.bitset();
            e.hit(b"AAAAGGGG", b"", &mut bits, None);
            assert_eq!(bits.ones().collect::<Vec<_>>(), [0, 1], "{name}");
        }
    }

    #[test]
    fn test_every_pattern_reports_spans() {
        // both patterns hit, and both are located (regex does not stop at the first)
        let sets = sets(&[], &[], &["GATTA", "TTACA"]);
        for (name, mut e) in engines(&sets) {
            let (bits, spans) = hit(&mut e, "GATTACA", "");
            let mut primary = spans.primary;
            primary.sort_unstable();
            assert_eq!(bits, [0, 1], "{name}");
            assert_eq!(primary, [(0, 5), (2, 7)], "{name}");
        }
    }

    #[test]
    fn test_repeated_hits_cover_the_occurrences() {
        // regex reports non-overlapping hits and aho-corasick overlapping ones
        // (the colour writer merges them); sassy may report even fewer
        let sets = sets(&["AA"], &[], &[]);
        for (name, mut e) in engines(&sets) {
            let (bits, spans) = hit(&mut e, "AAAATT", "");
            let mut covered = [false; 6];
            spans
                .primary
                .iter()
                .for_each(|&(start, end)| covered[start..end].fill(true));
            assert_eq!(bits, [0], "{name}");
            assert_eq!(covered[2..], [true, true, false, false], "{name}");
        }
    }

    #[test]
    fn test_regex_patterns() {
        let sets = sets(&["A.GT", "^ACGT$"], &[], &[]);
        let mut e = Engine::regex(&sets).unwrap();
        assert_eq!(hit(&mut e, "TTACGTTT", "").0, [0]);
        assert_eq!(hit(&mut e, "ACGT", "").0, [0, 1]);
    }

    #[test]
    fn test_aho_corasick_is_case_sensitive() {
        let sets = sets(&["AAAA"], &[], &[]);
        let mut e = Engine::aho_corasick(&sets, false).unwrap();
        assert_eq!(hit(&mut e, "aaaa", "").0, [] as [usize; 0]);
        assert_eq!(hit(&mut e, "AAAA", "").0, [0]);
    }

    #[test]
    fn test_invalid_regex_is_an_error() {
        assert!(Engine::regex(&sets(&["A("], &[], &[])).is_err());
    }

    #[cfg(feature = "fuzzy")]
    mod fuzzy {
        use super::{hit, sets, Engine};

        #[test]
        fn test_edit_distance() {
            let sets = sets(&["ACGTACGT"], &[], &[]);
            let mut exact = Engine::fuzzy(&sets, 0, false, None).unwrap();
            let mut one_off = Engine::fuzzy(&sets, 1, false, None).unwrap();
            assert_eq!(hit(&mut exact, "TTACGAACGTTT", "").0, [] as [usize; 0]);
            assert_eq!(hit(&mut one_off, "TTACGAACGTTT", "").0, [0]);
        }

        #[test]
        fn test_inexact_ignores_exact_hits() {
            let sets = sets(&["ACGTACGT"], &[], &[]);
            let mut inexact = Engine::fuzzy(&sets, 1, true, None).unwrap();
            // only the exact occurrence is present, and it does not count
            assert_eq!(
                hit(&mut inexact, "TTTTACGTACGTTTTT", "")
                    .1
                    .primary
                    .iter()
                    .filter(|s| **s == (4, 12))
                    .count(),
                0
            );
            let mut default = Engine::fuzzy(&sets, 1, false, None).unwrap();
            assert!(hit(&mut default, "TTTTACGTACGTTTTT", "")
                .1
                .primary
                .contains(&(4, 12)));
        }

        // Mirrors https://github.com/RagnarGrootKoerkamp/sassy/issues/66: the Iupac
        // profile treats `N` as a wildcard, so without an N-fraction filter a needle
        // matches a haystack made entirely of `N`s.
        #[test]
        fn test_default_max_n_frac_rejects_all_n_match() {
            let sets = sets(&["ACGTACGTACGT"], &[], &[]);
            let mut e = Engine::fuzzy(&sets, 1, false, None).unwrap();
            assert_eq!(
                hit(&mut e, "NNNNNNNNNNNNNNNNNN", "").0,
                [] as [usize; 0],
                "default max_n_frac (k/pattern_len) should reject an all-N match"
            );
        }

        #[test]
        fn test_max_n_frac_override_allows_all_n_match() {
            let sets = sets(&["ACGTACGTACGT"], &[], &[]);
            let mut e = Engine::fuzzy(&sets, 1, false, Some(1.0)).unwrap();
            assert_eq!(
                hit(&mut e, "NNNNNNNNNNNNNNNNNN", "").0,
                [0],
                "max_n_frac=1.0 should disable the N-fraction filter"
            );
        }

        // sassy's `Searcher::encode_patterns` panics (`assert!`) when a pattern set
        // contains mixed lengths; these tests confirm we catch that up front and
        // return an `Err` instead of letting the panic reach the caller.
        #[test]
        fn test_rejects_mismatched_pattern_lengths_in_any_set() {
            let mixed = ["AAAA", "AAAAA"];
            for sets in [
                sets(&mixed, &[], &[]),
                sets(&[], &mixed, &[]),
                sets(&[], &[], &mixed),
            ] {
                assert!(
                    Engine::fuzzy(&sets, 1, false, None).is_err(),
                    "mismatched pattern lengths should error, not panic"
                );
            }
        }

        #[test]
        fn test_accepts_uniform_pattern_lengths() {
            let sets = sets(&["AAAA", "TTTT", "CCCC"], &[], &[]);
            assert!(Engine::fuzzy(&sets, 1, false, None).is_ok());
        }

        /// Different sets may use different lengths.
        #[test]
        fn test_pattern_lengths_are_per_set() {
            let sets = sets(&["AAAA"], &["CCCCCC"], &["GGGGGGGG"]);
            assert!(Engine::fuzzy(&sets, 1, false, None).is_ok());
        }
    }
}
