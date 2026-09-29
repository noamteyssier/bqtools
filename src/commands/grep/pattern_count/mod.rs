use fixedbitset::FixedBitSet;

use super::{Engine, PatternSets};

mod processor;
pub use processor::PatternCountProcessor;

/// Counts, per pattern, the records that contain it (or, inverted, do not).
///
/// Counts are indexed over the concatenated (primary-only, extended-only,
/// either) pattern sets.
#[derive(Clone)]
pub struct PatternCounter {
    engine: Engine,
    bits: FixedBitSet,
    invert: bool,
    names: Vec<String>,
}
impl PatternCounter {
    pub fn new(engine: Engine, patterns: &PatternSets, invert: bool) -> Self {
        Self {
            bits: engine.bitset(),
            engine,
            invert,
            names: patterns.names(),
        }
    }

    /// Increments the count of every pattern found in the primary or secondary
    /// sequence (every pattern not found, when inverted). A pattern counts once
    /// per record however often it occurs.
    pub fn count_patterns(
        &mut self,
        primary: &[u8],
        secondary: &[u8],
        pattern_count: &mut [usize],
    ) {
        self.bits.clear();
        self.engine.hit(primary, secondary, &mut self.bits, None);
        if self.invert {
            self.bits.zeroes().for_each(|idx| pattern_count[idx] += 1);
        } else {
            self.bits.ones().for_each(|idx| pattern_count[idx] += 1);
        }
    }

    pub fn num_patterns(&self) -> usize {
        self.names.len()
    }

    /// Pattern names (FASTA headers if present, otherwise the pattern strings).
    pub fn pattern_names(&self) -> &[String] {
        &self.names
    }
}

#[cfg(test)]
mod pattern_count_tests {
    use anyhow::Result;

    use super::PatternCounter;
    use crate::commands::grep::{Engine, Pattern, PatternCollection, PatternSets};

    fn counter(
        engine: impl FnOnce(&PatternSets) -> Result<Engine>,
        pat1: PatternCollection,
        pat2: PatternCollection,
        pat: PatternCollection,
        invert: bool,
    ) -> Result<PatternCounter> {
        let sets = PatternSets { pat1, pat2, pat };
        Ok(PatternCounter::new(engine(&sets)?, &sets, invert))
    }

    fn regex_counter(
        pat1: PatternCollection,
        pat2: PatternCollection,
        pat: PatternCollection,
        invert: bool,
    ) -> Result<PatternCounter> {
        counter(Engine::regex, pat1, pat2, pat, invert)
    }

    fn ac_counter(
        pat1: PatternCollection,
        pat2: PatternCollection,
        pat: PatternCollection,
        no_dfa: bool,
        invert: bool,
    ) -> Result<PatternCounter> {
        counter(|s| Engine::aho_corasick(s, no_dfa), pat1, pat2, pat, invert)
    }

    #[cfg(feature = "fuzzy")]
    fn fuzzy_counter(
        pat1: PatternCollection,
        pat2: PatternCollection,
        pat: PatternCollection,
        k: usize,
        inexact: bool,
        invert: bool,
        max_n_frac: Option<f32>,
    ) -> Result<PatternCounter> {
        counter(
            |s| Engine::fuzzy(s, k, inexact, max_n_frac),
            pat1,
            pat2,
            pat,
            invert,
        )
    }

    fn pc(patterns: &[&[u8]]) -> PatternCollection {
        PatternCollection(
            patterns
                .iter()
                .map(|p| Pattern {
                    name: None,
                    sequence: p.to_vec(),
                })
                .collect(),
        )
    }

    #[test]
    fn test_regex_pattern_counter_single_pattern() {
        let mut counter = regex_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), false).unwrap();

        assert_eq!(counter.num_patterns(), 1);

        let primary = b"GGGGAAAATTTT";
        let secondary = b"GGGGCCCCTTTT";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Pattern should be found in primary");
    }

    #[test]
    fn test_regex_pattern_counter_multiple_patterns() {
        let mut counter =
            regex_counter(pc(&[b"AAAA", b"TTTT", b"CCCC"]), pc(&[]), pc(&[]), false).unwrap();

        assert_eq!(counter.num_patterns(), 3);

        let primary = b"AAAAGGGGTTTT";
        let secondary = b"GGGGCCCCGGGG";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "First pattern found");
        assert_eq!(counts[1], 1, "Second pattern found");
        assert_eq!(counts[2], 0, "Third pattern not found in primary");
    }

    #[test]
    fn test_regex_pattern_counter_secondary() {
        let mut counter = regex_counter(pc(&[]), pc(&[b"TTTT"]), pc(&[]), false).unwrap();

        let primary = b"GGGGAAAACCCC";
        let secondary = b"GGGGTTTTCCCC";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Pattern should be found in secondary");
    }

    #[test]
    fn test_regex_pattern_counter_either() {
        let mut counter = regex_counter(pc(&[]), pc(&[]), pc(&[b"CCCC"]), false).unwrap();

        // Test match in primary
        let primary1 = b"GGGGCCCCTTTT";
        let secondary1 = b"GGGGAAAATTTT";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1);

        // Test match in secondary
        let primary2 = b"GGGGAAAATTTT";
        let secondary2 = b"GGGGCCCCTTTT";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 1);

        // Test match in both (should still count as 1)
        let primary3 = b"GGGGCCCCTTTT";
        let secondary3 = b"GGGGCCCCTTTT";
        let mut counts3 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary3, secondary3, &mut counts3);
        assert_eq!(counts3[0], 1);
    }

    #[test]
    fn test_regex_pattern_counter_no_match() {
        let mut counter = regex_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), false).unwrap();

        let primary = b"GGGGCCCCTTTT";
        let secondary = b"GGGGCCCCTTTT";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 0, "Pattern should not be found");
    }

    #[test]
    fn test_regex_pattern_counter_invert() {
        let mut counter = regex_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), true).unwrap();

        // Sequence without pattern (should count when inverted)
        let primary1 = b"GGGGCCCCTTTT";
        let secondary1 = b"GGGGCCCCTTTT";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1, "Should count when pattern not found");

        // Sequence with pattern (should not count when inverted)
        let primary2 = b"GGGGAAAATTTT";
        let secondary2 = b"GGGGCCCCTTTT";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 0, "Should not count when pattern found");
    }

    #[test]
    fn test_regex_pattern_counter_combined_patterns() {
        let mut counter =
            regex_counter(pc(&[b"AAAA"]), pc(&[b"TTTT"]), pc(&[b"CCCC"]), false).unwrap();

        assert_eq!(counter.num_patterns(), 3);

        let primary = b"AAAACCCCGGGG";
        let secondary = b"GGGGTTTTCCCC";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Primary pattern found");
        assert_eq!(counts[1], 1, "Secondary pattern found");
        assert_eq!(counts[2], 1, "Either pattern found");
    }

    #[test]
    fn test_regex_pattern_counter_pattern_names() {
        let counter = regex_counter(pc(&[b"AAAA", b"TTTT"]), pc(&[]), pc(&[]), false).unwrap();

        let patterns = counter.pattern_names();
        assert_eq!(patterns.len(), 2);
        assert_eq!(patterns[0], "AAAA");
        assert_eq!(patterns[1], "TTTT");
    }

    #[test]
    fn test_regex_pattern_counter_multiple_records() {
        let mut counter = regex_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), false).unwrap();

        let mut counts = vec![0; counter.num_patterns()];

        // Process multiple records
        let records = vec![
            (b"GGGGAAAATTTT" as &[u8], b"" as &[u8]),
            (b"GGGGCCCCTTTT", b""),
            (b"AAAACCCCGGGG", b""),
            (b"GGGGCCCCGGGG", b""),
            (b"AAAAAAAAAAAA", b""),
        ];

        for (primary, secondary) in records {
            counter.count_patterns(primary, secondary, &mut counts);
        }

        assert_eq!(
            counts[0], 3,
            "Pattern should be found in 3 out of 5 records"
        );
    }

    #[test]
    fn test_regex_pattern_counter_empty_sequence() {
        let mut counter = regex_counter(pc(&[]), pc(&[b"AAAA"]), pc(&[]), false).unwrap();

        let primary = b"GGGGAAAATTTT";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 0, "Should not count in empty secondary");
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_single_pattern() {
        let mut counter =
            fuzzy_counter(pc(&[b"AAAAAAAA"]), pc(&[]), pc(&[]), 1, false, false, None).unwrap();

        assert_eq!(counter.num_patterns(), 1);

        let primary = b"GGGGAAAAAAAATTTT";
        let secondary = b"GGGGCCCCTTTT";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Pattern should be found in primary");
    }

    // Mirrors https://github.com/RagnarGrootKoerkamp/sassy/issues/66: the Iupac
    // profile treats `N` as a wildcard, so without an N-fraction filter a needle
    // matches a haystack made entirely of `N`s.
    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_default_max_n_frac_rejects_all_n_match() {
        let mut counter = fuzzy_counter(
            pc(&[b"ACGTACGTACGT"]),
            pc(&[]),
            pc(&[]),
            1,
            false,
            false,
            None,
        )
        .unwrap();

        let primary = b"NNNNNNNNNNNNNNNNNN";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];
        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(
            counts[0], 0,
            "default max_n_frac (k/pattern_len) should reject an all-N match"
        );
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_max_n_frac_override_allows_all_n_match() {
        let mut counter = fuzzy_counter(
            pc(&[b"ACGTACGTACGT"]),
            pc(&[]),
            pc(&[]),
            1,
            false,
            false,
            Some(1.0),
        )
        .unwrap();

        let primary = b"NNNNNNNNNNNNNNNNNN";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];
        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(
            counts[0], 1,
            "max_n_frac=1.0 should disable the N-fraction filter"
        );
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_max_n_frac_explicit_zero_rejects_any_n() {
        let mut counter = fuzzy_counter(
            pc(&[b"AAAAAAAA"]),
            pc(&[]),
            pc(&[]),
            1,
            false,
            false,
            Some(0.0),
        )
        .unwrap();

        // A single N substituted for an A is within edit distance 1, and would
        // pass the default k/pattern_len threshold (1/8), but max_n_frac=0.0
        // should reject any match containing an N.
        let primary_with_n = b"GGGGAAAAAAANTTTT";
        let mut counts = vec![0; counter.num_patterns()];
        counter.count_patterns(primary_with_n, b"", &mut counts);
        assert_eq!(
            counts[0], 0,
            "max_n_frac=0.0 should reject a match containing any N"
        );

        // A match with no N's at all should still be counted.
        let primary_no_n = b"GGGGAAAAAAAATTTT";
        let mut counts_clean = vec![0; counter.num_patterns()];
        counter.count_patterns(primary_no_n, b"", &mut counts_clean);
        assert_eq!(counts_clean[0], 1);
    }

    // sassy's `Searcher::encode_patterns` panics (`assert!`) when a pattern set
    // contains mixed lengths; these tests confirm we catch that up front and
    // return an `Err` instead of letting the panic reach the caller.
    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_rejects_mismatched_pattern_lengths_primary() {
        let result = fuzzy_counter(
            pc(&[b"AAAA", b"AAAAA"]),
            pc(&[]),
            pc(&[]),
            1,
            false,
            false,
            None,
        );
        assert!(
            result.is_err(),
            "mismatched primary pattern lengths should error, not panic"
        );
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_rejects_mismatched_pattern_lengths_secondary() {
        let result = fuzzy_counter(
            pc(&[]),
            pc(&[b"AAAA", b"AAAAA"]),
            pc(&[]),
            1,
            false,
            false,
            None,
        );
        assert!(
            result.is_err(),
            "mismatched secondary pattern lengths should error, not panic"
        );
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_rejects_mismatched_pattern_lengths_either() {
        let result = fuzzy_counter(
            pc(&[]),
            pc(&[]),
            pc(&[b"AAAA", b"AAAAA"]),
            1,
            false,
            false,
            None,
        );
        assert!(
            result.is_err(),
            "mismatched either-set pattern lengths should error, not panic"
        );
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_accepts_uniform_pattern_lengths() {
        let result = fuzzy_counter(
            pc(&[b"AAAA", b"TTTT", b"CCCC"]),
            pc(&[]),
            pc(&[]),
            1,
            false,
            false,
            None,
        );
        assert!(result.is_ok(), "uniform pattern lengths should not error");
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_with_mismatches() {
        let mut counter =
            fuzzy_counter(pc(&[b"AAAAAAAA"]), pc(&[]), pc(&[]), 2, false, false, None).unwrap();

        // Exact match
        let primary1 = b"GGGGAAAAAAAATTTT";
        let secondary1 = b"";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1);

        // One mismatch
        let primary2 = b"GGGGAAAAACAATTTT";
        let secondary2 = b"";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 1);

        // Two mismatches
        let primary3 = b"GGGGAAAACCAATTTT";
        let secondary3 = b"";
        let mut counts3 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary3, secondary3, &mut counts3);
        assert_eq!(counts3[0], 1);
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_inexact_only() {
        let mut counter =
            fuzzy_counter(pc(&[b"AAAAAAAA"]), pc(&[]), pc(&[]), 2, true, false, None).unwrap();

        // Exact match (should not count with inexact_only)
        let primary1 = b"GGGGAAAAAAAATTTT";
        let secondary1 = b"";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(
            counts1[0], 0,
            "Exact match should not count with inexact_only"
        );

        // Inexact match (should count)
        let primary2 = b"GGGGAAAAACAATTTT";
        let secondary2 = b"";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 1, "Inexact match should count");
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_invert() {
        let mut counter =
            fuzzy_counter(pc(&[b"AAAAAAAA"]), pc(&[]), pc(&[]), 1, false, true, None).unwrap();

        // Sequence without pattern (should count when inverted)
        let primary1 = b"GGGGCCCCTTTT";
        let secondary1 = b"";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1, "Should count when pattern not found");

        // Sequence with pattern (should not count when inverted)
        let primary2 = b"GGGGAAAAAAAATTTT";
        let secondary2 = b"";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 0, "Should not count when pattern found");
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_multiple_patterns() {
        let mut counter = fuzzy_counter(
            pc(&[b"AAAAAAAA", b"TTTTTTTT", b"CCCCCCCC"]),
            pc(&[]),
            pc(&[]),
            1,
            false,
            false,
            None,
        )
        .unwrap();

        assert_eq!(counter.num_patterns(), 3);

        let primary = b"AAAAAAAATTTTTTTT";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "First pattern found");
        assert_eq!(counts[1], 1, "Second pattern found");
        assert_eq!(counts[2], 0, "Third pattern not found");
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_secondary() {
        let mut counter =
            fuzzy_counter(pc(&[]), pc(&[b"TTTTTTTT"]), pc(&[]), 1, false, false, None).unwrap();

        let primary = b"GGGGAAAACCCC";
        let secondary = b"GGGGTTTTTTTTCCCC";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Pattern should be found in secondary");
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_either() {
        let mut counter =
            fuzzy_counter(pc(&[]), pc(&[]), pc(&[b"CCCCCCCC"]), 1, false, false, None).unwrap();

        // Test match in primary
        let primary1 = b"GGGGCCCCCCCCTTTT";
        let secondary1 = b"GGGGAAAATTTT";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1);

        // Test match in secondary
        let primary2 = b"GGGGAAAATTTT";
        let secondary2 = b"GGGGCCCCCCCCTTTT";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 1);
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_pattern_names() {
        let counter = fuzzy_counter(
            pc(&[b"AAAAAAAA", b"TTTTTTTT"]),
            pc(&[]),
            pc(&[]),
            1,
            false,
            false,
            None,
        )
        .unwrap();

        let patterns = counter.pattern_names();
        assert_eq!(patterns.len(), 2);
        assert_eq!(patterns[0], "AAAAAAAA");
        assert_eq!(patterns[1], "TTTTTTTT");
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_fuzzy_pattern_counter_edit_distance_zero() {
        let mut counter =
            fuzzy_counter(pc(&[b"AAAAAAAA"]), pc(&[]), pc(&[]), 0, false, false, None).unwrap();

        // Exact match (should count)
        let primary1 = b"GGGGAAAAAAAATTTT";
        let secondary1 = b"";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1);

        // One mismatch (should not count with k=0)
        let primary2 = b"GGGGAAAAACAATTTT";
        let secondary2 = b"";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 0);
    }

    #[test]
    fn test_aho_corasick_pattern_counter_single_pattern() {
        let mut counter = ac_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), false, false).unwrap();

        assert_eq!(counter.num_patterns(), 1);

        let primary = b"GGGGAAAATTTT";
        let secondary = b"GGGGCCCCTTTT";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Pattern should be found in primary");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_multiple_patterns() {
        let mut counter = ac_counter(
            pc(&[b"AAAA", b"TTTT", b"CCCC"]),
            pc(&[]),
            pc(&[]),
            false,
            false,
        )
        .unwrap();

        assert_eq!(counter.num_patterns(), 3);

        let primary = b"AAAAGGGGTTTT";
        let secondary = b"GGGGCCCCGGGG";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "First pattern found");
        assert_eq!(counts[1], 1, "Second pattern found");
        assert_eq!(counts[2], 0, "Third pattern not found in primary");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_secondary() {
        let mut counter = ac_counter(pc(&[]), pc(&[b"TTTT"]), pc(&[]), false, false).unwrap();

        let primary = b"GGGGAAAACCCC";
        let secondary = b"GGGGTTTTCCCC";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Pattern should be found in secondary");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_either() {
        let mut counter = ac_counter(pc(&[]), pc(&[]), pc(&[b"CCCC"]), false, false).unwrap();

        // Test match in primary
        let primary1 = b"GGGGCCCCTTTT";
        let secondary1 = b"GGGGAAAATTTT";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1);

        // Test match in secondary
        let primary2 = b"GGGGAAAATTTT";
        let secondary2 = b"GGGGCCCCTTTT";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 1);

        // Test match in both (should still count as 1)
        let primary3 = b"GGGGCCCCTTTT";
        let secondary3 = b"GGGGCCCCTTTT";
        let mut counts3 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary3, secondary3, &mut counts3);
        assert_eq!(counts3[0], 1);
    }

    #[test]
    fn test_aho_corasick_pattern_counter_no_match() {
        let mut counter = ac_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), false, false).unwrap();

        let primary = b"GGGGCCCCTTTT";
        let secondary = b"GGGGCCCCTTTT";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 0, "Pattern should not be found");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_invert() {
        let mut counter = ac_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), false, true).unwrap();

        // Sequence without pattern (should count when inverted)
        let primary1 = b"GGGGCCCCTTTT";
        let secondary1 = b"GGGGCCCCTTTT";
        let mut counts1 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary1, secondary1, &mut counts1);
        assert_eq!(counts1[0], 1, "Should count when pattern not found");

        // Sequence with pattern (should not count when inverted)
        let primary2 = b"GGGGAAAATTTT";
        let secondary2 = b"GGGGCCCCTTTT";
        let mut counts2 = vec![0; counter.num_patterns()];
        counter.count_patterns(primary2, secondary2, &mut counts2);
        assert_eq!(counts2[0], 0, "Should not count when pattern found");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_combined_patterns() {
        let mut counter =
            ac_counter(pc(&[b"AAAA"]), pc(&[b"TTTT"]), pc(&[b"CCCC"]), false, false).unwrap();

        assert_eq!(counter.num_patterns(), 3);

        let primary = b"AAAACCCCGGGG";
        let secondary = b"GGGGTTTTCCCC";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "Primary pattern found");
        assert_eq!(counts[1], 1, "Secondary pattern found");
        assert_eq!(counts[2], 1, "Either pattern found");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_pattern_names() {
        let counter = ac_counter(pc(&[b"AAAA", b"TTTT"]), pc(&[]), pc(&[]), false, false).unwrap();

        let patterns = counter.pattern_names();
        assert_eq!(patterns.len(), 2);
        assert_eq!(patterns[0], "AAAA");
        assert_eq!(patterns[1], "TTTT");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_multiple_records() {
        let mut counter = ac_counter(pc(&[b"AAAA"]), pc(&[]), pc(&[]), false, false).unwrap();

        let mut counts = vec![0; counter.num_patterns()];

        // Process multiple records
        let records = vec![
            (b"GGGGAAAATTTT" as &[u8], b"" as &[u8]),
            (b"GGGGCCCCTTTT", b""),
            (b"AAAACCCCGGGG", b""),
            (b"GGGGCCCCGGGG", b""),
            (b"AAAAAAAAAAAA", b""),
        ];

        for (primary, secondary) in records {
            counter.count_patterns(primary, secondary, &mut counts);
        }

        assert_eq!(
            counts[0], 3,
            "Pattern should be found in 3 out of 5 records"
        );
    }

    #[test]
    fn test_aho_corasick_pattern_counter_empty_sequence() {
        let mut counter = ac_counter(pc(&[]), pc(&[b"AAAA"]), pc(&[]), false, false).unwrap();

        let primary = b"GGGGAAAATTTT";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 0, "Should not count in empty secondary");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_overlapping_patterns() {
        let mut counter =
            ac_counter(pc(&[b"AAA", b"AAAA"]), pc(&[]), pc(&[]), false, false).unwrap();

        let primary = b"GGGGAAAAATTTT";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 1, "AAA pattern found");
        assert_eq!(counts[1], 1, "AAAA pattern found");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_multiple_occurrences() {
        let mut counter =
            ac_counter(pc(&[b"AAA", b"AAAA"]), pc(&[]), pc(&[]), false, false).unwrap();

        let primary = b"GGGGAAAAATTTT";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        // Both patterns should be found (overlapping matches)
        assert_eq!(counts[0], 1, "AAA pattern found");
        assert_eq!(counts[1], 1, "AAAA pattern found");
    }

    #[test]
    fn test_aho_corasick_pattern_counter_case_sensitive() {
        let mut counter = ac_counter(pc(&[b"aaaa"]), pc(&[]), pc(&[]), false, false).unwrap();

        // Different case should not match
        let primary = b"GGGGAAAATTTT";
        let secondary = b"";
        let mut counts = vec![0; counter.num_patterns()];

        counter.count_patterns(primary, secondary, &mut counts);

        assert_eq!(counts[0], 0, "Case sensitive pattern should not match");

        // Same case should match
        let primary2 = b"ggggaaaatttt";
        let secondary2 = b"";
        let mut counts2 = vec![0; counter.num_patterns()];

        counter.count_patterns(primary2, secondary2, &mut counts2);

        assert_eq!(counts2[0], 1, "Same case should match");
    }
}
