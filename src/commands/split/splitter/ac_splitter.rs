use aho_corasick::{AhoCorasick, AhoCorasickBuilder, AhoCorasickKind};
use anyhow::Result;
use fixedbitset::FixedBitSet;

use crate::commands::grep::PatternCollection;

/// Fixed-string matching with Aho-Corasick over the three pattern sets
/// (primary-only, secondary-only, either).
#[derive(Clone)]
pub struct AhoCorasickSplitter {
    state1: AhoCorasick,
    state2: AhoCorasick,
    state: AhoCorasick,
}

impl AhoCorasickSplitter {
    pub fn new(
        pat1: &PatternCollection,
        pat2: &PatternCollection,
        pat: &PatternCollection,
        no_dfa: bool,
    ) -> Result<Self> {
        Ok(Self {
            state1: corasick_builder(&pat1.bytes(), no_dfa)?,
            state2: corasick_builder(&pat2.bytes(), no_dfa)?,
            state: corasick_builder(&pat.bytes(), no_dfa)?,
        })
    }

    /// Sets a bit in `bits` for every pattern found in the sequences.
    pub fn hit(&self, primary: &[u8], secondary: &[u8], bits: &mut FixedBitSet) {
        let n1 = self.state1.patterns_len();
        let n2 = self.state2.patterns_len();
        match_patterns(&self.state1, bits, &[primary], 0);
        match_patterns(&self.state2, bits, &[secondary], n1);
        match_patterns(&self.state, bits, &[primary, secondary], n1 + n2);
    }
}

fn match_patterns(patterns: &AhoCorasick, bitset: &mut FixedBitSet, seqs: &[&[u8]], offset: usize) {
    if patterns.patterns_len() == 0 {
        return;
    }
    for seq in seqs.iter().filter(|s| !s.is_empty()) {
        for m in patterns.find_overlapping_iter(seq) {
            bitset.insert(offset + m.pattern().as_usize());
        }
    }
}

fn corasick_builder(patterns: &[Vec<u8>], no_dfa: bool) -> Result<AhoCorasick> {
    Ok(AhoCorasickBuilder::new()
        .ascii_case_insensitive(false)
        .kind(if no_dfa {
            None
        } else {
            Some(AhoCorasickKind::DFA)
        })
        .build(patterns)?)
}
