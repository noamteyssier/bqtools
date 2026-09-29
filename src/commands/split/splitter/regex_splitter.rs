use anyhow::Result;
use fixedbitset::FixedBitSet;

use crate::commands::grep::PatternCollection;

type Expressions = Vec<regex::bytes::Regex>;

/// Regular-expression matching over the three pattern sets
/// (primary-only, secondary-only, either).
#[derive(Clone)]
pub struct RegexSplitter {
    re1: Expressions,
    re2: Expressions,
    re: Expressions,
}

impl RegexSplitter {
    pub fn new(
        pat1: &PatternCollection,
        pat2: &PatternCollection,
        pat: &PatternCollection,
    ) -> Result<Self> {
        Ok(Self {
            re1: pat1.regexes()?,
            re2: pat2.regexes()?,
            re: pat.regexes()?,
        })
    }

    /// Sets a bit in `bits` for every pattern found in the sequences.
    pub fn hit(&self, primary: &[u8], secondary: &[u8], bits: &mut FixedBitSet) {
        let n1 = self.re1.len();
        let n2 = self.re2.len();
        match_patterns(&self.re1, bits, &[primary], 0);
        match_patterns(&self.re2, bits, &[secondary], n1);
        match_patterns(&self.re, bits, &[primary, secondary], n1 + n2);
    }
}

fn match_patterns(patterns: &Expressions, bitset: &mut FixedBitSet, seqs: &[&[u8]], offset: usize) {
    for seq in seqs.iter().filter(|s| !s.is_empty()) {
        for (idx, reg) in patterns.iter().enumerate() {
            if reg.is_match(seq) {
                bitset.insert(offset + idx);
            }
        }
    }
}
