mod ac_splitter;
#[cfg(feature = "fuzzy")]
mod fuzzy_splitter;
mod processor;
mod regex_splitter;

use fixedbitset::FixedBitSet;
use hashbrown::HashMap;

pub use ac_splitter::AhoCorasickSplitter;
#[cfg(feature = "fuzzy")]
pub use fuzzy_splitter::FuzzySplitter;
pub use processor::SplitProcessor;
pub use regex_splitter::RegexSplitter;

use crate::commands::grep::PatternCollection;

/// Maps hits over the concatenated pattern sets (primary, secondary, either)
/// onto unique aliases, and resolves a record to a single output bin.
#[derive(Clone)]
struct AliasBins {
    /// unique aliases across all pattern sets, ordered by bin index
    aliases: Vec<String>,

    /// alias (bin) index of the pattern at each global pattern index
    alias_idx: Vec<usize>,

    /// bitset over all patterns
    all_bits: FixedBitSet,

    /// bitset over all unique aliases
    unique_bits: FixedBitSet,
}
impl AliasBins {
    fn new(sets: [&PatternCollection; 3]) -> Self {
        let mut aliases = Vec::new();
        let mut alias_idx = Vec::new();
        let mut map = HashMap::new();
        for name in sets.iter().flat_map(|s| s.names()) {
            let idx = *map.entry(name.clone()).or_insert_with(|| {
                aliases.push(name);
                aliases.len() - 1
            });
            alias_idx.push(idx);
        }
        Self {
            all_bits: FixedBitSet::with_capacity(alias_idx.len()),
            unique_bits: FixedBitSet::with_capacity(aliases.len()),
            aliases,
            alias_idx,
        }
    }

    /// Returns the bin index when the pattern hits in `all_bits` resolve to
    /// exactly one unique alias, otherwise `None`.
    fn resolve(&mut self) -> Option<usize> {
        self.unique_bits.clear();
        for idx in self.all_bits.ones() {
            if let Some(u_idx) = self.alias_idx.get(idx) {
                self.unique_bits.insert(*u_idx);
            }
        }
        let mut hits = self.unique_bits.ones();
        hits.next().filter(|_| hits.next().is_none())
    }
}

/// The matching strategy behind a [`Splitter`].
#[derive(Clone)]
pub enum Matcher {
    AhoCorasick(AhoCorasickSplitter),
    Regex(RegexSplitter),
    #[cfg(feature = "fuzzy")]
    Fuzzy(Box<FuzzySplitter>),
}

/// Resolves a record's (primary, secondary) sequences to a single output bin.
///
/// Patterns are matched against the primary sequence (first set), the secondary
/// sequence (second set), and either sequence (third set). A record is assigned
/// to a bin only when its matches resolve to exactly one unique alias.
#[derive(Clone)]
pub struct Splitter {
    bins: AliasBins,
    matcher: Matcher,
}
impl Splitter {
    pub fn new(
        matcher: Matcher,
        pat1: &PatternCollection,
        pat2: &PatternCollection,
        pat: &PatternCollection,
    ) -> Self {
        Self {
            bins: AliasBins::new([pat1, pat2, pat]),
            matcher,
        }
    }

    /// Returns the bin index a record belongs to, or `None` when the record's
    /// matches do not resolve to exactly one unique alias.
    pub fn split_idx(&mut self, primary: &[u8], secondary: &[u8]) -> Option<usize> {
        let bits = &mut self.bins.all_bits;
        bits.clear();
        match &mut self.matcher {
            Matcher::AhoCorasick(m) => m.hit(primary, secondary, bits),
            Matcher::Regex(m) => m.hit(primary, secondary, bits),
            #[cfg(feature = "fuzzy")]
            Matcher::Fuzzy(m) => m.hit(primary, secondary, bits),
        }
        self.bins.resolve()
    }

    /// The unique aliases records can be split into, ordered by bin index.
    pub fn aliases(&self) -> &[String] {
        &self.bins.aliases
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::commands::grep::Pattern;

    fn pc(pats: &[(&str, &str)]) -> PatternCollection {
        PatternCollection(
            pats.iter()
                .map(|(n, s)| Pattern {
                    name: Some((*n).to_string()),
                    sequence: s.as_bytes().to_vec(),
                })
                .collect(),
        )
    }

    fn ac(p1: &PatternCollection, p2: &PatternCollection, p: &PatternCollection) -> Splitter {
        let m = AhoCorasickSplitter::new(p1, p2, p, false).unwrap();
        Splitter::new(Matcher::AhoCorasick(m), p1, p2, p)
    }

    #[test]
    fn test_unique_alias_hit_semantics() {
        let none = pc(&[]);
        // two patterns share alias "a"; "c" is a separate alias
        let p = pc(&[("a", "AAAA"), ("a", "CCCC"), ("c", "GGGG")]);
        let mut s = ac(&none, &none, &p);
        assert_eq!(s.aliases(), ["a", "c"]);
        // no hit
        assert_eq!(s.split_idx(b"TTTTTT", b""), None);
        // one pattern hit
        assert_eq!(s.split_idx(b"TAAAAT", b""), Some(0));
        // two patterns of the same alias still resolve to that alias
        assert_eq!(s.split_idx(b"AAAACCCC", b""), Some(0));
        // hits from either sequence resolve the same way
        assert_eq!(s.split_idx(b"TTTT", b"GGGG"), Some(1));
        // two distinct aliases are ambiguous
        assert_eq!(s.split_idx(b"AAAAGGGG", b""), None);
        assert_eq!(s.split_idx(b"AAAA", b"GGGG"), None);
    }

    #[test]
    fn test_empty_sequences_skipped() {
        let none = pc(&[]);
        let p = pc(&[("a", "AAAA")]);
        let mut s = ac(&none, &none, &p);
        assert_eq!(s.split_idx(b"", b""), None);
        assert_eq!(s.split_idx(b"", b"AAAA"), Some(0));
        assert_eq!(s.split_idx(b"AAAA", b""), Some(0));
        let mut re = Splitter::new(
            Matcher::Regex(RegexSplitter::new(&none, &none, &p).unwrap()),
            &none,
            &none,
            &p,
        );
        assert_eq!(re.split_idx(b"", b""), None);
        assert_eq!(re.split_idx(b"", b"AAAA"), Some(0));
    }

    #[test]
    fn test_set_offsets() {
        let none = pc(&[]);
        let p1 = pc(&[("x", "AAAA")]);
        let p2 = pc(&[("y", "CCCC")]);
        let p = pc(&[("z", "GGGG")]);
        let mut s = ac(&p1, &p2, &p);
        assert_eq!(s.aliases(), ["x", "y", "z"]);
        // primary-only pattern must not match the secondary sequence
        assert_eq!(s.split_idx(b"", b"AAAA"), None);
        assert_eq!(s.split_idx(b"AAAA", b""), Some(0));
        assert_eq!(s.split_idx(b"", b"CCCC"), Some(1));
        assert_eq!(s.split_idx(b"CCCC", b""), None);
        assert_eq!(s.split_idx(b"", b"GGGG"), Some(2));
        let _ = none;
    }
}
