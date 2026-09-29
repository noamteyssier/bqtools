use anyhow::{bail, Result};
use log::debug;

use crate::cli::Mate;

/// Returns true if the pattern is a fixed DNA string (only ACGT).
pub(crate) fn is_fixed(pattern: &[u8]) -> bool {
    !pattern.is_empty()
        && pattern
            .iter()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

/// A pattern with an optional name (from FASTA headers).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pattern {
    pub name: Option<String>,
    pub sequence: Vec<u8>,
}
impl Pattern {
    /// Compile the sequence into a regex.
    pub fn to_regex(&self) -> Result<regex::bytes::Regex> {
        let seq_str = std::str::from_utf8(&self.sequence)?;
        Ok(regex::bytes::Regex::new(seq_str)?)
    }

    /// Reverse complement the pattern's sequence in place.
    ///
    /// Errors if the sequence is not a fixed literal ACGT string, since
    /// reverse-complementing a regex (or an IUPAC-ambiguous pattern) is undefined.
    pub fn reverse_complement(&mut self) -> Result<()> {
        if !is_fixed(&self.sequence) {
            anyhow::bail!(
                "Cannot reverse complement pattern '{}': --rc only supports fixed ACGT patterns, not regex",
                String::from_utf8_lossy(&self.sequence)
            );
        }
        self.sequence.reverse();
        self.sequence.iter_mut().for_each(|b| {
            *b = match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => unreachable!("is_fixed guarantees only ACGT bytes"),
            }
        });
        Ok(())
    }
}

/// A collection of patterns with convenience methods for type conversions.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PatternCollection(pub Vec<Pattern>);
impl PatternCollection {
    pub fn bytes(&self) -> Vec<Vec<u8>> {
        self.0.iter().map(|p| p.sequence.clone()).collect()
    }

    pub fn regexes(&self) -> Result<Vec<regex::bytes::Regex>> {
        self.0.iter().map(Pattern::to_regex).collect()
    }

    /// Reverse complement every pattern in the collection, in place.
    ///
    /// Errors if any pattern is not a fixed literal ACGT string.
    pub fn reverse_complement(&mut self) -> Result<()> {
        for pattern in &mut self.0 {
            pattern.reverse_complement()?;
        }
        Ok(())
    }

    pub fn names(&self) -> Vec<String> {
        self.0
            .iter()
            .map(|p| {
                p.name.clone().unwrap_or_else(|| {
                    std::str::from_utf8(&p.sequence)
                        .expect("Non-UTF8 sequence in pattern")
                        .to_string()
                })
            })
            .collect()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, Pattern> {
        self.0.iter()
    }

    /// Takes all patterns from `other` and moves them into this collection.
    pub fn ingest(&mut self, other: &mut Self) {
        self.0.append(&mut other.0);
    }

    /// Clears all patterns from this collection.
    pub fn clear(&mut self) {
        self.0.clear();
    }
}
/// The primary-only, extended-only and either-sequence pattern sets.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PatternSets {
    pub pat1: PatternCollection,
    pub pat2: PatternCollection,
    pub pat: PatternCollection,
}
impl PatternSets {
    fn each_mut(&mut self) -> [&mut PatternCollection; 3] {
        [&mut self.pat1, &mut self.pat2, &mut self.pat]
    }

    /// Total number of patterns across all three sets.
    pub fn len(&self) -> usize {
        self.pat1.len() + self.pat2.len() + self.pat.len()
    }

    /// Reverse complement every pattern; errors on any non-ACGT pattern.
    pub fn reverse_complement(&mut self) -> Result<()> {
        self.each_mut()
            .into_iter()
            .try_for_each(PatternCollection::reverse_complement)
    }

    /// Fixed-string (Aho-Corasick) matching applies when forced, or when every
    /// pattern is a plain uppercase ACGT string.
    pub fn use_fixed(&self, forced: bool) -> bool {
        let fixed = forced
            || [&self.pat1, &self.pat2, &self.pat]
                .into_iter()
                .flat_map(PatternCollection::iter)
                .all(|p| is_fixed(&p.sequence));
        if fixed && !forced {
            debug!("All patterns are fixed strings — auto-selecting Aho-Corasick");
        }
        fixed
    }

    /// Applies `--mate`: restricts every pattern to the selected mate, so no
    /// match can occur on an ignored mate.
    pub fn redistribute(&mut self, mate: Mate) -> Result<()> {
        let (keep, drop, n) = match mate {
            Mate::Both => return Ok(()),
            Mate::One => (&mut self.pat1, &mut self.pat2, 1),
            Mate::Two => (&mut self.pat2, &mut self.pat1, 2),
        };
        drop.clear();
        keep.ingest(&mut self.pat);
        if keep.is_empty() {
            bail!("No patterns provided for mate {n}");
        }
        Ok(())
    }
}

impl IntoIterator for PatternCollection {
    type Item = Pattern;
    type IntoIter = std::vec::IntoIter<Pattern>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

#[cfg(test)]
mod reverse_complement_tests {
    use super::Pattern;

    fn pattern(seq: &[u8]) -> Pattern {
        Pattern {
            name: None,
            sequence: seq.to_vec(),
        }
    }

    #[test]
    fn test_reverse_complement_basic() {
        let mut p = pattern(b"ACGT");
        p.reverse_complement().unwrap();
        assert_eq!(p.sequence, b"ACGT");

        let mut p = pattern(b"AACCGGTT");
        p.reverse_complement().unwrap();
        assert_eq!(p.sequence, b"AACCGGTT");

        let mut p = pattern(b"AAAA");
        p.reverse_complement().unwrap();
        assert_eq!(p.sequence, b"TTTT");

        let mut p = pattern(b"ACGTACGT");
        p.reverse_complement().unwrap();
        assert_eq!(p.sequence, b"ACGTACGT");

        let mut p = pattern(b"AGGT");
        p.reverse_complement().unwrap();
        assert_eq!(p.sequence, b"ACCT");
    }

    #[test]
    fn test_reverse_complement_preserves_name() {
        let mut p = Pattern {
            name: Some("my_pattern".to_string()),
            sequence: b"AGGT".to_vec(),
        };
        p.reverse_complement().unwrap();
        assert_eq!(p.name, Some("my_pattern".to_string()));
        assert_eq!(p.sequence, b"ACCT");
    }

    #[test]
    fn test_reverse_complement_rejects_regex() {
        assert!(pattern(b"AC.GT").reverse_complement().is_err());
        assert!(pattern(b"AC[GT]").reverse_complement().is_err());
        assert!(pattern(b"A{3}").reverse_complement().is_err());
        assert!(pattern(b"^ACGT").reverse_complement().is_err());
    }

    #[test]
    fn test_reverse_complement_rejects_iupac_ambiguity() {
        assert!(pattern(b"ACGN").reverse_complement().is_err());
    }

    #[test]
    fn test_reverse_complement_rejects_empty() {
        assert!(pattern(b"").reverse_complement().is_err());
    }
}

#[cfg(test)]
mod pattern_sets_tests {
    use super::{is_fixed, Mate, Pattern, PatternCollection, PatternSets};

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

    fn sets(pat1: &[&[u8]], pat2: &[&[u8]], pat: &[&[u8]]) -> PatternSets {
        PatternSets {
            pat1: pc(pat1),
            pat2: pc(pat2),
            pat: pc(pat),
        }
    }

    #[test]
    fn test_is_fixed() {
        for p in [&b"ACGTACGT"[..], b"AAAAAAAAAA", b"ACGT"] {
            assert!(is_fixed(p));
        }
        // empty, IUPAC, lowercase and regex syntax are not fixed
        for p in [
            &b""[..],
            b"ACGTNRYW",
            b"ACGN",
            b"acgt",
            b"AC.GT",
            b"AC[GT]",
            b"A{3}",
            b"^ACGT",
            b"ACG|TGA",
            b"(ACG)",
            b"AC\\dGT",
        ] {
            assert!(!is_fixed(p), "{}", String::from_utf8_lossy(p));
        }
    }

    #[test]
    fn test_use_fixed() {
        assert!(sets(&[b"ACGT", b"TTTT"], &[b"GGGG"], &[]).use_fixed(false));
        assert!(!sets(&[b"ACGT", b"AC.GT"], &[b"GGGG"], &[]).use_fixed(false));
        assert!(sets(&[], &[], &[]).use_fixed(false));
        // forcing wins over autodetection
        assert!(sets(&[b"AC.GT"], &[], &[]).use_fixed(true));
    }

    #[test]
    fn test_redistribution_noop() {
        let mut s = sets(&[b"ACGT", b"TTTT"], &[b"GGGG"], &[b"AC.GT"]);
        let before = s.clone();
        s.redistribute(Mate::Both).unwrap();
        assert_eq!(s, before);
    }

    #[test]
    fn test_redistribution_m1() {
        let mut s = sets(&[b"ACGT", b"TTTT"], &[b"GGGG"], &[b"AC.GT"]);
        s.redistribute(Mate::One).unwrap();
        assert_eq!(s, sets(&[b"ACGT", b"TTTT", b"AC.GT"], &[], &[]));
    }

    #[test]
    fn test_redistribution_m2() {
        let mut s = sets(&[b"ACGT", b"TTTT"], &[b"GGGG"], &[b"AC.GT"]);
        s.redistribute(Mate::Two).unwrap();
        assert_eq!(s, sets(&[], &[b"GGGG", b"AC.GT"], &[]));
    }

    #[test]
    fn test_redistribution_errors_without_patterns() {
        assert!(sets(&[], &[b"GGGG"], &[]).redistribute(Mate::One).is_err());
        assert!(sets(&[b"GGGG"], &[], &[]).redistribute(Mate::Two).is_err());
    }
}
