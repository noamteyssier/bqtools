mod splitter;

use anyhow::Result;
use binseq::BinseqReader;

use splitter::{SplitProcessor, Splitter};

use crate::{
    cli::SplitCommand,
    commands::{
        grep::{Engine, PatternSets},
        utils::builder_from_reader,
    },
};

/// Loads the primary-only, secondary-only and either-sequence pattern sets.
fn load_patterns(args: &SplitCommand) -> Result<PatternSets> {
    let mut patterns = args.patterns.load_all_patterns()?;
    if args.split.rc {
        patterns.reverse_complement()?;
    }
    Ok(patterns)
}

/// Selects and builds the splitter backend.
///
/// Fuzzy matching (`-z/--fuzzy`) takes priority when enabled. Otherwise,
/// fixed-string pattern sets use the Aho-Corasick backend (auto-detected, or
/// forced with `-x/--fixed`); anything else falls back to the regex backend.
fn build_splitter(args: &SplitCommand) -> Result<Splitter> {
    let patterns = load_patterns(args)?;

    #[cfg(feature = "fuzzy")]
    if args.fuzzy_args.fuzzy {
        log::trace!(
            "Using fuzzy splitter backend (k={}, inexact={}, backend=sassy)",
            args.fuzzy_args.distance,
            args.fuzzy_args.inexact,
        );
        let engine = Engine::fuzzy(
            &patterns,
            args.fuzzy_args.distance,
            args.fuzzy_args.inexact,
            args.fuzzy_args.max_n_frac,
        )?;
        return Ok(Splitter::new(engine, &patterns));
    }

    let engine = if patterns.use_fixed(args.split.fixed) {
        log::trace!(
            "Using Aho-Corasick splitter backend (dfa={})",
            !args.split.no_dfa,
        );
        Engine::aho_corasick(&patterns, args.split.no_dfa)?
    } else {
        log::trace!("Using regex splitter backend");
        Engine::regex(&patterns)?
    };
    Ok(Splitter::new(engine, &patterns))
}

pub fn run(args: &SplitCommand) -> Result<()> {
    let splitter = build_splitter(args)?;
    let reader = BinseqReader::new(args.input.path())?;
    let range = args.input.range(reader.num_records()?)?;
    let builder = builder_from_reader(&reader);
    std::fs::create_dir_all(&args.split.basepath)?;
    let mut proc = SplitProcessor::new(splitter, &builder, args)?;
    reader.process_parallel_range(proc.clone(), args.split.threads, range)?;
    proc.finish()?;
    if !args.split.quiet {
        proc.pprint_counts()?;
    }
    if args.split.min_records > 0 {
        let removed = proc.prune_below(args.split.min_records)?;
        if removed > 0 {
            log::debug!("Removed {removed} output file(s) below the record threshold");
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::Write as _;

    use anyhow::Result;
    use clap::Parser;
    use tempfile::NamedTempFile;

    use crate::cli::BinseqMode;
    use crate::testutils::{count_binseq, write_fastx, DEFAULT_NUM_RECORDS};

    fn encode(in_path: &std::path::Path, out_path: &std::path::Path) -> Result<()> {
        let cmd = crate::cli::EncodeCommand::try_parse_from([
            "encode",
            in_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])?;
        crate::commands::encode::run(&cmd)
    }

    /// Write a plain-text pattern file (one pattern per line).
    fn write_patterns(patterns: &[&str]) -> Result<NamedTempFile> {
        let tmp = NamedTempFile::with_suffix(".txt")?;
        std::fs::write(tmp.path(), patterns.join("\n"))?;
        Ok(tmp)
    }

    /// Parse `split <args...>` and run it.
    fn run_split(args: &[&str]) -> Result<()> {
        let cmd = crate::cli::SplitCommand::try_parse_from(
            std::iter::once("split").chain(args.iter().copied()),
        )?;
        super::run(&cmd)
    }

    fn path(p: &std::path::Path) -> &str {
        p.to_str().unwrap()
    }

    /// Sum binseq record counts across every file in `dir` with `extension`.
    fn count_all_in_dir(dir: &std::path::Path, extension: &str) -> Result<usize> {
        let ext = extension.trim_start_matches('.');
        std::fs::read_dir(dir)?
            .filter_map(std::result::Result::ok)
            .filter(|e| {
                e.path()
                    .extension()
                    .and_then(|x| x.to_str())
                    .is_some_and(|x| x == ext)
            })
            .map(|e| count_binseq(&e.path()))
            .try_fold(0usize, |acc, r| r.map(|n| acc + n))
    }

    /// The total records across all split output files must equal the input count,
    /// regardless of how records distribute across patterns.
    #[test]
    fn test_split_record_conservation() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let in_tmp = write_fastx().call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            // Two patterns: no matter how records split, the total must equal input.
            let pat_file = write_patterns(&["AAAA", "CCCC"])?;
            let out_dir = tempfile::tempdir()?;

            run_split(&[
                path(bq_tmp.path()),
                "--file",
                path(pat_file.path()),
                "--basepath",
                path(out_dir.path()),
                "--min-records",
                "0",
                "--quiet",
            ])?; // keep empty files so everything is counted

            let total = count_all_in_dir(out_dir.path(), mode.extension())?;
            assert_eq!(
                total, DEFAULT_NUM_RECORDS,
                "split total count wrong for {mode:?}"
            );
        }
        Ok(())
    }

    /// `--span` used to be accepted but ignored.
    #[test]
    fn test_split_respects_span() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let pat_file = write_patterns(&["AAAA", "CCCC"])?;
        let out_dir = tempfile::tempdir()?;
        run_split(&[
            path(bq_tmp.path()),
            "--file",
            path(pat_file.path()),
            "--basepath",
            path(out_dir.path()),
            "--min-records",
            "0",
            "--span",
            "10..30",
            "--quiet",
        ])?;
        assert_eq!(count_all_in_dir(out_dir.path(), ".cbq")?, 20);
        Ok(())
    }

    /// An alias equal to the unmatched basename would open the same file twice.
    #[test]
    fn test_split_rejects_alias_colliding_with_unmatched() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let pat_file = NamedTempFile::with_suffix(".fa")?;
        std::fs::write(pat_file.path(), ">unmatched\nAAAA\n")?;
        let out_dir = tempfile::tempdir()?;
        assert!(run_split(&[
            path(bq_tmp.path()),
            "--file",
            path(pat_file.path()),
            "--basepath",
            path(out_dir.path()),
            "--quiet"
        ])
        .is_err());

        run_split(&[
            path(bq_tmp.path()),
            "--file",
            path(pat_file.path()),
            "--basepath",
            path(out_dir.path()),
            "--unmatched-basename",
            "rest",
            "--quiet",
        ])?;
        Ok(())
    }

    /// Using --skip-unmatched: no unmatched file should be created.
    /// With a universal pattern ("A"), every record should match, so the
    /// matched file contains all records.
    #[test]
    fn test_split_skip_unmatched() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        // "A" is a single-character pattern that matches any sequence containing
        // an 'A'. After BQ encoding (Ns replaced), all 100-base sequences will
        // contain at least one A with overwhelming probability.
        let pat_file = write_patterns(&["A"])?;
        let out_dir = tempfile::tempdir()?;

        run_split(&[
            path(bq_tmp.path()),
            "--file",
            path(pat_file.path()),
            "--basepath",
            path(out_dir.path()),
            "--skip-unmatched",
            "--quiet",
        ])?;

        // Only the "A.cbq" file should exist; verify its count.
        let matched_path = out_dir.path().join("A.cbq");
        assert!(matched_path.exists(), "expected A.cbq in output dir");
        assert_eq!(count_binseq(&matched_path)?, DEFAULT_NUM_RECORDS);

        // With --skip-unmatched and a universal pattern, nothing else should be there.
        let file_count = std::fs::read_dir(out_dir.path())?.count();
        assert_eq!(file_count, 1, "expected exactly one output file");

        Ok(())
    }

    /// Named patterns (FASTA headers) become output file aliases.
    #[test]
    fn test_split_named_patterns() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        // FASTA pattern: header becomes the alias → output file name
        let pat_file = {
            let tmp = NamedTempFile::with_suffix(".fasta")?;
            let mut f = std::fs::File::create(tmp.path())?;
            writeln!(f, ">universal_pattern")?;
            writeln!(f, "A")?;
            tmp
        };
        let out_dir = tempfile::tempdir()?;

        run_split(&[
            path(bq_tmp.path()),
            "--file",
            path(pat_file.path()),
            "--basepath",
            path(out_dir.path()),
            "--skip-unmatched",
            "--quiet",
        ])?;

        // The alias "universal_pattern" → "universal_pattern.cbq"
        let matched_path = out_dir.path().join("universal_pattern.cbq");
        assert!(matched_path.exists(), "expected universal_pattern.cbq");
        assert_eq!(count_binseq(&matched_path)?, DEFAULT_NUM_RECORDS);

        Ok(())
    }

    /// `--rc` should reverse complement file patterns before matching: a
    /// pattern file containing the reverse complement of a known substring
    /// should match reads containing that substring directly.
    #[test]
    fn test_split_rc_matches_reverse_complement() -> Result<()> {
        // "GATTACA" is not a palindrome; its reverse complement is "TGTAATC".
        let seq = "ACGTACGTGATTACAACGTACGT";
        let in_tmp = NamedTempFile::with_suffix(".fastq")?;
        {
            let mut f = std::fs::File::create(in_tmp.path())?;
            writeln!(f, "@read1")?;
            writeln!(f, "{seq}")?;
            writeln!(f, "+")?;
            writeln!(f, "{}", "I".repeat(seq.len()))?;
        }
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let pat_file = write_patterns(&["TGTAATC"])?;
        let out_dir = tempfile::tempdir()?;

        run_split(&[
            path(bq_tmp.path()),
            "--file",
            path(pat_file.path()),
            "--basepath",
            path(out_dir.path()),
            "--skip-unmatched",
            "--quiet",
            "--rc",
        ])?;

        // The pattern is reverse complemented to "GATTACA" before matching,
        // and the output alias reflects the RC'd sequence.
        let matched_path = out_dir.path().join("GATTACA.cbq");
        assert!(matched_path.exists(), "expected GATTACA.cbq after --rc");
        assert_eq!(count_binseq(&matched_path)?, 1);

        Ok(())
    }

    /// `--rc` must reject regex patterns, since reverse complementing a
    /// regex is undefined.
    #[test]
    fn test_split_rc_rejects_regex_pattern() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let pat_file = write_patterns(&["AC.GT"])?;
        let out_dir = tempfile::tempdir()?;

        assert!(
            run_split(&[
                path(bq_tmp.path()),
                "--file",
                path(pat_file.path()),
                "--basepath",
                path(out_dir.path()),
                "--quiet",
                "--rc",
            ])
            .is_err(),
            "--rc should reject regex patterns"
        );

        Ok(())
    }
}
