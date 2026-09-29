mod engine;
mod filter;
mod pattern_count;
mod patterns;
mod range;

pub use engine::{Engine, Spans};
use filter::FilterProcessor;
use pattern_count::{PatternCountProcessor, PatternCounter};
pub use patterns::{Pattern, PatternCollection, PatternSets};
pub use range::SimpleRange;

use super::decode::build_writer;
use crate::{
    cli::{FileFormat, GrepCommand, Mate},
    commands::decode::SplitWriter,
};

use anyhow::{bail, Result};
use binseq::prelude::*;

fn load_patterns(args: &GrepCommand, paired: bool) -> Result<PatternSets> {
    let mut patterns = PatternSets {
        pat1: args.grep.patterns_m1()?,
        pat2: args.grep.patterns_m2()?,
        pat: args.grep.patterns()?,
    };
    // `--mate` is meaningless on single-end files (and `-m 2` would route every
    // pattern to the empty extended sequence), and extended-only patterns can
    // never match, so single-end input never carries an extended pattern set.
    if paired {
        patterns.redistribute(args.output.mate)?;
    } else if !patterns.pat2.is_empty() {
        bail!("-R/--xfile patterns require paired input");
    }
    if args.grep.rc {
        patterns.reverse_complement()?;
    }
    Ok(patterns)
}

/// `and_logic` (all patterns must hit) only applies to plain grep, not `-P`.
fn build_engine(args: &GrepCommand, patterns: &PatternSets, and_logic: bool) -> Result<Engine> {
    #[cfg(feature = "fuzzy")]
    if args.grep.fuzzy_args.fuzzy {
        return Engine::fuzzy(
            patterns,
            args.grep.fuzzy_args.distance,
            args.grep.fuzzy_args.inexact,
            args.grep.fuzzy_args.max_n_frac,
        );
    }
    if patterns.use_fixed(args.grep.fixed) {
        if and_logic {
            // ponytail: the regex crate's literal search beats aho-corasick's
            // overlapping scan for the few patterns AND is used with; switch
            // back if AND over hundreds of fixed patterns matters
            Engine::regex(&patterns.escaped()?)
        } else {
            Engine::aho_corasick(patterns, args.grep.no_dfa)
        }
    } else {
        Engine::regex(patterns)
    }
}

fn run_pattern_count(args: &GrepCommand, reader: BinseqReader) -> Result<()> {
    let patterns = load_patterns(args, reader.is_paired())?;
    let counter = PatternCounter::new(
        build_engine(args, &patterns, false)?,
        &patterns,
        args.grep.invert,
    );
    let proc = PatternCountProcessor::new(
        counter,
        args.grep.range.unwrap_or_default(),
        args.grep.header,
    );
    if let Some(span) = args.input.span {
        let num_records = reader.num_records()?;
        reader.process_parallel_range(
            proc.clone(),
            args.output.threads(),
            span.get_range(num_records)?,
        )?;
    } else {
        reader.process_parallel(proc.clone(), args.output.threads())?;
    }
    proc.pprint_pattern_counts()?;
    Ok(())
}

fn run_grep(
    args: &GrepCommand,
    reader: BinseqReader,
    writer: SplitWriter,
    format: FileFormat,
    mate: Option<Mate>,
) -> Result<()> {
    let count = args.grep.count || args.grep.frac;
    let patterns = load_patterns(args, reader.is_paired())?;
    // AND vs OR logic is only meaningful when combining 2+ patterns
    let and_logic = args.grep.and_logic() && patterns.len() > 1;
    let proc = FilterProcessor::new(
        build_engine(args, &patterns, and_logic)?,
        and_logic,
        args.grep.invert,
        count,
        args.grep.frac,
        args.grep.range.unwrap_or_default(),
        args.grep.header,
        writer,
        format,
        mate,
        args.should_color(),
    );

    if let Some(span) = args.input.span {
        let num_records = reader.num_records()?;
        reader.process_parallel_range(
            proc.clone(),
            args.output.threads(),
            span.get_range(num_records)?,
        )?;
    } else {
        reader.process_parallel(proc.clone(), args.output.threads())?;
    }
    if count {
        proc.pprint_counts();
    }

    Ok(())
}

pub fn run(args: &GrepCommand) -> Result<()> {
    args.grep.validate()?;
    let reader = BinseqReader::new(args.input.path())?;
    if args.grep.pattern_count {
        return run_pattern_count(args, reader);
    }
    let format = args.output.format()?;
    let writer = build_writer(&args.output, format, reader.is_paired())?;
    let mate = if reader.is_paired() {
        Some(args.output.mate)
    } else {
        None
    };
    run_grep(args, reader, writer, format, mate)
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use clap::Parser;
    use itertools::iproduct;
    use tempfile::NamedTempFile;

    use crate::cli::{BinseqMode, FileFormat};
    use crate::testutils::{count_fastx_records, write_fastx, DEFAULT_NUM_RECORDS};

    fn encode(in_path: &std::path::Path, out_path: &std::path::Path) -> Result<()> {
        let cmd = crate::cli::EncodeCommand::try_parse_from([
            "encode",
            in_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])?;
        crate::commands::encode::run(&cmd)
    }

    fn grep_count(bq_path: &std::path::Path, pattern: &str, invert: bool) -> Result<usize> {
        // Build a grep command that writes to a temp file and count the result.
        let out_tmp = NamedTempFile::with_suffix(".fastq")?;
        let mut args = vec![
            "grep",
            bq_path.to_str().unwrap(),
            pattern,
            "-o",
            out_tmp.path().to_str().unwrap(),
        ];
        if invert {
            args.push("-v");
        }
        let cmd = crate::cli::GrepCommand::try_parse_from(args)?;
        super::run(&cmd)?;
        count_fastx_records(out_tmp.path())
    }

    /// grep returns a count ≤ total records and > 0 for a short common pattern.
    #[test]
    fn test_grep_basic_count() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let in_tmp = write_fastx().call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let count = grep_count(bq_tmp.path(), "A", false)?;
            assert!(count > 0, "grep count should be > 0 for {mode:?}");
            assert!(
                count <= DEFAULT_NUM_RECORDS,
                "grep count {count} exceeds total for {mode:?}"
            );
        }
        Ok(())
    }

    /// A single fixed-string pattern under the default AND logic must not
    /// panic (Aho-Corasick doesn't support AND) and must match the count
    /// produced with explicit OR logic, since AND/OR are equivalent with
    /// only one pattern.
    #[test]
    fn test_grep_single_fixed_pattern_with_default_and_logic() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let in_tmp = write_fastx().call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let and_count = grep_count(bq_tmp.path(), "AAAA", false)?;

            let out_tmp = NamedTempFile::with_suffix(".fastq")?;
            let cmd = crate::cli::GrepCommand::try_parse_from([
                "grep",
                bq_tmp.path().to_str().unwrap(),
                "AAAA",
                "-o",
                out_tmp.path().to_str().unwrap(),
                "--or-logic",
            ])?;
            super::run(&cmd)?;
            let or_count = count_fastx_records(out_tmp.path())?;

            assert_eq!(
                and_count, or_count,
                "AND and OR logic should match for a single pattern, mode={mode:?}"
            );
        }
        Ok(())
    }

    /// forward matches + inverted matches must equal the total record count exactly.
    #[test]
    fn test_grep_invert_complementary() -> Result<()> {
        for (mode, fmt) in iproduct!(BinseqMode::enum_iter(), FileFormat::fastx_iter()) {
            let in_tmp = write_fastx().format(fmt).call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let fwd = grep_count(bq_tmp.path(), "AAAA", false)?;
            let inv = grep_count(bq_tmp.path(), "AAAA", true)?;
            assert_eq!(
                fwd + inv,
                DEFAULT_NUM_RECORDS,
                "fwd({fwd}) + inv({inv}) != {DEFAULT_NUM_RECORDS} for {mode:?} {fmt:?}"
            );
        }
        Ok(())
    }

    /// grep writes matching records to a file across all (mode, format) combinations.
    #[test]
    fn test_grep_all_modes_and_formats() -> Result<()> {
        for (mode, fmt) in iproduct!(BinseqMode::enum_iter(), FileFormat::fastx_iter()) {
            let in_tmp = write_fastx().format(fmt).call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let out_tmp = NamedTempFile::with_suffix(fmt.fastx_suffix())?;
            let cmd = crate::cli::GrepCommand::try_parse_from([
                "grep",
                bq_tmp.path().to_str().unwrap(),
                "A",
                "-o",
                out_tmp.path().to_str().unwrap(),
            ])?;
            super::run(&cmd)
                .map_err(|e| anyhow::anyhow!("grep failed for {mode:?} {fmt:?}: {e}"))?;

            let count = count_fastx_records(out_tmp.path())?;
            assert!(
                count <= DEFAULT_NUM_RECORDS,
                "grep output count {count} > total for {mode:?} {fmt:?}"
            );
        }
        Ok(())
    }

    /// `--rc` should reverse complement the pattern before matching: searching
    /// for the reverse complement of a known substring should match exactly
    /// as if the substring itself had been used directly (without --rc).
    #[test]
    fn test_grep_rc_matches_reverse_complement() -> Result<()> {
        use std::io::Write as _;

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

        // Without --rc, searching for the RC pattern directly should not match.
        let direct = grep_count(bq_tmp.path(), "TGTAATC", false)?;
        assert_eq!(direct, 0, "TGTAATC should not appear literally in the read");

        // With --rc, the pattern is reverse complemented back to "GATTACA",
        // which does appear in the read.
        let out_tmp = NamedTempFile::with_suffix(".fastq")?;
        let cmd = crate::cli::GrepCommand::try_parse_from([
            "grep",
            bq_tmp.path().to_str().unwrap(),
            "TGTAATC",
            "-o",
            out_tmp.path().to_str().unwrap(),
            "--rc",
        ])?;
        super::run(&cmd)?;
        let rc_count = count_fastx_records(out_tmp.path())?;
        assert_eq!(rc_count, 1, "--rc should match the reverse complement");
        Ok(())
    }

    /// `--rc` must reject regex patterns, since reverse complementing a
    /// regex is undefined.
    #[test]
    fn test_grep_rc_rejects_regex_pattern() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let out_tmp = NamedTempFile::with_suffix(".fastq")?;
        let cmd = crate::cli::GrepCommand::try_parse_from([
            "grep",
            bq_tmp.path().to_str().unwrap(),
            "AC.GT",
            "-o",
            out_tmp.path().to_str().unwrap(),
            "--rc",
        ])?;
        assert!(
            super::run(&cmd).is_err(),
            "--rc should reject regex patterns"
        );
        Ok(())
    }

    /// `--header` should match against the record header (`seq.{idx}`, per
    /// `write_fastx`) rather than the sequence. The literal "seq.1" can never
    /// appear in a random ACGT(N) sequence, so it isolates header matching:
    /// it should hit headers seq.1 and seq.10-seq.19 (11 records out of 100)
    /// when searching headers, and 0 records when searching sequences.
    #[test]
    fn test_grep_header_matches_header_not_sequence() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let header_out = NamedTempFile::with_suffix(".fastq")?;
        let header_cmd = crate::cli::GrepCommand::try_parse_from([
            "grep",
            bq_tmp.path().to_str().unwrap(),
            "seq.1",
            "-o",
            header_out.path().to_str().unwrap(),
            "--header",
            "-x",
        ])?;
        super::run(&header_cmd)?;
        let header_count = count_fastx_records(header_out.path())?;
        assert_eq!(
            header_count, 11,
            "expected seq.1 and seq.10-seq.19 to match on header"
        );

        let seq_out = NamedTempFile::with_suffix(".fastq")?;
        let seq_cmd = crate::cli::GrepCommand::try_parse_from([
            "grep",
            bq_tmp.path().to_str().unwrap(),
            "seq.1",
            "-o",
            seq_out.path().to_str().unwrap(),
            "-x",
        ])?;
        super::run(&seq_cmd)?;
        let seq_count = count_fastx_records(seq_out.path())?;
        assert_eq!(
            seq_count, 0,
            "'seq.1' cannot appear in an ACGTN sequence, so sequence-mode grep should find nothing"
        );
        Ok(())
    }

    /// Writes a minimal paired FASTQ record set with the given headers and
    /// arbitrary-but-distinct ACGT sequences.
    fn write_paired_fastq(path: &std::path::Path, headers: &[&str]) -> Result<()> {
        use std::io::Write as _;
        let bases = *b"ACGT";
        let mut f = std::fs::File::create(path)?;
        for (idx, header) in headers.iter().enumerate() {
            let seq: String = (0..12).map(|i| bases[(idx + i) % 4] as char).collect();
            writeln!(f, "@{header}")?;
            writeln!(f, "{seq}")?;
            writeln!(f, "+")?;
            writeln!(f, "{}", "I".repeat(seq.len()))?;
        }
        Ok(())
    }

    /// `-r`/`-R`/positional `--header` patterns must match the primary
    /// header, extended header, and either header respectively — and must
    /// NOT cross-match the other mate's header.
    #[test]
    fn test_grep_header_respects_primary_extended_either() -> Result<()> {
        let r1_tmp = NamedTempFile::with_suffix(".fastq")?;
        let r2_tmp = NamedTempFile::with_suffix(".fastq")?;
        write_paired_fastq(r1_tmp.path(), &["primary_A", "primary_B", "primary_C"])?;
        write_paired_fastq(r2_tmp.path(), &["extended_X", "extended_Y", "extended_Z"])?;

        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        let encode_cmd = crate::cli::EncodeCommand::try_parse_from([
            "encode",
            r1_tmp.path().to_str().unwrap(),
            r2_tmp.path().to_str().unwrap(),
            "-o",
            bq_tmp.path().to_str().unwrap(),
        ])?;
        crate::commands::encode::run(&encode_cmd)?;

        // Counts matched pairs. `--mate` (unset here, defaults to `both`)
        // also controls which mate's patterns are searched (see
        // `redistribute_patterns`), so it must stay unset to test `-r`/`-R`
        // independently; instead, matched pairs are written interleaved
        // (both mates per match) and the record count is halved.
        let run = |extra: &[&str]| -> Result<usize> {
            let out_tmp = NamedTempFile::with_suffix(".fastq")?;
            let mut args = vec![
                "grep",
                bq_tmp.path().to_str().unwrap(),
                "-o",
                out_tmp.path().to_str().unwrap(),
                "--header",
                "-x",
            ];
            args.extend_from_slice(extra);
            let cmd = crate::cli::GrepCommand::try_parse_from(args)?;
            super::run(&cmd)?;
            let records = count_fastx_records(out_tmp.path())?;
            assert_eq!(records % 2, 0, "expected interleaved R1+R2 pairs");
            Ok(records / 2)
        };

        // `-r` searches the primary header only.
        assert_eq!(
            run(&["-r", "primary_B"])?,
            1,
            "-r should match the primary header it names"
        );
        assert_eq!(
            run(&["-r", "extended_X"])?,
            0,
            "-r must not match the extended header"
        );

        // `-R` searches the extended header only.
        assert_eq!(
            run(&["-R", "extended_Z"])?,
            1,
            "-R should match the extended header it names"
        );
        assert_eq!(
            run(&["-R", "primary_A"])?,
            0,
            "-R must not match the primary header"
        );

        // A positional (either) pattern matches against whichever header
        // contains it.
        assert_eq!(
            run(&["primary_C"])?,
            1,
            "positional pattern should match via the primary header"
        );
        assert_eq!(
            run(&["extended_Y"])?,
            1,
            "positional pattern should match via the extended header"
        );

        Ok(())
    }

    /// `--header` must reject `--rc`, since reverse complementing header text
    /// is undefined.
    #[test]
    fn test_grep_header_rejects_rc() {
        let result = crate::cli::GrepCommand::try_parse_from([
            "grep", "input.bq", "seq.1", "--header", "--rc",
        ]);
        assert!(result.is_err(), "--header should conflict with --rc");
    }

    /// `--header` must reject `--range`, since a coordinate range is
    /// meaningless against header text.
    #[test]
    fn test_grep_header_rejects_range() {
        let result = crate::cli::GrepCommand::try_parse_from([
            "grep", "input.bq", "seq.1", "--header", "--range", "0..10",
        ]);
        assert!(result.is_err(), "--header should conflict with --range");
    }

    /// `-m 2` on a single-end file used to move every pattern to the empty
    /// extended sequence, so nothing matched.
    #[test]
    fn test_grep_single_end_ignores_mate() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        let out_tmp = NamedTempFile::with_suffix(".fastq")?;
        let cmd = crate::cli::GrepCommand::try_parse_from([
            "grep",
            bq_tmp.path().to_str().unwrap(),
            "A",
            "-m",
            "2",
            "-o",
            out_tmp.path().to_str().unwrap(),
        ])?;
        super::run(&cmd)?;
        assert!(count_fastx_records(out_tmp.path())? > 0);
        Ok(())
    }

    /// Count modes never write records, so `-o/-p` are rejected up front.
    #[test]
    fn test_grep_count_modes_reject_output() {
        for flag in ["-C", "-F", "-P"] {
            assert!(
                crate::cli::GrepCommand::try_parse_from(["grep", "x.cbq", "A", flag, "-o", "o.fq"])
                    .is_err(),
                "{flag} should conflict with -o"
            );
        }
    }

    #[cfg(feature = "fuzzy")]
    #[test]
    fn test_grep_fuzzy_options_require_fuzzy() {
        for extra in [&["-k", "2"][..], &["-i"], &["--max-n-frac", "0.5"]] {
            let mut argv = vec!["grep", "x.cbq", "ACGT"];
            argv.extend_from_slice(extra);
            assert!(
                crate::cli::GrepCommand::try_parse_from(&argv).is_err(),
                "{extra:?} should require -z"
            );
            argv.push("-z");
            assert!(crate::cli::GrepCommand::try_parse_from(&argv).is_ok());
        }
    }
}
