mod processor;

use anyhow::{bail, Result};
use binseq::BinseqReader;
use clap::ValueEnum;
use log::warn;
use serde::Serialize;

use crate::cli::{Mate, VerifyCommand, VerifyOptions};
use processor::{FieldMask, VerifyProcessor};

/// Whether `reader`'s underlying file actually stores per-record headers.
///
/// BQ never supports headers at all. For VBQ/CBQ files that don't store
/// them, [`BinseqRecord::sheader`](binseq::BinseqRecord::sheader)/`xheader`
/// still return *something* - a string synthesized from the record's
/// position in the file - for use by commands like `decode` that need a
/// name to print. That fallback must never be hashed by `verify`: it isn't
/// real header data, and it would leak record order into a checksum that's
/// supposed to be order-independent.
fn reader_has_headers(reader: &BinseqReader) -> bool {
    match reader {
        BinseqReader::Bq(_) => false,
        BinseqReader::Vbq(reader) => reader.header().headers,
        BinseqReader::Cbq(reader) => reader.header().has_headers(),
    }
}

fn field_mask(opts: &VerifyOptions) -> Result<FieldMask> {
    let fields = FieldMask {
        seq: !opts.skip_seq,
        qual: !opts.skip_qual,
        headers: !opts.skip_headers,
        flags: !opts.skip_flags,
    };
    if !(fields.seq || fields.qual || fields.headers || fields.flags) {
        bail!(
            "At least one field must be included in the checksum \
             (seq, qual, headers, flags cannot all be skipped)"
        );
    }
    Ok(fields)
}

fn field_labels(f: FieldMask) -> Vec<&'static str> {
    [
        ("seq", f.seq),
        ("qual", f.qual),
        ("headers", f.headers),
        ("flags", f.flags),
    ]
    .into_iter()
    .filter(|&(_, on)| on)
    .map(|(name, _)| name)
    .collect()
}

/// Runs the checksum computation without printing, so it can be reused by tests.
fn compute(args: &VerifyCommand) -> Result<VerifyReport> {
    let mut fields = field_mask(&args.opts)?;

    let reader = BinseqReader::new(args.input.path())?;
    if args.opts.mate == Mate::Two && !reader.is_paired() {
        bail!(
            "`--mate/-M 2` was requested but `{}` is single-channel (no extended/mate-2 \
             sequence); the checksum would be computed over no fields",
            args.input.path()
        );
    }

    let (stored_qual, stored_flags) = match &reader {
        BinseqReader::Bq(_) => (false, false),
        BinseqReader::Vbq(r) => (r.header().qual, r.header().flags),
        BinseqReader::Cbq(r) => (r.header().has_qualities(), r.header().has_flags()),
    };
    if fields.headers && !reader_has_headers(&reader) {
        warn!(
            "`{}` has no header data; excluding headers from the checksum",
            args.input.path()
        );
    }
    fields.headers &= reader_has_headers(&reader);
    fields.qual &= stored_qual;
    fields.flags &= stored_flags;
    if !(fields.seq || fields.qual || fields.headers || fields.flags) {
        bail!(
            "`{}` stores none of the selected fields; the checksum would be computed over no fields",
            args.input.path()
        );
    }

    let processor = VerifyProcessor::new(fields, args.opts.mate);

    let range = args.input.range(reader.num_records()?)?;
    reader.process_parallel_range(processor.clone(), args.opts.threads, range)?;

    Ok(VerifyReport {
        path: args.input.path().to_string(),
        algorithm: "xxh3-64/wrapping-sum",
        fields: field_labels(fields),
        mate: args
            .opts
            .mate
            .to_possible_value()
            .unwrap()
            .get_name()
            .to_string(),
        num_records: processor.num_records(),
        checksum: format!("{:016x}", processor.checksum()),
    })
}

#[derive(Serialize)]
struct VerifyReport {
    path: String,
    algorithm: &'static str,
    fields: Vec<&'static str>,
    mate: String,
    num_records: usize,
    checksum: String,
}

pub fn run(args: &VerifyCommand) -> Result<()> {
    let report = compute(args)?;

    if args.opts.json {
        println!("{}", serde_json::to_string_pretty(&report)?);
    } else {
        println!(
            "{}\t{}\t{}",
            report.checksum, report.num_records, report.path
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use clap::Parser;
    use itertools::iproduct;
    use tempfile::NamedTempFile;

    use crate::cli::BinseqMode;
    use crate::testutils::write_fastx;

    fn encode(in_path: &std::path::Path, out_path: &std::path::Path, extra: &[&str]) -> Result<()> {
        let cmd = crate::cli::EncodeCommand::try_parse_from(
            [
                "encode",
                in_path.to_str().unwrap(),
                "-o",
                out_path.to_str().unwrap(),
            ]
            .into_iter()
            .chain(extra.iter().copied()),
        )?;
        crate::commands::encode::run(&cmd)
    }

    fn checksum(path: &std::path::Path, extra: &[&str]) -> Result<u64> {
        let mut cmd_args = vec!["verify".to_string(), path.to_str().unwrap().to_string()];
        cmd_args.extend(extra.iter().map(std::string::ToString::to_string));
        let cmd = crate::cli::VerifyCommand::try_parse_from(cmd_args)?;
        Ok(u64::from_str_radix(&super::compute(&cmd)?.checksum, 16)?)
    }

    /// Re-encoding the same input twice (independent parallel runs, so record
    /// order between the two outputs is not guaranteed) must produce the
    /// same checksum.
    #[test]
    fn test_verify_stable_across_independent_encodes() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let in_tmp = write_fastx().call()?;

            let bq_a = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_a.path(), &[])?;
            let bq_b = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_b.path(), &[])?;

            let checksum_a = checksum(bq_a.path(), &[])?;
            let checksum_b = checksum(bq_b.path(), &[])?;
            assert_eq!(
                checksum_a, checksum_b,
                "checksum differed across independent encodes for {mode:?}"
            );
        }
        Ok(())
    }

    /// A corrupted byte in the encoded payload must change the checksum
    /// (when the file still parses after the corruption).
    #[test]
    fn test_verify_detects_content_change() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path(), &[])?;
        let original = checksum(bq_tmp.path(), &[])?;

        let mut bytes = std::fs::read(bq_tmp.path())?;
        let mid = bytes.len() / 2;
        bytes[mid] ^= 0xFF;
        std::fs::write(bq_tmp.path(), &bytes)?;

        if let Ok(changed) = checksum(bq_tmp.path(), &[]) {
            assert_ne!(original, changed, "bit flip was not detected");
        }
        Ok(())
    }

    /// `--skip-*` flags must actually change which fields feed the checksum.
    #[test]
    fn test_verify_skip_flags_change_checksum() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path(), &[])?;

        let full = checksum(bq_tmp.path(), &[])?;
        let skip_headers = checksum(bq_tmp.path(), &["--skip-headers"])?;
        let skip_qual = checksum(bq_tmp.path(), &["--skip-qual"])?;

        assert_ne!(full, skip_headers);
        assert_ne!(full, skip_qual);
        assert_ne!(skip_headers, skip_qual);
        Ok(())
    }

    /// `bqtools encode` never writes per-record flags, so `--skip-flags` has
    /// nothing to exclude on those files - it must be a no-op, not just
    /// "different from the other skip flags".
    #[test]
    fn test_verify_skip_flags_is_noop_without_flag_data() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path(), &[])?;

        let full = checksum(bq_tmp.path(), &[])?;
        let skip_flags = checksum(bq_tmp.path(), &["--skip-flags"])?;

        assert_eq!(full, skip_flags);
        Ok(())
    }

    /// A file encoded with `--skip-headers` (no header data at all) must
    /// hash the same whether or not `verify --skip-headers` is passed.
    #[test]
    fn test_verify_skip_headers_is_noop_without_header_data() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path(), &["--skip-headers"])?;

        let full = checksum(bq_tmp.path(), &[])?;
        let skip_headers = checksum(bq_tmp.path(), &["--skip-headers"])?;

        assert_eq!(full, skip_headers);
        Ok(())
    }

    /// Regression test: for a headerless file, `BinseqRecord::sheader()`
    /// falls back to a string synthesized from the record's position, since
    /// it's also used by commands like `decode` to print a name when none is
    /// stored. Hashing that fallback would leak record order into the
    /// checksum - defeating the entire point of `verify` - for any file
    /// large/parallel enough that record order isn't already incidentally
    /// stable across encodes. This must not happen: independent encodes of
    /// the same headerless input must produce the same checksum.
    #[test]
    fn test_verify_stable_across_independent_encodes_without_headers() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            // `include_n(false)`: bq/vbq redraw a random replacement for `N`
            // on every encode, which would (correctly) change the checksum
            // for an unrelated reason and mask the invariant under test.
            let in_tmp = write_fastx().nrec(20_000).include_n(false).call()?;

            let bq_a = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_a.path(), &["--skip-headers"])?;
            let bq_b = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_b.path(), &["--skip-headers"])?;

            let checksum_a = checksum(bq_a.path(), &[])?;
            let checksum_b = checksum(bq_b.path(), &[])?;
            assert_eq!(
                checksum_a, checksum_b,
                "headerless checksum differed across independent encodes for {mode:?}"
            );
        }
        Ok(())
    }

    /// If a file has no header data and the user also skips seq/qual/flags,
    /// there's nothing left to hash once headers are excluded - this must
    /// hard-error rather than silently produce a checksum over no fields.
    #[test]
    fn test_verify_rejects_no_fields_left_after_dropping_headers() -> Result<()> {
        use crate::cli::FileFormat;

        let in_tmp = write_fastx().format(FileFormat::Fasta).call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path(), &["--skip-headers"])?;

        let err = checksum(
            bq_tmp.path(),
            &["--skip-seq", "--skip-qual", "--skip-flags"],
        )
        .unwrap_err();
        assert!(err.to_string().contains("stores none"));
        Ok(())
    }

    /// `.bq` stores only sequences, so skipping seq leaves nothing to hash.
    #[test]
    fn test_verify_rejects_no_stored_fields_left() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".bq")?;
        encode(in_tmp.path(), bq_tmp.path(), &[])?;

        let err = checksum(bq_tmp.path(), &["--skip-seq", "--skip-qual"]).unwrap_err();
        assert!(err.to_string().contains("stores none"));
        Ok(())
    }

    /// Skipping every field is rejected up front.
    #[test]
    fn test_verify_rejects_empty_field_selection() {
        let cmd = crate::cli::VerifyCommand::try_parse_from([
            "verify",
            "input.cbq",
            "--skip-seq",
            "--skip-qual",
            "--skip-headers",
            "--skip-flags",
        ])
        .unwrap();
        assert!(super::field_mask(&cmd.opts).is_err());
    }

    /// Paired files restricted to a single mate must differ from the
    /// checksum over both mates.
    #[test]
    fn test_verify_mate_selection_changes_checksum() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let r1 = write_fastx().call()?;
            let r2 = write_fastx().call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            let cmd = crate::cli::EncodeCommand::try_parse_from([
                "encode",
                r1.path().to_str().unwrap(),
                r2.path().to_str().unwrap(),
                "-o",
                bq_tmp.path().to_str().unwrap(),
            ])?;
            crate::commands::encode::run(&cmd)?;

            let both = checksum(bq_tmp.path(), &[])?;
            let mate1 = checksum(bq_tmp.path(), &["-M", "1"])?;
            let mate2 = checksum(bq_tmp.path(), &["-M", "2"])?;

            assert_ne!(both, mate1, "mode={mode:?}");
            assert_ne!(both, mate2, "mode={mode:?}");
            assert_ne!(mate1, mate2, "mode={mode:?}");
        }
        Ok(())
    }

    /// Requesting mate 2 on a single-channel (unpaired) file must hard-error
    /// rather than silently produce a checksum over no fields. Mate 1 and
    /// mate "both" are unaffected, since the primary channel always exists.
    #[test]
    fn test_verify_rejects_mate_two_on_single_channel_file() -> Result<()> {
        let in_tmp = write_fastx().call()?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path(), &[])?;

        let err = checksum(bq_tmp.path(), &["-M", "2"]).unwrap_err();
        assert!(err.to_string().contains("--mate/-M 2"));

        assert!(checksum(bq_tmp.path(), &["-M", "1"]).is_ok());
        assert!(checksum(bq_tmp.path(), &["-M", "both"]).is_ok());
        Ok(())
    }

    /// `verify::run` must not error for any mode or output option.
    #[test]
    fn test_verify_run_all_modes() -> Result<()> {
        for (mode, json_flag) in iproduct!(BinseqMode::enum_iter(), [&[][..], &["--json"]]) {
            let in_tmp = write_fastx().call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path(), &[])?;

            let mut args = vec![
                "verify".to_string(),
                bq_tmp.path().to_str().unwrap().to_string(),
            ];
            args.extend(json_flag.iter().map(std::string::ToString::to_string));
            let cmd = crate::cli::VerifyCommand::try_parse_from(args)?;
            super::run(&cmd)?;
        }
        Ok(())
    }
}
