mod processor;

use anyhow::Result;
use binseq::BinseqReader;
use log::{info, warn};

use crate::{
    cli::{Mate, RevcompCommand},
    commands::utils::builder_from_reader,
};
use processor::RevCompProcessor;

pub fn run(args: &RevcompCommand) -> Result<()> {
    let reader = BinseqReader::new(args.input.path())?;
    let mate = if !reader.is_paired() && args.mate != Mate::Both {
        warn!("Ignoring `--mate/-M` flag as only single channel found in file");
        Mate::Both
    } else {
        args.mate
    };

    let builder = builder_from_reader(&reader);
    let ohandle = args.output.as_writer(args.input.mode()?)?;
    let writer = builder.build(ohandle)?;
    let mut processor = RevCompProcessor::new(writer, mate)?;

    let num_records = reader.num_records()?;
    let range = match args.input.span {
        Some(span) => span.get_range(num_records)?,
        None => 0..num_records,
    };
    reader.process_parallel_range(processor.clone(), args.output.threads(), range)?;
    processor.finish()?;

    info!(
        "Wrote {} reverse complemented records",
        processor.get_global_record_count()
    );

    Ok(())
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use clap::Parser;
    use itertools::iproduct;
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

    fn decode_to_fasta(bq_path: &std::path::Path, out_path: &std::path::Path) -> Result<()> {
        let cmd = crate::cli::DecodeCommand::try_parse_from([
            "decode",
            bq_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
            "-f",
            "a",
        ])?;
        crate::commands::decode::run(&cmd)
    }

    fn revcomp(
        in_path: &std::path::Path,
        out_path: &std::path::Path,
        extra: &[&str],
    ) -> Result<()> {
        let mut cmd_args = vec![
            "revcomp".to_string(),
            in_path.to_str().unwrap().to_string(),
            "-o".to_string(),
            out_path.to_str().unwrap().to_string(),
        ];
        cmd_args.extend(extra.iter().map(std::string::ToString::to_string));
        let cmd = crate::cli::RevcompCommand::try_parse_from(cmd_args)?;
        super::run(&cmd)
    }

    /// Extracts just the sequence lines from a FASTA file, sorted, so
    /// comparisons are insensitive to reordering from parallel processing.
    fn sorted_sequences(path: &std::path::Path) -> Result<Vec<String>> {
        let content = std::fs::read_to_string(path)?;
        let mut seqs: Vec<String> = content
            .lines()
            .filter(|l| !l.starts_with('>'))
            .map(std::string::ToString::to_string)
            .collect();
        seqs.sort_unstable();
        Ok(seqs)
    }

    /// Round-tripping revcomp twice must recover the original sequences.
    #[test]
    fn test_revcomp_double_application_is_identity() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let in_tmp = write_fastx().call()?;
            let bq_tmp = NamedTempFile::with_suffix(mode.extension())?;
            encode(in_tmp.path(), bq_tmp.path())?;

            let rc_once = NamedTempFile::with_suffix(mode.extension())?;
            revcomp(bq_tmp.path(), rc_once.path(), &[])?;

            let rc_twice = NamedTempFile::with_suffix(mode.extension())?;
            revcomp(rc_once.path(), rc_twice.path(), &[])?;

            let original_fa = NamedTempFile::with_suffix(".fasta")?;
            decode_to_fasta(bq_tmp.path(), original_fa.path())?;
            let roundtrip_fa = NamedTempFile::with_suffix(".fasta")?;
            decode_to_fasta(rc_twice.path(), roundtrip_fa.path())?;

            assert_eq!(
                sorted_sequences(original_fa.path())?,
                sorted_sequences(roundtrip_fa.path())?,
                "double revcomp should be identity for {mode:?}"
            );

            assert_eq!(
                count_binseq(rc_once.path())?,
                DEFAULT_NUM_RECORDS,
                "revcomp record count wrong for {mode:?}"
            );
        }
        Ok(())
    }

    /// Reverse complementing a known sequence should produce the expected result.
    /// On single-end files `-M` is ignored, so `-M 2` must behave the same.
    #[test]
    fn test_revcomp_known_sequence() -> Result<()> {
        let seq = "ACGTACGTGATTACAACGTACGT";
        let expected: String = seq
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                'T' => 'A',
                other => other,
            })
            .collect();
        let in_tmp = NamedTempFile::with_suffix(".fastq")?;
        std::fs::write(
            in_tmp.path(),
            format!("@read1\n{seq}\n+\n{}\n", "I".repeat(seq.len())),
        )?;
        let bq_tmp = NamedTempFile::with_suffix(".cbq")?;
        encode(in_tmp.path(), bq_tmp.path())?;

        for extra in [&[][..], &["-M", "2"][..]] {
            let rc_tmp = NamedTempFile::with_suffix(".cbq")?;
            revcomp(bq_tmp.path(), rc_tmp.path(), extra)?;

            let out_fa = NamedTempFile::with_suffix(".fasta")?;
            decode_to_fasta(rc_tmp.path(), out_fa.path())?;
            let content = std::fs::read_to_string(out_fa.path())?;
            assert!(
                content.contains(&expected),
                "expected {expected} in output for args {extra:?}: {content}"
            );
        }
        Ok(())
    }

    /// With `-M 1`/`-M 2`, only the targeted mate should be reverse complemented.
    #[test]
    fn test_revcomp_paired_single_mate() -> Result<()> {
        for (mode, mate_flag) in iproduct!(BinseqMode::enum_iter(), ["1", "2"]) {
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

            let rc_tmp = NamedTempFile::with_suffix(mode.extension())?;
            revcomp(bq_tmp.path(), rc_tmp.path(), &["-M", mate_flag])?;

            assert_eq!(
                count_binseq(rc_tmp.path())?,
                DEFAULT_NUM_RECORDS,
                "revcomp single-mate record count wrong for {mode:?} mate={mate_flag}"
            );
        }
        Ok(())
    }
}
