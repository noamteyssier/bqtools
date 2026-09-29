use std::{fs::File, io::Write};

use anyhow::{anyhow, bail, ensure, Result};
use binseq::{bq, cbq, vbq, BinseqReader, BinseqWriter, BinseqWriterBuilder, ParallelReader};
use memmap2::MmapOptions;

use crate::{
    cli::{BinseqMode, CatCommand},
    commands::encode::processor::Encoder,
};

/// Returns the header of the first file, bailing if any other file differs.
fn same_header<H: PartialEq>(paths: &[String], get: impl Fn(&str) -> Result<H>) -> Result<H> {
    let mut it = paths.iter();
    let first = get(it.next().ok_or_else(|| anyhow!("No input files."))?)?;
    for path in it {
        if get(path)? != first {
            bail!("Inconsistent header found for path: {path}");
        }
    }
    Ok(first)
}

fn determine_mode(paths: &[String]) -> Result<BinseqMode> {
    let modes = paths
        .iter()
        .map(|path| {
            Ok(match BinseqReader::new(path)? {
                BinseqReader::Bq(_) => BinseqMode::Bq,
                BinseqReader::Vbq(_) => BinseqMode::Vbq,
                BinseqReader::Cbq(_) => BinseqMode::Cbq,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    ensure!(
        modes.windows(2).all(|w| w[0] == w[1]),
        "Inconsistent modes found, expecting the same BINSEQ mode for all input files."
    );
    modes
        .first()
        .copied()
        .ok_or_else(|| anyhow!("No input files."))
}

fn run_bq(args: CatCommand) -> Result<()> {
    let header = same_header(&args.input.input, |p| Ok(bq::MmapReader::new(p)?.header()))?;
    let mut out_handle = args.output.as_writer(BinseqMode::Bq)?;

    header.write_bytes(&mut out_handle)?;
    for path in args.input.input {
        let file = File::open(path)?;
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        out_handle.write_all(&mmap[bq::SIZE_HEADER..])?;
    }
    out_handle.flush()?;

    Ok(())
}

fn run_cat(args: CatCommand, writer: BinseqWriter<Box<dyn Write + Send>>) -> Result<()> {
    let mut processor = Encoder::new(writer)?;
    for path in args.input.input {
        let reader = BinseqReader::new(&path)?;
        reader.process_parallel(processor.clone(), args.output.threads())?;
    }
    processor.finish()?;
    Ok(())
}

pub fn run(args: CatCommand) -> Result<()> {
    let paths = &args.input.input;
    match determine_mode(paths)? {
        BinseqMode::Bq => run_bq(args),
        BinseqMode::Vbq => {
            let out = args.output.as_writer(BinseqMode::Vbq)?;
            let header = same_header(paths, |p| Ok(vbq::MmapReader::new(p)?.header()))?;
            let writer = BinseqWriterBuilder::from_vbq_header(header).build(out)?;
            run_cat(args, writer)
        }
        BinseqMode::Cbq => {
            let out = args.output.as_writer(BinseqMode::Cbq)?;
            let header = same_header(paths, |p| Ok(cbq::MmapReader::new(p)?.header()))?;
            let writer = BinseqWriterBuilder::from_cbq_header(header).build(out)?;
            run_cat(args, writer)
        }
    }
}

#[cfg(test)]
mod tests {
    use anyhow::Result;
    use clap::Parser;
    use itertools::iproduct;
    use tempfile::NamedTempFile;

    use crate::cli::{BinseqMode, FileFormat};
    use crate::testutils::{count_binseq, write_fastx, Compression, DEFAULT_NUM_RECORDS};

    fn encode(in_path: &std::path::Path, out_path: &std::path::Path) -> Result<()> {
        let cmd = crate::cli::EncodeCommand::try_parse_from([
            "encode",
            in_path.to_str().unwrap(),
            "-o",
            out_path.to_str().unwrap(),
        ])?;
        crate::commands::encode::run(&cmd)
    }

    fn cat(in_paths: &[&std::path::Path], out_path: &std::path::Path) -> Result<()> {
        let mut args = vec!["cat".to_string()];
        for p in in_paths {
            args.push(p.to_str().unwrap().to_string());
        }
        args.extend(["-o".to_string(), out_path.to_str().unwrap().to_string()]);
        let cmd = crate::cli::CatCommand::try_parse_from(args)?;
        super::run(cmd)
    }

    /// Concatenating two N-record files must produce 2*N records.
    #[test]
    fn test_cat_two_files() -> Result<()> {
        for (mode, fmt) in iproduct!(BinseqMode::enum_iter(), FileFormat::fastx_iter()) {
            let in1 = write_fastx().format(fmt).call()?;
            let in2 = write_fastx().format(fmt).call()?;
            let bq1 = NamedTempFile::with_suffix(mode.extension())?;
            let bq2 = NamedTempFile::with_suffix(mode.extension())?;
            encode(in1.path(), bq1.path())?;
            encode(in2.path(), bq2.path())?;

            let out = NamedTempFile::with_suffix(mode.extension())?;
            cat(&[bq1.path(), bq2.path()], out.path())?;

            assert_eq!(
                count_binseq(out.path())?,
                DEFAULT_NUM_RECORDS * 2,
                "cat 2-file count wrong for {mode:?} {fmt:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_cat_three_files() -> Result<()> {
        for mode in BinseqMode::enum_iter() {
            let parts: Vec<_> = (0..3)
                .map(|_| {
                    let in_tmp = write_fastx().call()?;
                    let bq = NamedTempFile::with_suffix(mode.extension())?;
                    encode(in_tmp.path(), bq.path())?;
                    Ok::<_, anyhow::Error>((in_tmp, bq))
                })
                .collect::<Result<_>>()?;

            let bq_paths: Vec<_> = parts.iter().map(|(_, bq)| bq.path()).collect();
            let out = NamedTempFile::with_suffix(mode.extension())?;
            cat(&bq_paths, out.path())?;

            assert_eq!(
                count_binseq(out.path())?,
                DEFAULT_NUM_RECORDS * 3,
                "cat 3-file count wrong for {mode:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn test_cat_compressed_inputs() -> Result<()> {
        for (mode, comp) in iproduct!(BinseqMode::enum_iter(), Compression::all()) {
            let in1 = write_fastx().comp(comp).call()?;
            let in2 = write_fastx().comp(comp).call()?;
            let bq1 = NamedTempFile::with_suffix(mode.extension())?;
            let bq2 = NamedTempFile::with_suffix(mode.extension())?;
            encode(in1.path(), bq1.path())?;
            encode(in2.path(), bq2.path())?;

            let out = NamedTempFile::with_suffix(mode.extension())?;
            cat(&[bq1.path(), bq2.path()], out.path())?;

            assert_eq!(
                count_binseq(out.path())?,
                DEFAULT_NUM_RECORDS * 2,
                "cat compressed count wrong for {mode:?} {comp:?}"
            );
        }
        Ok(())
    }
}
