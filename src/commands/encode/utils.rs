use std::{path::PathBuf, sync::LazyLock};

use anyhow::{bail, Context, Result};
use hashbrown::HashMap;
use log::warn;
use paraseq::{fastx, Record};
use regex::Regex;

#[cfg(feature = "htslib")]
use paraseq::rust_htslib::{self, bam::Read as BamRead};

use crate::types::BoxedReader;

/// Sequence lengths of (up to) the first `n` records of a reader.
///
/// The result is empty if the reader has no data, and shorter than `n` if it has fewer records.
fn first_lens(reader: &mut fastx::Reader<BoxedReader>, n: usize) -> Result<Vec<u32>> {
    let mut rset = reader.new_record_set_with_size(n);
    if !rset.fill(reader)? {
        return Ok(Vec::new());
    }
    let lens = rset
        .iter()
        .take(n)
        .map(|record| Ok(record?.seq().len() as u32))
        .collect::<Result<Vec<_>>>()?;
    reader.reload(&mut rset)?;
    Ok(lens)
}

pub fn get_sequence_len(reader: &mut fastx::Reader<BoxedReader>) -> Result<u32> {
    first_lens(reader, 1)?
        .first()
        .copied()
        .context("Input file is empty - cannot convert")
}

#[cfg(feature = "htslib")]
pub fn get_sequence_len_htslib(path: &str, paired: bool) -> Result<(u32, u32)> {
    let mut reader = rust_htslib::bam::Reader::from_path(path)?;
    let lens = reader
        .rc_records()
        .take(1 + usize::from(paired))
        .map(|res| res.map(|rec| rec.seq_len() as u32))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((
        lens.first().copied().unwrap_or(0),
        lens.get(1).copied().unwrap_or(0),
    ))
}

pub fn get_interleaved_sequence_len(reader: &mut fastx::Reader<BoxedReader>) -> Result<(u32, u32)> {
    match first_lens(reader, 2)?[..] {
        [slen, xlen] => Ok((slen, xlen)),
        [] => bail!("Input file (interleaved) is missing R2 - cannot convert"),
        _ => bail!("Input file is empty - cannot convert"),
    }
}

/// Pairs R1/R2 files from a list of file paths efficiently using a `HashMap`
/// Returns a vector of pairs, where each pair is [`R1_file`, `R2_file`]
pub fn pair_r1_r2_files(files: &[PathBuf]) -> Result<Vec<Vec<PathBuf>>> {
    let pair_regex = Regex::new(r"^(.+)_R([12])(_[^.]*)?\.(?:fastq|fq|fasta|fa)(?:\.gz|\.zst)?$")?;

    // HashMap to store files by their pairing key (base + suffix)
    let mut r1_files: HashMap<String, PathBuf> = HashMap::new();
    let mut r2_files: HashMap<String, PathBuf> = HashMap::new();

    // Single pass through files to categorize them
    for file in files {
        let file_str = file.to_str().unwrap();

        if let Some(caps) = pair_regex.captures(file_str) {
            let base = &caps[1];
            let read_num = &caps[2];
            let suffix = caps.get(3).map_or("", |m| m.as_str());

            // Create a unique key for pairing: base + suffix
            let pair_key = format!("{base}{suffix}");

            match read_num {
                "1" => {
                    r1_files.insert(pair_key, file.clone());
                }
                "2" => {
                    r2_files.insert(pair_key, file.clone());
                }
                _ => unreachable!(), // regex only matches 1 or 2
            }
        }
    }

    // Create pairs by finding matching keys
    let mut pairs = Vec::new();
    for (pair_key, r1_file) in &r1_files {
        if let Some(r2_file) = r2_files.get(pair_key) {
            pairs.push(vec![r1_file.to_owned(), r2_file.to_owned()]);
        } else {
            warn!("No R2 pair found for {} (skipping)", r1_file.display());
        }
    }

    // Check for orphaned R2 files
    for (pair_key, r2_file) in &r2_files {
        if !r1_files.contains_key(pair_key) {
            warn!("No R1 pair found for {} (skipping)", r2_file.display());
        }
    }

    // Sort pairs by R1 filename for consistent output
    pairs.sort_by(|a, b| a[0].cmp(&b[0]));

    Ok(pairs)
}

/// Generates a unique output filename based on the input file(s)
/// For single files: removes the original extension and replaces with new extension
/// For paired files: extracts the base name + suffix (everything except _R[12]) and adds new extension
pub fn generate_output_name(input_files: &[PathBuf], new_extension: &str) -> Result<String> {
    static PAIR: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"^(.+)_R[12](_[^.]*)?\.(?:fastq|fq|fasta|fa)(?:\.gz|\.zst)?$").unwrap()
    });
    static FASTX_EXT: LazyLock<Regex> =
        LazyLock::new(|| Regex::new(r"\.(?:fastq|fq|fasta|fa)(?:\.gz|\.zst)?$").unwrap());
    static ANY_EXT: LazyLock<Regex> = LazyLock::new(|| {
        Regex::new(r"\.(?:fastq|fq|fasta|fa|sam|bam|cram)(?:\.gz|\.zst)?$").unwrap()
    });

    let paired = match input_files.len() {
        1 => false,
        2 => true,
        n => bail!("Invalid number of input files: {n}"),
    };
    let input_path = input_files[0].to_str().unwrap();
    let output_name = if let Some(caps) = PAIR.captures(input_path).filter(|_| paired) {
        let suffix = caps.get(2).map_or("", |m| m.as_str());
        format!("{}{suffix}{new_extension}", &caps[1])
    } else {
        let regex = if paired { &FASTX_EXT } else { &ANY_EXT };
        regex.replace(input_path, new_extension).into_owned()
    };
    if output_name == input_path {
        bail!("Unable to autodetermine the output filename for {input_path}");
    }
    Ok(output_name)
}

pub fn pull_single_files(input_files: &[PathBuf]) -> Vec<Vec<PathBuf>> {
    static PAIR_LIKE: LazyLock<Regex> = LazyLock::new(|| Regex::new(r".+_R[12].+").unwrap());
    let num_suspect = input_files
        .iter()
        .filter(|file| PAIR_LIKE.is_match(file.to_str().unwrap()))
        .count();
    if num_suspect > 0 {
        warn!(
            "Found {num_suspect} files that may be paired but are not. If this is not intentional, consider adding the `--paired` flag."
        );
    }
    input_files.iter().map(|file| vec![file.clone()]).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_pair_r1_r2_files_sequential() {
        let files = vec![
            PathBuf::from("sample_0_R1_001.fq"),
            PathBuf::from("sample_0_R2_001.fq"),
            PathBuf::from("sample_1_R1.fastq"),
            PathBuf::from("sample_1_R2.fastq"),
        ];

        let pairs = pair_r1_r2_files(&files).unwrap();
        assert_eq!(pairs.len(), 2);

        assert_eq!(pairs[0][0], PathBuf::from("sample_0_R1_001.fq"));
        assert_eq!(pairs[0][1], PathBuf::from("sample_0_R2_001.fq"));

        assert_eq!(pairs[1][0], PathBuf::from("sample_1_R1.fastq"));
        assert_eq!(pairs[1][1], PathBuf::from("sample_1_R2.fastq"));
    }

    #[test]
    fn test_pair_r1_r2_files_non_sequential() {
        // This is the key test case you mentioned
        let files = vec![
            PathBuf::from("library_A_R1_lane1.fastq"),
            PathBuf::from("library_A_R1_lane2.fastq"),
            PathBuf::from("library_A_R2_lane1.fastq"),
            PathBuf::from("library_A_R2_lane2.fastq"),
        ];

        let pairs = pair_r1_r2_files(&files).unwrap();
        assert_eq!(pairs.len(), 2);

        // Should pair lane1 with lane1, lane2 with lane2
        assert_eq!(pairs[0][0], PathBuf::from("library_A_R1_lane1.fastq"));
        assert_eq!(pairs[0][1], PathBuf::from("library_A_R2_lane1.fastq"));

        assert_eq!(pairs[1][0], PathBuf::from("library_A_R1_lane2.fastq"));
        assert_eq!(pairs[1][1], PathBuf::from("library_A_R2_lane2.fastq"));
    }

    #[test]
    fn test_pair_r1_r2_files_missing_pairs() {
        let files = vec![
            PathBuf::from("sample_1_R1.fastq"),
            PathBuf::from("sample_2_R2.fastq"), // Missing R1
            PathBuf::from("sample_3_R1.fastq"), // Missing R2
        ];

        let pairs = pair_r1_r2_files(&files).unwrap();
        assert_eq!(pairs.len(), 0); // No complete pairs
    }

    #[test]
    fn test_large_file_list_performance() {
        // Generate a large list to ensure O(n) performance
        let mut files = Vec::new();
        for i in 0..10000 {
            files.push(PathBuf::from(format!("sample_{i:04}_R1_lane1.fastq")));
            files.push(PathBuf::from(format!("sample_{i:04}_R2_lane1.fastq")));
        }

        let pairs = pair_r1_r2_files(&files).unwrap();
        assert_eq!(pairs.len(), 10000);
    }

    #[test]
    fn test_generate_output_name_single_file() {
        let files = vec![PathBuf::from("sample_001.fastq")];
        let output = generate_output_name(&files, ".encoded").unwrap();
        assert_eq!(output, "sample_001.encoded");
    }

    #[test]
    fn test_generate_output_name_single_file_compressed() {
        let files = vec![PathBuf::from("sample_001.fastq.gz")];
        let output = generate_output_name(&files, ".encoded").unwrap();
        assert_eq!(output, "sample_001.encoded");
    }

    #[test]
    fn test_generate_output_name_paired_files() {
        let files = vec![
            PathBuf::from("sample_001_R1.fastq"),
            PathBuf::from("sample_001_R2.fastq"),
        ];
        let output = generate_output_name(&files, ".encoded").unwrap();
        assert_eq!(output, "sample_001.encoded");
    }

    #[test]
    fn test_generate_output_name_paired_files_with_suffix() {
        let files = vec![
            PathBuf::from("library_A_R1_lane1.fastq"),
            PathBuf::from("library_A_R2_lane1.fastq"),
        ];
        let output = generate_output_name(&files, ".encoded").unwrap();
        assert_eq!(output, "library_A_lane1.encoded");
    }

    #[test]
    fn test_generate_output_name_paired_files_complex_suffix() {
        let files = vec![
            PathBuf::from("sample_0_R1_001.fq.gz"),
            PathBuf::from("sample_0_R2_001.fq.gz"),
        ];
        let output = generate_output_name(&files, ".encoded").unwrap();
        assert_eq!(output, "sample_0_001.encoded");
    }

    #[test]
    fn test_generate_output_name_different_lane_numbers() {
        // This test shows that different lanes get different output names
        let files1 = vec![
            PathBuf::from("library_A_R1_lane1.fastq"),
            PathBuf::from("library_A_R2_lane1.fastq"),
        ];
        let files2 = vec![
            PathBuf::from("library_A_R1_lane2.fastq"),
            PathBuf::from("library_A_R2_lane2.fastq"),
        ];

        let output1 = generate_output_name(&files1, ".encoded").unwrap();
        let output2 = generate_output_name(&files2, ".encoded").unwrap();

        assert_eq!(output1, "library_A_lane1.encoded");
        assert_eq!(output2, "library_A_lane2.encoded");
        assert_ne!(output1, output2); // Ensure they're different
    }

    /// Expected output name, or `Err(message)` for names that cannot be determined.
    fn check_names(cases: &[(&[&str], &str, Result<&str, &str>)]) {
        for (inputs, ext, expected) in cases {
            let files: Vec<PathBuf> = inputs.iter().map(PathBuf::from).collect();
            let got = generate_output_name(&files, ext).map_err(|e| e.to_string());
            let expected = expected.map(str::to_string).map_err(str::to_string);
            assert_eq!(got, expected, "inputs={inputs:?} ext={ext}");
        }
    }

    #[test]
    fn test_generate_output_name_single_cases() {
        check_names(&[
            (&["a.fq"], ".cbq", Ok("a.cbq")),
            (&["a.fasta"], ".vbq", Ok("a.vbq")),
            (&["a.fa"], ".bq", Ok("a.bq")),
            (&["a.fa.zst"], ".bq", Ok("a.bq")),
            (&["a.fastq.gz"], ".cbq", Ok("a.cbq")),
            (&["dir.v1/a.b.fq.gz"], ".cbq", Ok("dir.v1/a.b.cbq")),
            // `_R1` is not stripped for single inputs
            (&["a_R1.fq"], ".cbq", Ok("a_R1.cbq")),
            // the extension is only replaced at the end of the name
            (&["a.fq.txt.fq"], ".cbq", Ok("a.fq.txt.cbq")),
            (&["a.sam"], ".cbq", Ok("a.cbq")),
            (&["a.bam"], ".cbq", Ok("a.cbq")),
            (&["a.cram"], ".cbq", Ok("a.cbq")),
            (&["a.bam.gz"], ".cbq", Ok("a.cbq")),
            // unchanged: nothing to replace
            (
                &["a.txt"],
                ".cbq",
                Err("Unable to autodetermine the output filename for a.txt"),
            ),
            (
                &["a.FQ"],
                ".cbq",
                Err("Unable to autodetermine the output filename for a.FQ"),
            ),
            (
                &["a.fq.bz2"],
                ".cbq",
                Err("Unable to autodetermine the output filename for a.fq.bz2"),
            ),
            (
                &["fq"],
                ".cbq",
                Err("Unable to autodetermine the output filename for fq"),
            ),
            (
                &["a.cbq"],
                ".cbq",
                Err("Unable to autodetermine the output filename for a.cbq"),
            ),
        ]);
    }

    #[test]
    fn test_generate_output_name_paired_cases() {
        check_names(&[
            (&["s_R1.fq", "s_R2.fq"], ".cbq", Ok("s.cbq")),
            (&["s_R2.fq", "s_R1.fq"], ".cbq", Ok("s.cbq")),
            (
                &["d/s_R1.fastq.gz", "d/s_R2.fastq.gz"],
                ".vbq",
                Ok("d/s.vbq"),
            ),
            (&["s_R1.fa.zst", "s_R2.fa.zst"], ".bq", Ok("s.bq")),
            (&["s_R1_001.fq", "s_R2_001.fq"], ".cbq", Ok("s_001.cbq")),
            (&["s_R1_a_b.fq", "s_R2_a_b.fq"], ".cbq", Ok("s_a_b.cbq")),
            // the last `_R[12]` that leaves a dot-free suffix is the one stripped
            (&["a_R1_R2.fq", "a_R2_R2.fq"], ".cbq", Ok("a_R1.cbq")),
            (&["a_R1_x_R2_y.fq", "b"], ".cbq", Ok("a_R1_x_y.cbq")),
            // only the first file is inspected
            (&["s_R1.fq", "other.txt"], ".cbq", Ok("s.cbq")),
            // no `_R[12]` (or a dotted suffix): fall back to swapping the extension
            (&["s_1.fq", "s_2.fq"], ".cbq", Ok("s_1.cbq")),
            (&["s.fq.gz", "t.fq.gz"], ".cbq", Ok("s.cbq")),
            (&["a_R1_x.y.fq", "a_R2_x.y.fq"], ".cbq", Ok("a_R1_x.y.cbq")),
            (
                &["/t/x_R1/a.fq", "/t/x_R2/a.fq"],
                ".cbq",
                Ok("/t/x_R1/a.cbq"),
            ),
            // the fallback does not know about sam/bam/cram
            (
                &["s_R1.bam", "s_R2.bam"],
                ".cbq",
                Err("Unable to autodetermine the output filename for s_R1.bam"),
            ),
            (
                &["s.sam", "t.sam"],
                ".cbq",
                Err("Unable to autodetermine the output filename for s.sam"),
            ),
            (
                &["s_R1.txt", "s_R2.txt"],
                ".cbq",
                Err("Unable to autodetermine the output filename for s_R1.txt"),
            ),
            (
                &["s_R1.FQ", "s_R2.FQ"],
                ".cbq",
                Err("Unable to autodetermine the output filename for s_R1.FQ"),
            ),
            // the stripped name equals the input
            (
                &["s_R1.fq", "s_R2.fq"],
                "_R1.fq",
                Err("Unable to autodetermine the output filename for s_R1.fq"),
            ),
        ]);
    }

    #[test]
    fn test_generate_output_name_invalid_count() {
        for n in [0, 3] {
            let files = vec![PathBuf::from("a.fq"); n];
            let err = generate_output_name(&files, ".cbq").unwrap_err();
            assert_eq!(
                err.to_string(),
                format!("Invalid number of input files: {n}")
            );
        }
    }
}
