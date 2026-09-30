#![cfg(feature = "htslib")]
use std::{
    io::Write,
    process::{Command, Stdio},
};

#[test]
fn encode_sam_from_stdin() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("o.cbq");
    let mut child = Command::new(env!("CARGO_BIN_EXE_bqtools"))
        // -T 1: htslib encode hangs with more than one thread (also from a file path)
        .args(["encode", "-fb", "-T", "1", "-o"])
        .arg(&out)
        .stdin(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(b"@HD\tVN:1.6\nr0\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n")
        .unwrap();
    assert!(child.wait().unwrap().success());
    assert!(out.metadata().unwrap().len() > 0);
}
