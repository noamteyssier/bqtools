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
        .args(["encode", "-fb", "-o"])
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

/// regression for #107: multi-threaded htslib SAM encode used to hang
#[test]
fn encode_sam_multithreaded_from_file() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("r.sam");
    let out = dir.path().join("o.cbq");
    let mut s = String::from("@HD\tVN:1.6\n");
    for i in 0..5000 {
        s += &format!("r{i}\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGT\tIIIIIIII\n");
    }
    std::fs::write(&sam, s).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_bqtools"))
        .args(["encode", "-fb", "-T", "2", "-o"])
        .arg(&out)
        .arg(&sam)
        .spawn()
        .unwrap();
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(60);
    let status = loop {
        if let Some(s) = child.try_wait().unwrap() {
            break s;
        }
        if std::time::Instant::now() > deadline {
            child.kill().unwrap();
            panic!("encode hung with -T 2");
        }
        std::thread::sleep(std::time::Duration::from_millis(100));
    };
    assert!(status.success());
    assert!(out.metadata().unwrap().len() > 0);
}
