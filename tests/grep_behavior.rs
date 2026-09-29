//! Black-box guardrails for `bqtools grep`: run the real binary against a tiny
//! hand-written fixture and pin exact stdout. Fixed-string cases in the matrix
//! run through every backend (Aho-Corasick, regex, fuzzy with `-k 0`), which
//! must agree. Only behavior that is deliberate and consistent is pinned; see
//! the notes at the bottom for known backend divergences left unpinned.

use std::{fmt::Write as _, path::PathBuf, process::Command};

use tempfile::TempDir;

// Paired fixture (R1 / R2). Headers are `r0..r3` on both mates.
const R1: [&str; 4] = ["AACCGGTTAA", "GATTACAGAT", "TTTTTTTTTT", "ACGTACGTAC"];
const R2: [&str; 4] = ["GGGGGGGGGG", "CCCCAAAAAA", "GATTACAGAT", "TTTTCCCCGG"];

const RED: &str = "\x1b[31;1m";
const RESET: &str = "\x1b[0m";

struct Fixture {
    dir: TempDir,
    paired: PathBuf,
    single: PathBuf,
}
impl Fixture {
    fn new() -> Self {
        let dir = TempDir::new().unwrap();
        let fastq = |name: &str, seqs: &[&str; 4]| {
            let mut s = String::new();
            for (i, seq) in seqs.iter().enumerate() {
                writeln!(s, "@r{i}\n{seq}\n+\n{}", "I".repeat(seq.len())).unwrap();
            }
            let path = dir.path().join(name);
            std::fs::write(&path, s).unwrap();
            path
        };
        let (r1, r2) = (fastq("a_1.fq", &R1), fastq("a_2.fq", &R2));
        let (paired, single) = (dir.path().join("p.cbq"), dir.path().join("s.cbq"));
        for (out, ins) in [(&paired, vec![&r1, &r2]), (&single, vec![&r1])] {
            let ok = Command::new(env!("CARGO_BIN_EXE_bqtools"))
                .arg("encode")
                .args(ins)
                .arg("-o")
                .arg(out)
                .status()
                .unwrap();
            assert!(ok.success());
        }
        Self {
            dir,
            paired,
            single,
        }
    }

    fn file(&self, name: &str, contents: &str) -> String {
        let path = self.dir.path().join(name);
        std::fs::write(&path, contents).unwrap();
        path.to_str().unwrap().to_string()
    }

    /// Runs `bqtools grep <input> <args>`; returns (success, stdout).
    fn run(&self, single: bool, args: &[String]) -> (bool, String) {
        let input = if single { &self.single } else { &self.paired };
        let out = Command::new(env!("CARGO_BIN_EXE_bqtools"))
            .arg("grep")
            .arg(input)
            .args(args)
            .output()
            .unwrap();
        (out.status.success(), String::from_utf8(out.stdout).unwrap())
    }

    fn ok(&self, single: bool, args: &[&str]) -> String {
        let args: Vec<String> = args.iter().map(ToString::to_string).collect();
        let (ok, stdout) = self.run(single, &args);
        assert!(ok, "grep {args:?} failed");
        stdout
    }
}

/// Expected TSV lines (sorted) for both mates of the given records.
fn pair(recs: &[usize]) -> Vec<String> {
    let mut v: Vec<_> = recs
        .iter()
        .flat_map(|&i| [format!("r{i}\t{}", R1[i]), format!("r{i}\t{}", R2[i])])
        .collect();
    v.sort();
    v
}
fn mate(seqs: &[&str; 4], recs: &[usize]) -> Vec<String> {
    recs.iter().map(|&i| format!("r{i}\t{}", seqs[i])).collect()
}
fn sorted_lines(s: &str) -> Vec<String> {
    let mut v: Vec<_> = s.lines().map(str::to_string).collect();
    v.sort();
    v
}

/// Expands `@PATTERN` args for a backend: 0 = auto (Aho-Corasick), 1 = regex
/// (last base wrapped in a class so it is no longer a fixed string), 2 = fuzzy
/// with `-k 0` (exact).
fn expand(args: &[&str], backend: usize) -> Vec<String> {
    let mut out: Vec<String> = args
        .iter()
        .map(|a| match a.strip_prefix('@') {
            Some(p) if backend == 1 => format!("{}[{}]", &p[..p.len() - 1], &p[p.len() - 1..]),
            Some(p) => p.to_string(),
            None => a.to_string(),
        })
        .collect();
    if backend == 2 {
        out.extend(["-z", "-k", "0"].map(String::from));
    }
    out
}

fn backends() -> &'static [usize] {
    if cfg!(feature = "fuzzy") {
        &[0, 1, 2]
    } else {
        &[0, 1]
    }
}

/// Runs a matrix case on every backend; `norm` maps stdout to comparable lines.
fn matrix(fx: &Fixture, single: bool, args: &[&str], want: &[String], sort: bool) {
    for &b in backends() {
        let (ok, stdout) = fx.run(single, &expand(args, b));
        assert!(ok, "backend {b}: grep {args:?} failed");
        // `-P` prints pattern text as the row name; undo the regex wrapping
        let stdout = if b == 1 && !stdout.contains('\x1b') {
            stdout.replace(['[', ']'], "")
        } else {
            stdout
        };
        let got: Vec<String> = if sort {
            sorted_lines(&stdout)
        } else {
            stdout.lines().map(str::to_string).collect()
        };
        assert_eq!(got, want, "backend {b}: grep {args:?}");
    }
}

#[test]
fn record_selection() {
    let fx = Fixture::new();
    let (p, s) = (false, true);
    let m = |single, args: &[&str], want: Vec<String>| matrix(&fx, single, args, &want, true);

    // pattern sets: either / primary (-r) / extended (-R)
    m(p, &["@GATTA"], pair(&[1, 2]));
    m(p, &["-r", "@GATTA"], pair(&[1]));
    m(p, &["-R", "@GATTA"], pair(&[2]));
    m(p, &["@CCCCCCCC"], vec![]);

    // AND (default) vs OR across sets and within a set
    m(p, &["-r", "@GATTA", "-R", "@CCCC"], pair(&[1]));
    m(
        p,
        &["-r", "@GATTA", "-R", "@CCCC", "--or-logic"],
        pair(&[1, 3]),
    );
    m(p, &["@GATT", "@CCCC"], pair(&[1]));
    m(p, &["@GATT", "@CCCC", "--or-logic"], pair(&[1, 2, 3]));

    // invert, including inverted AND / OR
    m(p, &["-v", "@GATTA"], pair(&[0, 3]));
    m(p, &["-v", "@GATT", "@CCCC"], pair(&[0, 2, 3]));
    m(p, &["-v", "@GATT", "@CCCC", "--or-logic"], pair(&[0]));

    // --mate restricts patterns AND output to that mate
    m(p, &["-m", "1", "@TTTT"], mate(&R1, &[2]));
    m(p, &["-m", "2", "@TTTT"], mate(&R2, &[3]));
    m(p, &["-m", "both", "@TTTT"], pair(&[2, 3]));

    // --range is 0-based, end-exclusive, clamped
    m(p, &["--range", "0..5", "@GATTA"], pair(&[1, 2]));
    m(p, &["--range", "1..", "@GATTA"], vec![]);
    m(p, &["--range", "..4", "@GATTA"], vec![]);
    m(p, &["--range", "2..7", "@TTACA"], pair(&[1, 2]));
    m(p, &["--range", "3..", "@TTACA"], vec![]);
    m(p, &["--range", "..100", "@GATTA"], pair(&[1, 2]));

    // single-end: --mate is ignored
    m(s, &["@GATTA"], mate(&R1, &[1]));
    m(s, &["-m", "2", "@GATTA"], mate(&R1, &[1]));
    m(s, &["-v", "@GATTA"], mate(&R1, &[0, 2, 3]));
    m(s, &["@AC", "@GT"], mate(&R1, &[0, 3]));
    m(s, &["@AC", "@GT", "--or-logic"], mate(&R1, &[0, 1, 3]));

    // knobs that must not change results
    m(p, &["--no-dfa", "@GATTA"], pair(&[1, 2]));
    m(p, &["-T", "4", "@GATTA"], pair(&[1, 2]));
    m(p, &["--span", "2..4", "@GATTA"], pair(&[2]));
}

#[test]
fn count_modes() {
    let fx = Fixture::new();
    let (p, s) = (false, true);
    let m = |single, args: &[&str], want: &str| {
        let want: Vec<String> = want.lines().map(str::to_string).collect();
        matrix(&fx, single, args, &want, false);
    };

    m(p, &["-C", "@GATTA"], "2");
    m(p, &["-C", "-v", "@GATTA"], "2");
    m(p, &["-C", "@CCCCCCCC"], "0");
    m(p, &["-C", "--span", "0..1", "@CCCCCCCC"], "0");
    m(p, &["-C", "--span", "2..4", "@GATTA"], "1");
    m(s, &["-C", "@GATTA"], "1");
    m(p, &["-F", "@GATTA"], "count\ttotal\tfrac\n2\t4\t0.5000");
    m(p, &["-F", "@CCCCCCCC"], "count\ttotal\tfrac\n0\t4\t0.0000");
    m(
        p,
        &["-F", "-v", "@GATT", "@CCCC"],
        "count\ttotal\tfrac\n3\t4\t0.7500",
    );

    // -P: one row per pattern, ordered primary (-r), extended (-R), either
    let head = "name\tcount\tfrac_total\n";
    m(
        p,
        &["-P", "@AAAA", "@CCCC"],
        &format!("{head}AAAA\t1\t0.25\nCCCC\t2\t0.5"),
    );
    m(
        p,
        &["-P", "-v", "@AAAA", "@CCCC"],
        &format!("{head}AAAA\t3\t0.75\nCCCC\t2\t0.5"),
    );
    m(
        p,
        &["-P", "-r", "@AAAA", "-R", "@AAAA"],
        &format!("{head}AAAA\t0\t0.0\nAAAA\t1\t0.25"),
    );
    m(
        p,
        &["-P", "-r", "@GATTA", "-R", "@CCCC", "@TTTT"],
        &format!("{head}GATTA\t1\t0.25\nCCCC\t2\t0.5\nTTTT\t2\t0.5"),
    );
    m(
        p,
        &["-P", "--range", "0..5", "@GATTA"],
        &format!("{head}GATTA\t2\t0.5"),
    );
    m(
        s,
        &["-P", "@GATT", "@TTTT"],
        &format!("{head}GATT\t1\t0.25\nTTTT\t1\t0.25"),
    );
}

#[test]
fn color() {
    let fx = Fixture::new();
    let (p, s) = (false, true);
    let hit = |pre: &str, m: &str, post: &str| format!("{pre}{RED}{m}{RESET}{post}");
    let m = |single, args: &[&str], want: Vec<String>| {
        let mut args = args.to_vec();
        args.extend(["--color", "always"]);
        matrix(&fx, single, &args, &want, true);
    };

    m(
        s,
        &["@GATTA"],
        vec![format!("r1\t{}", hit("", "GATTA", "CAGAT"))],
    );
    // highlight positions stay absolute under --range
    m(
        s,
        &["--range", "2..", "@TTACA"],
        vec![format!("r1\t{}", hit("GA", "TTACA", "GAT"))],
    );
    // only the mate containing the match is highlighted
    let mut want = vec![
        format!("r1\t{}", hit("", "GATTA", "CAGAT")),
        format!("r1\t{}", R2[1]),
        format!("r2\t{}", R1[2]),
        format!("r2\t{}", hit("", "GATTA", "CAGAT")),
    ];
    want.sort();
    m(p, &["@GATTA"], want);
    // inverted matches highlight nothing
    m(s, &["-v", "@GATTA"], mate(&R1, &[0, 2, 3]));

    // overlapping hits merge into one span (fixed strings only: the regex
    // backend stops at the first matching pattern under OR, see below)
    assert_eq!(
        fx.ok(true, &["GATTA", "TTACA", "--or-logic", "--color", "always"]),
        format!("r1\t{}\n", hit("", "GATTACA", "GAT"))
    );

    // fastq colors sequence and quality with the same spans
    assert_eq!(
        fx.ok(true, &["GATTA", "-f", "q", "--color", "always"]),
        format!(
            "@r1\n{}\n+\n{}\n",
            hit("", "GATTA", "CAGAT"),
            hit("", "IIIII", "IIIII")
        )
    );

    // no color unless asked for when stdout is not a terminal
    assert_eq!(fx.ok(true, &["GATTA"]), "r1\tGATTACAGAT\n");
    assert_eq!(
        fx.ok(true, &["GATTA", "--color", "never"]),
        "r1\tGATTACAGAT\n"
    );
}

#[test]
fn regex_semantics() {
    let fx = Fixture::new();
    let g = |args: &[&str]| sorted_lines(&fx.ok(false, args));

    assert_eq!(g(&["A{3}"]), pair(&[1]));
    assert_eq!(g(&["^GATTA"]), pair(&[1, 2]));
    assert_eq!(g(&["GAT$"]), pair(&[1, 2]));
    // anchors apply to the sliced range, not the whole sequence
    assert_eq!(g(&["^ATTA", "--range", "1.."]), pair(&[1, 2]));
    assert_eq!(g(&["^ATTA"]), Vec::<String>::new());
    // -x forces fixed-string matching even for non-ACGT text
    assert_eq!(g(&["-x", "GATTA"]), pair(&[1, 2]));
}

#[test]
fn headers() {
    let fx = Fixture::new();
    let g = |args: &[&str]| sorted_lines(&fx.ok(false, args));

    assert_eq!(g(&["-H", "r1"]), pair(&[1]));
    assert_eq!(g(&["-H", "-v", "r1"]), pair(&[0, 2, 3]));
    assert_eq!(g(&["-H", "-r", "r[23]"]), pair(&[2, 3]));
    assert_eq!(fx.ok(false, &["-H", "-C", "r[12]"]), "2\n");
    assert_eq!(
        fx.ok(false, &["-H", "-P", "r1", "r2"]),
        "name\tcount\tfrac_total\nr1\t1\t0.25\nr2\t1\t0.25\n"
    );
    // header text is never colored
    assert_eq!(g(&["-H", "r1", "--color", "always"]), pair(&[1]));
}

#[test]
fn reverse_complement() {
    let fx = Fixture::new();
    // GATTACA <-> TGTAATC
    assert_eq!(fx.ok(true, &["--rc", "TGTAATC"]), "r1\tGATTACAGAT\n");
    assert_eq!(
        sorted_lines(&fx.ok(false, &["--rc", "TGTAATC"])),
        pair(&[1, 2])
    );
    assert_eq!(fx.ok(true, &["TGTAATC"]), "");
    assert_eq!(fx.ok(true, &["--rc", "-C", "AAAA"]), "1\n"); // TTTT in r2
}

#[test]
fn pattern_files() {
    let fx = Fixture::new();
    let fa = fx.file("pats.fa", ">alpha\nGATTA\n>beta\nAAAA\n");
    let tsv = fx.file("pats.tsv", "x\tGATTA\ny\tCCCC\n");
    let txt = fx.file("pats.txt", "GATTA\nCCCC\n");
    let head = "name\tcount\tfrac_total\n";

    // pattern files always combine with OR (AND would give only r1)
    assert_eq!(sorted_lines(&fx.ok(false, &["--file", &fa])), pair(&[1, 2]));
    assert_eq!(fx.ok(false, &["-C", "--file", &txt]), "3\n");
    // -P reports FASTA/TSV names, or the sequence when unnamed
    assert_eq!(
        fx.ok(false, &["-P", "--file", &fa]),
        format!("{head}alpha\t2\t0.5\nbeta\t1\t0.25\n")
    );
    assert_eq!(
        fx.ok(false, &["-P", "--file", &tsv]),
        format!("{head}x\t2\t0.5\ny\t2\t0.5\n")
    );
    assert_eq!(
        fx.ok(false, &["-P", "--file", &txt]),
        format!("{head}GATTA\t2\t0.5\nCCCC\t2\t0.5\n")
    );
    // --sfile / --xfile restrict to one mate
    assert_eq!(
        fx.ok(false, &["-P", "--sfile", &fa]),
        format!("{head}alpha\t1\t0.25\nbeta\t0\t0.0\n")
    );
    assert_eq!(
        fx.ok(false, &["-P", "--xfile", &fa]),
        format!("{head}alpha\t1\t0.25\nbeta\t1\t0.25\n")
    );
    // CLI and file patterns are merged into one set
    assert_eq!(
        fx.ok(false, &["-P", "--file", &fa, "TTTT"]),
        format!("{head}TTTT\t2\t0.5\nalpha\t2\t0.5\nbeta\t1\t0.25\n")
    );
}

#[test]
fn rejected_invocations() {
    let fx = Fixture::new();
    let fails = |single, args: &[&str]| {
        let args: Vec<String> = args.iter().map(ToString::to_string).collect();
        assert!(!fx.run(single, &args).0, "grep {args:?} should fail");
    };
    fails(false, &[]); // no patterns
    fails(false, &["-m", "1", "-R", "GATTA"]); // nothing left for mate 1
    fails(false, &["-m", "2", "-r", "GATTA"]); // nothing left for mate 2
    fails(true, &["-R", "GATTA"]); // no extended mate to search
    fails(true, &["-r", "GATTA", "-R", "CCCC"]);
    fails(true, &["-m", "1", "-R", "GATTA"]);
    let xfile = fx.file("x.txt", "GATTA\n");
    fails(true, &["--xfile", &xfile]);
    // primary and either patterns stay valid on single-end input
    assert_eq!(fx.ok(true, &["--sfile", &xfile]), "r1\tGATTACAGAT\n");
    fails(false, &["A("]); // invalid regex
    fails(false, &["--rc", "AC.GT"]); // rc needs fixed ACGT
    fails(false, &["-C", "-P", "GATTA"]); // exclusive modes
    #[cfg(feature = "fuzzy")]
    fails(false, &["-z", "GATTA", "CC"]); // mixed pattern lengths
}

// Known backend divergences, deliberately NOT pinned (revisit when unifying):
// * regex under OR stops at the first matching pattern, so only that
//   pattern's hits are highlighted; Aho-Corasick and fuzzy highlight all.
