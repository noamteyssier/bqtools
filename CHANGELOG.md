# Changelog

## 0.6.0

A dependency-refresh and internal cleanup release, with a grep engine consolidation, a few bug fixes, and one new flag (`pipe --span`). Also, `encode` no longer logs an `error!` line before bailing, so those errors now print once, via `anyhow`. The work landed as stacked PRs (#7–#80) on `dev-0.6.0`.

### Features

- `pipe` supports `--span` to restrict processing to a range of records (#80)

### Fixes

- `grep` rejects extended patterns on single-end input instead of accepting them
- `grep` fuzzy matching skips empty sequences
- `pipe` errors when the output path already exists and is not a fifo. Symlinks are not followed (#76)

### Dependencies

- upgrade dependencies, with no breaking changes (#7)
  - `binseq` 0.9.6 → 0.10.0
  - `xxhash-rust` 0.8.18 → 0.8.19
- drop `is-terminal` in favor of `std::io::IsTerminal` (#8)
- drop `num_cpus` in favor of `std::thread::available_parallelism`, which falls back to 1 (#9)
- drop `parking_lot` in favor of `std::sync::Mutex` (#10)
- move `niffler` to `dev-dependencies`, since only test helpers use it (#11)
- drop `nix` in favor of `libc` for fifo creation (#76)
- drop `memmap2`, and with it an `unsafe` mmap. `cat` now copies the file tail with `seek` + `io::copy` (#28)

### Internal refactors

Behavior is unchanged. The changes remove duplicated code and simplify each subcommand.

- **shared plumbing**
  - inline the gzip and zstd passthrough helpers into `compress_passthrough` (#13)
  - simplify `match_output` (#14) and `create_fifos` (#16)
  - dedupe the `BoxedWriter` alias (#15)
  - use `fs::create_dir_all` instead of `make_directory` (#17)
  - inline single-use fuzzy helpers (#18)
  - simplify the output option accessors (#19) and input helpers (#20)
  - add a shared writer builder based on the input header (#34)
  - add a shared placeholder quality fill (#22)
  - add `SplitWriter::write_batch` (#21)
- **decode and sample**
  - simplify decode writers and options (#23)
  - build `sample` on the decode processor (#25)
- **cat**
  - share the header check across formats (#27)
- **info**
  - `--num` reads the record count directly instead of building the full info (#29)
  - add a `section()` banner helper (#29)
- **verify**
  - accumulate the checksum atomically (#31)
  - simplify record hashing (#32) and field and mate labels (#33)
- **revcomp**
  - simplify record building (#36)
- **pipe**
  - build record pairs once (#38)
  - simplify consumer spawning (#39) and the processor (#40)
  - trim helpers (#41)
- **encode**
  - dedupe `encode_collection` calls (#43)
  - simplify the encoder processor (#44)
  - drop the debug interval counter (#45)
  - share the first-records length probe (#46)
  - share the writer builder and finish step (#47)
  - simplify file queueing (#48) and output naming (#49)
  - drop duplicate error logs before bailing (#50)
- **split**
  - simplify setup (#52)
  - share alias binning across splitters (#53)
  - simplify the processor (#54)
- **grep** (#63–#77)
  - share one pattern engine between `-P` and `split`, and run the filter on it
  - share pattern sets and the aho-corasick builder between `grep` and `split`
  - write colored output through the shared record writer
  - store match ranges in a vec, use atomics for shared counts, and simplify range parsing and slicing
  - drop the pattern collection drain and unused `pattern_strings`
- **cli**
  - resolve `--span` through `InputBinseq::range`
- **qc**
  - build modules directly from `QcOptions` and delete `qc/config.rs` (#56)
  - share TSV writing and stats (#57)
  - share read-pair plumbing (#58)
  - collapse module dispatch (#59)
  - trim leftovers (#60)

### Tests and docs

- share test helpers for pipe (#42), split (#55), sample, and encode (`testutils`, #26)
- reuse the decode helper across tests (#24)
- merge revcomp tests (#37)
- fold encode specialization tests into one (#51)
- trim the info record count test (#30)
- add tests covering encode output naming
- add grep behavior guardrail tests and a `bench_grep` timing example
- use lowercase "htslib" in docs
