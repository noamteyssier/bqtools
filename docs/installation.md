# Installation

## From crates.io

```bash
cargo install bqtools
```

Install `cargo` with [rustup](https://www.rust-lang.org/tools/install).

## From source

```bash
git clone https://github.com/noamteyssier/bqtools.git
cd bqtools
cargo install --path .

bqtools --help  # verify
```

## Feature flags

| Flag | Enables | Default |
| --- | --- | --- |
| `htslib` | Reading SAM/BAM/CRAM via [`htslib`](https://docs.rs/rust-htslib/latest/rust_htslib/) | yes |
| `gcs` | Reading from Google Cloud Storage | no |
| `fuzzy` | Fuzzy matching in `grep` and `split` via [`sassy`](https://crates.io/crates/sassy) | no |

`fuzzy` requires building for the native CPU:

```bash
export RUSTFLAGS="-C target-cpu=native"

# from crates.io
cargo install bqtools -F fuzzy

# from source
cargo install --path . -F fuzzy
```

To change the defaults, disable them and list the flags you want:

```bash
# no htslib, with fuzzy matching
cargo install bqtools --no-default-features -F fuzzy

# no htslib, with fuzzy matching and gcs
cargo install bqtools --no-default-features -F fuzzy,gcs
```
