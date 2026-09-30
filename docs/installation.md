# Installation

## From Cargo

bqtools can be installed using `cargo`, the Rust package manager:

```bash
cargo install bqtools
```

To install `cargo` you can follow the instructions on the [official Rust website](https://www.rust-lang.org/tools/install).

## From Source

```bash
# Clone the repository
git clone https://github.com/noamteyssier/bqtools.git
cd bqtools

# Install
cargo install --path .

# Check installation
bqtools --help
```

## Feature Flags

bqtools supports the following feature flags:

- `htslib`: Enable support for reading SAM/BAM/CRAM files using the [`htslib`](https://docs.rs/rust-htslib/latest/rust_htslib/) library (default).
- `gcs`: Enable support for reading Google Cloud Storage files.
- `fuzzy`: Enable fuzzy matching in the `grep` and `split` commands using the [`sassy`](https://crates.io/crates/sassy) library

To enable fuzzy matching, `bqtools` must be compiled using a `native` target cpu:

```bash
# Install from source
export RUSTFLAGS="-C target-cpu=native"; cargo install --path . -F fuzzy;

# Or install from crates but enforce native target cpu
export RUSTFLAGS="-C target-cpu=native"; cargo install bqtools -F fuzzy;
```

To selectively enable/disable feature flags:

```bash
# (for fuzzy matching support sassy requires native target cpu)
export RUSTFLAGS="-C target-cpu=native";

# Install bqtools without htslib but with fuzzy matching
cargo install bqtools --no-default-features -F fuzzy
#
# Install bqtools without htslib but with fuzzy matching and gcs
cargo install bqtools --no-default-features -F fuzzy,gcs
```
