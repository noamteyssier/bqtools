# bqtools

[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.md)
[![Crates.io](https://img.shields.io/crates/d/bqtools?color=orange&label=crates.io)](https://crates.io/crates/bqtools)

A command-line tool for [BINSEQ](https://github.com/noamteyssier/binseq) files.

**Documentation: <https://noamteyssier.github.io/bqtools/>**

## BINSEQ

BINSEQ is a family of binary formats for fast processing of DNA sequences:

- **CBQ** (`*.cbq`): variable-length, optional quality scores and headers, 2-bit plus `N`. Lossless by default. Use this unless you have a reason not to.
- **BQ** (`*.bq`): fixed-length, no quality scores, 2-bit or 4-bit. Fastest, lossy by design.
- **VBQ** (`*.vbq`): variable-length, optional quality scores and headers, 2-bit or 4-bit. Deprecated in favor of CBQ.

All support single and paired reads. bqtools is built on [`binseq`](https://crates.io/crates/binseq), [`bitnuc`](https://crates.io/crates/bitnuc), and [`paraseq`](https://crates.io/crates/paraseq). See the [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014181) for the format family and benchmarks.

## Installation

```bash
cargo install bqtools
```

See [Installation](https://noamteyssier.github.io/bqtools/installation/) for feature flags (`fuzzy`, `gcs`, `htslib`).

## Usage

```bash
bqtools encode sample_R1.fastq.gz sample_R2.fastq.gz -o sample.cbq
bqtools decode sample.cbq --prefix sample -f q -c g
```

See [Usage](https://noamteyssier.github.io/bqtools/usage/encoding/) for examples and [Commands](https://noamteyssier.github.io/bqtools/commands/encode/) for the full reference.

## Citation

```
Teyssier N, Dobin A (2026) BINSEQ: A family of high-performance binary formats for nucleotide sequences. PLoS Comput Biol 22(5): e1014181. https://doi.org/10.1371/journal.pcbi.1014181
```
