# Formats

BINSEQ is a family of binary formats for fast processing of DNA sequences. There are currently three variants:

- **CBQ** (`*.cbq`): variable-length sequences with optional quality scores and headers. 2-bit encoding plus `N`.

- **BQ** (`*.bq`): fixed-length sequences, no quality scores. 2-bit or 4-bit encoding.

- **VBQ** (`*.vbq`): variable-length sequences with optional quality scores and headers. 2-bit or 4-bit encoding.

All three support single and paired reads. Nucleotides are packed with [`bitnuc`](https://crates.io/crates/bitnuc) and FASTX input is processed in parallel with [`paraseq`](https://crates.io/crates/paraseq).

For the format family, its applications, and benchmarks against other formats, see the [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014181).

## Choosing a variant

!!! tip
    Use `*.cbq` unless you have a reason not to.

`*.cbq` is lossless by default and handles variable-length reads. Blocked-columnar compression of sequence attributes gives it better compression than `*.vbq` and `*.bq`. Quality scores and headers are kept by default and can be excluded. See the [BINSEQ docs](https://docs.rs/binseq/latest/binseq/cbq/index.html) for details.

Use `*.bq` if you only need sequences and your reads are fixed-length. It is the fastest variant, but lossy by design.

!!! note
    `*.vbq` is deprecated. `*.cbq` compresses better, is lossless, and decodes faster.
