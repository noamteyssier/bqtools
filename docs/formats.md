# Formats

BINSEQ is a binary file format family designed for high-performance processing of DNA sequences.
It currently has three variants: BQ, VBQ, and CBQ.

- **BQ (\*.bq)**: Optimized for _fixed-length_ DNA sequences **without** quality scores (2bit/4bit).
- **VBQ (\*.vbq)**: Optimized for _variable-length_ DNA sequences **with optional** quality scores, headers with 2bit/4bit.
- **CBQ (\*.cbq)**: Optimized for _variable-length_ DNA sequences **with optional** quality scores, headers with 2bit + N.

All support single and paired sequences and make use of two-bit or four-bit encoding for efficient nucleotide packing using [`bitnuc`](https://crates.io/crates/bitnuc) and efficient parallel FASTX processing using [`paraseq`](https://crates.io/crates/paraseq).

For more information about BINSEQ, see our [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014181) where we describe the format family, applications, and benchmark against other sequencing formats.

## Choosing a variant

> TL;DR: `*.cbq` is the recommended format for most applications.

For most applications the BINSEQ variant of choice is `*.cbq`.
This format is lossless by default and supports variable-length sequences.
It achieves better compression than `*.vbq` and `*.bq` by using blocked-columnar compression of sequence attributes.
It can optionally exclude quality scores and headers (but they are included by default).
For an overview of the format check out the [BINSEQ docs](https://docs.rs/binseq/latest/binseq/cbq/index.html).

If your application _only requires sequences_ and has _fixed-length_ reads then `*.bq` is the best choice.
It is the _fastest_ variant but _is lossy_ by design.

> Note: `*.vbq` was originally designed for variable-length sequences with quality scores and headers, but it is now deprecated in favor of `*.cbq` which is more compressable, lossless, and has faster decoding.
