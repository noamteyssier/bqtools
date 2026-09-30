# bqtools

A command-line tool for [BINSEQ](https://github.com/noamteyssier/binseq) files (`*.bq`, `*.vbq`, `*.cbq`).

Start with [Installation](installation.md), then see [Formats](formats.md) to choose a variant and [Usage](usage/encoding.md) for worked examples. Every command has its own reference page.

## Commands

| Command | Description |
| --- | --- |
| [`encode`](commands/encode.md) | Convert FASTA, FASTQ, or SAM/BAM/CRAM to BINSEQ |
| [`decode`](commands/decode.md) | Convert BINSEQ to FASTA, FASTQ, or TSV |
| [`cat`](commands/cat.md) | Concatenate BINSEQ files |
| [`info`](commands/info.md) | Show statistics for BINSEQ files |
| [`grep`](commands/grep.md) | Search for fixed-string, regex, or fuzzy matches |
| [`sample`](commands/sample.md) | Randomly subsample to FASTA, FASTQ, or TSV |
| [`split`](commands/split.md) | Split records into files by matching pattern |
| [`pipe`](commands/pipe.md) | Stream records through named pipes to tools that don't read BINSEQ, optionally running and supervising them (`-x`/`-X`) |
| [`revcomp`](commands/revcomp.md) | Reverse complement sequences |
| [`verify`](commands/verify.md) | Compute an order-independent checksum |
| [`qc`](commands/qc.md) | Run FastQC-style quality control |
