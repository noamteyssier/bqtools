# bqtools

A command-line tool for [BINSEQ](https://github.com/noamteyssier/binseq) files (`*.bq`, `*.vbq`, `*.cbq`).

Start with [Installation](installation.md), then see [Formats](formats.md) to choose a variant. Every command has its own reference page.

## Commands

| Command | Description |
| --- | --- |
| [**Encode**](commands/encode.md) | Convert FASTA, FASTQ, or SAM/BAM/CRAM to BINSEQ |
| [**Decode**](commands/decode.md) | Convert BINSEQ to FASTA, FASTQ, or TSV |
| [**Cat**](commands/cat.md) | Concatenate BINSEQ files |
| [**Info**](commands/info.md) | Show statistics for BINSEQ files |
| [**Grep**](commands/grep.md) | Search for fixed-string, regex, or fuzzy matches |
| [**Sample**](commands/sample.md) | Randomly subsample to FASTA, FASTQ, or TSV |
| [**Split**](commands/split.md) | Split records into files by matching pattern |
| [**Pipe**](commands/pipe.md) | Stream records through named pipes to tools that don't read BINSEQ, optionally running and supervising them (`-x`/`-X`) |
| [**Revcomp**](commands/revcomp.md) | Reverse complement sequences |
| [**Verify**](commands/verify.md) | Compute an order-independent checksum |
| [**QC**](commands/qc.md) | Run FastQC-style quality control |
