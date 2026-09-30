# bqtools

A command-line utility for working with [BINSEQ](https://github.com/noamteyssier/binseq) files (`*.bq`, `*.vbq`, `*.cbq`).

See [Installation](installation.md) to get started and [Formats](formats.md) to pick a variant. Each command has a generated page under Reference.

## Commands

| Command | Description |
| --- | --- |
| [**Encode**](commands/encode.md) | Convert FASTA, FASTQ, or SAM/BAM/CRAM files to a BINSEQ format |
| [**Decode**](commands/decode.md) | Convert a BINSEQ file back to FASTA, FASTQ, or TSV format |
| [**Cat**](commands/cat.md) | Concatenate multiple BINSEQ files |
| [**Info**](commands/info.md) | Show information and statistics about one or more BINSEQ files. |
| [**Grep**](commands/grep.md) | Search for fixed-string, regex, or fuzzy matches in BINSEQ files. |
| [**Sample**](commands/sample.md) | Randomly subsample a BINSEQ file to FASTA, FASTQ, or TSV. |
| [**Split**](commands/split.md) | Split a BINSEQ file into multiple files based on matching patterns. |
| [**Pipe**](commands/pipe.md) | Create named-pipes for efficient data processing with legacy tools that don't support BINSEQ, optionally spawning and supervising the consumer commands directly (`-x`/`-X`). |
| [**Revcomp**](commands/revcomp.md) | Reverse complement the sequences in a BINSEQ file. |
| [**Verify**](commands/verify.md) | Compute an order-independent checksum over a BINSEQ file. |
| [**QC**](commands/qc.md) | Run FastQC-style quality control on a BINSEQ file. |
