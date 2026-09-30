# bqtools

A command-line utility for working with [BINSEQ](https://github.com/noamteyssier/binseq) files (`*.bq`, `*.vbq`, `*.cbq`).

See [Installation](installation.md) to get started and [Formats](formats.md) to pick a variant. Each command has a generated page under Reference.

## Commands

- **Encode**: Convert FASTA, FASTQ, or SAM/BAM/CRAM files to a BINSEQ format
- **Decode**: Convert a BINSEQ file back to FASTA, FASTQ, or TSV format
- **Cat**: Concatenate multiple BINSEQ files
- **Info**: Show information and statistics about one or more BINSEQ files.
- **Grep**: Search for fixed-string, regex, or fuzzy matches in BINSEQ files.
- **Sample**: Randomly subsample a BINSEQ file to FASTA, FASTQ, or TSV.
- **Split**: Split a BINSEQ file into multiple files based on matching patterns.
- **Pipe**: Create named-pipes for efficient data processing with legacy tools that don't support BINSEQ, optionally spawning and supervising the consumer commands directly (`-x`/`-X`).
- **Revcomp**: Reverse complement the sequences in a BINSEQ file.
- **Verify**: Compute an order-independent checksum over a BINSEQ file.
- **QC**: Run FastQC-style quality control on a BINSEQ file.
