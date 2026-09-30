# Encoding

Convert FASTQ, FASTA, or SAM/BAM/CRAM to BINSEQ. Input format and compression (gzip/zstd) are auto-detected.

## Single and paired files

```bash
# single-end
bqtools encode reads.fastq -o reads.cbq

# paired-end: pass R1 then R2
bqtools encode sample_R1.fastq.gz sample_R2.fastq.gz -o sample.cbq

# without -o, the name is derived from the input (here: sample.cbq)
bqtools encode sample_R1.fastq.gz sample_R2.fastq.gz
```

The variant comes from [`-m`](../commands/encode.md#encode--mode), else the `-o` extension, else `cbq`.

```bash
# extension picks the variant
bqtools encode reads.fastq -o reads.bq

# -m wins over the extension
bqtools encode reads.fastq -o reads.out -m cbq
```

## Streams

With no input path, `encode` reads stdin. `-o` is required.

```bash
# decompress with any tool and stream in
zstdcat reads.fastq.zst | bqtools encode -o reads.cbq

# interleaved FASTQ (alternating R1/R2) becomes a paired file
bqtools encode interleaved.fastq -I -o sample.cbq
```

## SAM/BAM/CRAM

```bash
# detected from the extension
bqtools encode aligned.bam -o reads.cbq

# non-standard extension: set the format with -fb
bqtools encode aligned.bam.tmp -fb -o reads.cbq

# paired alignments must be name-sorted; -I pairs consecutive records
samtools sort -n -o sorted.bam aligned.bam
bqtools encode sorted.bam -I -o sample.cbq
```

!!! note
    SAM/BAM/CRAM needs the `htslib` feature (on by default). Stdin works for `.cbq` and `.vbq`; `.bq` needs a file path.

## Headers, qualities, and genomes

```bash
# drop headers
bqtools encode reads.fastq -o reads.cbq -H

# drop quality scores
bqtools encode reads.fastq -o reads.cbq -Q

# genomes and long references: 4-bit vbq that keeps Ns and headers
bqtools encode genome.fa -o genome.vbq -A
```

!!! note
    [`-A`](../commands/encode.md#encode--archive) does not change the variant. Pair it with a `.vbq` output or `-m vbq`.

!!! tip
    Genomes have few, long records. A small [`-b`](../commands/encode.md#encode--batch-size) (e.g. `-b 2`) spreads them across threads.

## Many files

More than two inputs are each encoded to their own output, in parallel.

```bash
# one output per file
bqtools encode *.fastq.gz

# one output per R1/R2 pair
bqtools encode *.fastq.gz --paired

# everything into one file
bqtools encode *.fastq.gz --paired --collate -o run.cbq

# paths listed one per line
bqtools encode --manifest files.txt --paired
```

!!! warning
    Exactly two inputs are always treated as one R1/R2 pair. A glob that matches two single-end files encodes one paired file.

## Directories

```bash
# every FASTX file under ./run, written next to its input
bqtools encode --recursive ./run

# pairs only, at most two levels deep
bqtools encode --recursive --paired --depth 2 ./run
```

See the [encode reference](../commands/encode.md) for every flag.
