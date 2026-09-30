# Decoding

Get FASTQ, FASTA, or TSV back out, subsample, merge, and reverse complement.

## Decode

Format and compression are inferred from the `-o` extension. On stdout the default is TSV.

```bash
# FASTQ, gzip-compressed
bqtools decode reads.cbq -o reads.fastq.gz

# FASTA
bqtools decode reads.cbq -o reads.fa

# FASTQ to stdout
bqtools decode reads.cbq -f q | head
```

## Paired files

Paired records are interleaved unless split with [`--prefix`](../commands/decode.md#decode--prefix).

```bash
# writes sample_R1.fq.gz and sample_R2.fq.gz
bqtools decode sample.cbq --prefix sample -f q -c g

# only R1
bqtools decode sample.cbq -o sample_R1.fastq.gz -m 1
```

## Peek

```bash
# first 10 records as FASTQ
bqtools decode reads.cbq --span ..10 -f q

# records 1000..2000 (0-based, end-exclusive)
bqtools decode reads.cbq --span 1000..2000 -o slice.fq
```

## Subsample

Each record is kept with probability [`-F`](../commands/sample.md#sample--fraction), so output size is approximate. Output options match `decode`.

```bash
# ~10% of reads
bqtools sample reads.cbq -F 0.1 -o subset.fastq.gz

# same seed, same records, at any thread count
bqtools sample reads.cbq -F 0.1 -S 7 -o subset.fq

# paired into subset_R1.fq / subset_R2.fq
bqtools sample sample.cbq -F 0.5 --prefix subset -f q
```

## Concatenate

```bash
bqtools cat lane1.cbq lane2.cbq -o merged.cbq
```

!!! note
    All inputs must be the same variant with identical headers (paired, quality, headers, bitsize). For vbq/cbq, record order is not preserved.

## Reverse complement

The output keeps the input's variant, so the extension must match.

```bash
# both mates
bqtools revcomp sample.cbq -o sample.rc.cbq

# only R2
bqtools revcomp sample.cbq -o sample.rc.cbq -M 2
```

!!! tip
    `cat` and `revcomp` never write binary to stdout on their own. Use `--pipe` to stream.
