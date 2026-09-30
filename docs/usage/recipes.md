# Recipes

Short workflows that combine commands.

## Convert and confirm a round trip

```bash
bqtools encode sample_R1.fastq.gz sample_R2.fastq.gz -o sample.cbq
bqtools decode sample.cbq --prefix roundtrip -f q -c g
bqtools encode roundtrip_R1.fq.gz roundtrip_R2.fq.gz -o roundtrip.cbq

# same checksum = same records, in any order
bqtools verify sample.cbq
bqtools verify roundtrip.cbq
```

## Count reads across a run

```bash
bqtools encode ./run/*.fastq.gz --paired
bqtools info ./run/*.cbq --num | awk '{s += $1} END {print s}'
```

## Subsample, then search

```bash
bqtools sample sample.cbq -F 0.01 -S 42 -o sub.fq
bqtools encode sub.fq -o sub.cbq
bqtools grep sub.cbq --sfile patterns.tsv -P
```

!!! tip
    For a quick estimate, [`--span`](../commands/grep.md#grep--span) is cheaper than sampling: `bqtools grep sample.cbq --sfile patterns.tsv -P --span ..100000`. It takes the first records, which may not be representative.

## Demultiplex and tally

```bash
bqtools split sample.cbq --sfile patterns.tsv --basepath demux
bqtools info demux/*.cbq --num | sort -k1,1nr
```

## Feed a tool without decoding to disk

```bash
# instead of decode -> fastq -> tool
bqtools pipe sample.cbq -p 8 -x 'bowtie2 -x ref -1 {R1} -2 {R2} -S shard_{n}.sam'
```

## Encode from Google Cloud Storage

Needs the `gcs` feature. `gs://` paths work anywhere `encode` takes an input.

```bash
bqtools encode gs://bucket/sample_R1.fastq.gz gs://bucket/sample_R2.fastq.gz -o sample.cbq
```
