# Piping to other tools

Most tools don't read BINSEQ yet. `pipe` streams records as FASTQ (or FASTA) through named pipes (FIFOs), so tools read them like files with no intermediate files on disk.

Each FIFO carries a slice of the file, so tools that take several inputs or run one process per file get parallelism for free.

!!! note
    `pipe` is not available on Windows.

## Run a command per pipe

[`-x`](../commands/pipe.md#pipe--exec) runs one shell command per FIFO and waits for all of them. Tokens:

| Token | Expands to |
| --- | --- |
| `{}` | FIFO path (single-end) |
| `{R1}`, `{R2}` | R1/R2 FIFO paths (paired) |
| `{n}` | pipe index, for per-shard outputs |

```bash
# count reads in 4 shards
bqtools pipe reads.cbq -p 4 -x 'wc -l {} > shard_{n}.lines'

# compress 4 shards in parallel
bqtools pipe reads.cbq -p 4 -x 'gzip -c < {} > shard_{n}.fq.gz'

# paired: -p 8 gives 4 R1/R2 pairs
bqtools pipe sample.cbq -p 8 -x 'bowtie2 -x ref -1 {R1} -2 {R2} -S shard_{n}.sam'
```

!!! tip
    Reference only `{R1}` (or `{R2}`) to process one mate. The other mate's FIFOs are never created.

    ```bash
    bqtools pipe sample.cbq -p 4 -x 'gzip -c < {R1} > r1_{n}.fq.gz'
    ```

## Run one command on all pipes

[`-X`](../commands/pipe.md#pipe--exec-batch) runs a single command with every FIFO path substituted. Use it for tools that take many inputs and parallelize internally.

```bash
# single-end: {} is every FIFO path
bqtools pipe reads.cbq -p 4 -X 'cat {} | wc -l'

# paired: adjacent {R1} {R2} interleaves as r1_0 r2_0 r1_1 r2_1 ...
bqtools pipe sample.cbq -p 4 -X 'cat {R1} {R2} | wc -l'
```

!!! note
    `{R1}` and `{R2}` written apart each expand to their own list (`--in1 {R1} --in2 {R2}`). `{n}` is not expanded by `-X`.

`pipe` exits non-zero if any spawned command fails, so it works under `set -e` and in workflow managers.

## Manage FIFOs yourself

Without `-x`/`-X`, `pipe` creates the FIFOs and blocks until they are fully read.

```bash
# FIFOs: fifo_0.fq .. fifo_3.fq
bqtools pipe reads.cbq -p 4 -b fifo &
sleep 1  # let the FIFOs appear before globbing

ls fifo_*.fq | xargs -P 4 -I {} sh -c 'wc -l < "$1" > "$1.lines"' _ {}
wait
```

!!! warning
    Every FIFO must be read to the end or `pipe` blocks forever. Prefer `-x`/`-X`, which spawn the readers for you.
