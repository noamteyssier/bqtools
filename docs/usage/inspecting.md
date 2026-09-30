# Inspecting

Check what a file holds, confirm it is intact, and run quality control.

## info

```bash
# summary: variant, paired, quality, headers, blocks, records
bqtools info reads.cbq

# record count per file: <count>\t<path>
bqtools info *.cbq --num

# block index (vbq/cbq)
bqtools info reads.cbq --show-index
```

For scripts, use `--num` or `--json`. The default table prints counts with `_` separators.

```bash
# total reads across files
bqtools info *.cbq --num | awk '{s += $1} END {print s}'

# one field with jq
bqtools info reads.cbq --json | jq '.[0].num_records'
```

## verify

[`verify`](../commands/verify.md) prints `<checksum>\t<records>\t<path>`. The checksum ignores record order, so files written by parallel encoders still compare equal.

```bash
bqtools verify reads.cbq

# compare two files
[ "$(bqtools verify a.cbq | cut -f1)" = "$(bqtools verify b.cbq | cut -f1)" ] && echo same
```

```bash
# ignore header differences
bqtools verify reads.cbq --skip-headers

# only R1
bqtools verify sample.cbq -M 1
```

!!! note
    Encoding bq/vbq replaces `N` with a random base by default, so two encodes of the same input differ. Encode with `-p a` for a deterministic result. bq never stores headers.

!!! warning
    The checksum (`xxh3-64`) catches corruption and reordering, not deliberate tampering.

## qc

FastQC-style modules. Writes `summary.md` plus one TSV per module (and per mate) to `-o`.

```bash
# all modules into ./qc
bqtools qc reads.cbq -o qc

# first 100k records, skip the slow modules
bqtools qc reads.cbq -o qc --span ..100000 --skip-dup-levels --skip-overrepresented
```
