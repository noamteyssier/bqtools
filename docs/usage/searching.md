# Searching

Find, count, and demultiplex records by sequence or header with `grep` and `split`.

## Match records

Multiple patterns must all match (AND) unless [`--or-logic`](../commands/grep.md#grep--or-logic) is set. Matches go to stdout as TSV, or to `-o` in the format of its extension.

```bash
# fixed string in either mate
bqtools grep reads.cbq ACGTACGT

# regex
bqtools grep reads.cbq "ACGT[AC]TCCA"

# both patterns (AND)
bqtools grep reads.cbq ACGTACGT "AG(TTTT|CCCC)A"

# either pattern (OR)
bqtools grep reads.cbq ACGTACGT "AG(TTTT|CCCC)A" --or-logic

# records without a match
bqtools grep reads.cbq ACGTACGT -v -o nomatch.fastq.gz
```

## Target a mate or region

```bash
# -r searches R1, -R searches R2
bqtools grep sample.cbq -r ACGTACGT -R TTGGCCAA

# only bases 0..12 of each read (0-based, end-exclusive)
bqtools grep sample.cbq ACGTACGT --range ..12

# match and output R1 only
bqtools grep sample.cbq ACGTACGT -m 1 -o hits_R1.fq
```

!!! tip
    Anchoring a barcode with [`--range`](../commands/grep.md#grep--range) is faster and avoids hits deeper in the read.

## Reverse complement and headers

```bash
# also search the other strand: rerun with --rc (fixed ACGT patterns only)
bqtools grep reads.cbq ACGTACGT --rc

# match header text; -x treats it as a literal
bqtools grep reads.cbq "sample_alpha" --header -x
```

## Pattern files

[`--file`](../commands/grep.md#grep--file) accepts plain text (one per line), FASTA, or TSV (`alias<TAB>pattern`). FASTA headers and TSV aliases name the pattern in output. File patterns always use OR logic.

```bash
# either mate
bqtools grep reads.cbq --file patterns.txt

# R1 only / R2 only
bqtools grep sample.cbq --sfile barcodes.fa
bqtools grep sample.cbq --xfile guides.tsv
```

!!! note
    All-uppercase `ACGT` patterns use Aho-Corasick automatically. For other literals (lowercase, header text) pass `-x`.

## Count

Counting never writes records, so it can't be combined with `-o`/`-p`.

```bash
# number of matching records
bqtools grep reads.cbq ACGTACGT -C

# count, total, fraction
bqtools grep reads.cbq ACGTACGT -F

# per-pattern counts: name, count, fraction
bqtools grep reads.cbq --file patterns.tsv -P
```

!!! tip
    `-P` counts every barcode or guide in one pass over the file. A record counts at most once per pattern but may count toward several patterns.

## Fuzzy matching

Needs the `fuzzy` feature. All patterns in a set must be the same length, and regex is not supported.

```bash
# up to 1 edit (default)
bqtools grep reads.cbq ACGTACGT -z

# up to 2 edits, only inexact hits
bqtools grep reads.cbq ACGTACGT -z -k 2 -i

# reject any match containing an N
bqtools grep reads.cbq ACGTACGT -z --max-n-frac 0.0

# fuzzy per-pattern counts
bqtools grep reads.cbq --file patterns.fa -P -z
```

## Demultiplex with split

`split` writes each record to `<alias>.<ext>` in the input's variant. Records matching no pattern, or more than one, go to `unmatched`.

```bash
# one file per alias under ./split_outs
bqtools split sample.cbq --sfile patterns.tsv

# custom directory, keep only files with 100+ records
bqtools split sample.cbq --sfile patterns.tsv --basepath by_sample --min-records 100

# fuzzy, and don't write unmatched records
bqtools split sample.cbq --sfile patterns.tsv -z -k 1 --skip-unmatched
```

!!! note
    Empty outputs are removed by default (`--min-records 1`). Use `--min-records 0` to keep them.
