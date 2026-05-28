# rsomics-bed-cluster

Cluster overlapping BED intervals and append a cluster ID — `bedtools cluster` equivalent.

## Usage

```
rsomics-bed-cluster [OPTIONS]

Options:
  -i, --input <FILE>   BED input file (default: stdin)
  -o, --out <FILE>     Output file (default: stdout)
  -d, --dist <N>       Max gap (bp) between intervals in the same cluster [default: 0]
  -s, --strand         Strand-specific clustering (requires ≥6-column BED)
  -h, --help           Print help
  -V, --version        Print version
```

## Description

Reads a coordinate-sorted BED file and appends a 1-based cluster-ID column.
Two records belong to the same cluster when they overlap (or, with `-d N`,
when their gap is ≤ N bp).  The cluster ID increments each time a new cluster
starts.

Input must be sorted by `(chrom, start)`.  With `-s`, the input must also be
sorted by strand.

## Example

```
$ printf 'chr1\t1\t10\nchr1\t5\t15\nchr1\t20\t30\nchr2\t1\t5\n' \
    | rsomics-bed-cluster
chr1	1	10	1
chr1	5	15	1
chr1	20	30	2
chr2	1	5	3
```

## Install

```sh
cargo install rsomics-bed-cluster
```

## Origin

Rust reimplementation of `bedtools cluster`.  Algorithm informed by the
[bedtools2 source](https://github.com/arq5x/bedtools2) (MIT License) and the
[bedtools documentation](https://bedtools.readthedocs.io/en/latest/content/tools/cluster.html).

License: MIT OR Apache-2.0.
Upstream credit: bedtools2 <https://github.com/arq5x/bedtools2> (MIT License).
