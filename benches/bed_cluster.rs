use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use rsomics_bed_cluster::cluster;
use std::io::sink;

fn make_fixture(n: usize) -> Vec<u8> {
    use std::fmt::Write as _;
    let mut s = String::with_capacity(n * 32);
    let chroms = ["chr1", "chr2", "chr3", "chr4", "chr5"];
    let mut pos: u64 = 0;
    for i in 0..n {
        let chrom = chroms[i % chroms.len()];
        if i % chroms.len() == 0 {
            pos = 0;
        }
        let start = pos;
        let end = pos + 100;
        let _ = writeln!(s, "{chrom}\t{start}\t{end}");
        // Alternate between overlapping (+50) and non-overlapping (+200).
        pos += if i % 3 == 0 { 50 } else { 200 };
    }
    s.into_bytes()
}

fn bench_cluster(c: &mut Criterion) {
    let fixture = make_fixture(100_000);
    let mut group = c.benchmark_group("cluster");
    group.throughput(Throughput::Bytes(fixture.len() as u64));
    group.bench_function("100k_records", |b| {
        b.iter(|| {
            cluster(fixture.as_slice(), sink(), 0, false).unwrap();
        });
    });
    group.finish();
}

criterion_group!(benches, bench_cluster);
criterion_main!(benches);
