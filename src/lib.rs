//! Cluster overlapping BED intervals — `bedtools cluster` equivalent.
//!
//! Reads a coordinate-sorted BED stream and appends a 1-based cluster-ID
//! column to each record.  Two records belong to the same cluster when
//! they overlap (or, with `-d`, when their gap is ≤ `d` bp).
//!
//! Requirements (matching bedtools):
//! - Input sorted by (chrom, start).
//! - With `-s` (strand-specific), input must additionally be sorted by strand.
//! - Chromosomes that re-appear after another chromosome was seen are a
//!   sort violation and terminate with a loud error.

use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use rsomics_common::{Context, Result, RsomicsError};

/// Annotate a sorted BED stream `r` with cluster IDs, writing to `w`.
///
/// - `dist`: extend each interval by this many bases when checking overlap
///   (intervals with a gap ≤ `dist` join the same cluster). `0` = require
///   actual overlap.
/// - `strand_specific`: only cluster intervals on the same strand; the input
///   must be sorted by (chrom, strand, start).
pub fn cluster<R: Read, W: Write>(r: R, w: W, dist: i64, strand_specific: bool) -> Result<()> {
    let mut rdr = BufReader::new(r);
    let mut bw = BufWriter::new(w);
    let mut line: Vec<u8> = Vec::with_capacity(256);

    // Current cluster state.
    let mut cluster_id: u64 = 0;
    let mut cur_chrom: Vec<u8> = Vec::with_capacity(32);
    let mut cur_strand: u8 = b'.'; // only meaningful when strand_specific
    let mut cur_end: i64 = 0; // max end seen in the current cluster

    // Chromosomes already closed (for sort-violation detection).
    let mut closed: Vec<Vec<u8>> = Vec::new();

    let mut have_record = false;
    let mut lineno: usize = 0;

    loop {
        line.clear();
        if rdr.read_until(b'\n', &mut line).map_err(RsomicsError::Io)? == 0 {
            break;
        }
        lineno += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.is_empty()
            || line[0] == b'#'
            || line.starts_with(b"track")
            || line.starts_with(b"browser")
        {
            continue;
        }

        let (chrom, start, end, strand) = parse_bed3_strand(&line, strand_specific)
            .map_err(|e| RsomicsError::InvalidInput(format!("BED line {lineno}: {e}")))?;

        let same_group = if !have_record {
            false
        } else if chrom != cur_chrom.as_slice() {
            // New chromosome.
            if closed.iter().any(|c| c.as_slice() == chrom) {
                return Err(RsomicsError::InvalidInput(format!(
                    "BED line {lineno}: chromosome {} reappears after close — \
                     sort with `sort -k1,1 -k2,2n` first",
                    String::from_utf8_lossy(chrom)
                )));
            }
            closed.push(cur_chrom.clone());
            false
        } else if strand_specific && strand != cur_strand {
            // Same chromosome, different strand.
            false
        } else {
            // Same chromosome (and same strand if strand_specific).
            // Cluster if start <= cur_end + dist (gap ≤ dist).
            start <= cur_end + dist
        };

        if same_group {
            // Extend the cluster's reach.
            if end > cur_end {
                cur_end = end;
            }
        } else {
            // Open a new cluster.
            cluster_id += 1;
            cur_chrom.clear();
            cur_chrom.extend_from_slice(chrom);
            cur_strand = strand;
            cur_end = end;
        }
        have_record = true;

        bw.write_all(&line).rs_context("writing cluster record")?;
        write!(bw, "\t{cluster_id}").rs_context("writing cluster id")?;
        bw.write_all(b"\n").rs_context("writing newline")?;
    }
    bw.flush().map_err(RsomicsError::Io)?;
    Ok(())
}

/// Parse chrom, start, end, and (optionally) strand from a BED byte slice.
/// Returns `(chrom_bytes, start_i64, end_i64, strand_byte)`.
fn parse_bed3_strand(
    s: &[u8],
    need_strand: bool,
) -> std::result::Result<(&[u8], i64, i64, u8), String> {
    let mut it = s.splitn(7, |&c| c == b'\t');
    let chrom = it.next().ok_or("missing chrom")?;
    let start = parse_i64(it.next().ok_or("missing start")?)?;
    let end = parse_i64(it.next().ok_or("missing end")?)?;
    if start >= end {
        return Err(format!(
            "empty or inverted interval: start={start} >= end={end}"
        ));
    }
    let strand = if need_strand {
        // col 4 = name, col 5 = score, col 6 = strand
        let _name = it.next(); // col 4
        let _score = it.next(); // col 5
        let s = it
            .next()
            .ok_or("strand-specific mode requires ≥6 columns")?;
        match s {
            b"+" => b'+',
            b"-" => b'-',
            _ => b'.',
        }
    } else {
        b'.'
    };
    Ok((chrom, start, end, strand))
}

fn parse_i64(b: &[u8]) -> std::result::Result<i64, String> {
    if b.is_empty() {
        return Err("empty integer field".into());
    }
    let (neg, digits) = if b[0] == b'-' {
        (true, &b[1..])
    } else {
        (false, b)
    };
    let mut n: i64 = 0;
    for &c in digits {
        let d = c.wrapping_sub(b'0');
        if d > 9 {
            return Err(format!("bad integer {:?}", String::from_utf8_lossy(b)));
        }
        n = n
            .checked_mul(10)
            .and_then(|n| n.checked_add(i64::from(d)))
            .ok_or_else(|| format!("integer overflows i64: {:?}", String::from_utf8_lossy(b)))?;
    }
    Ok(if neg { -n } else { n })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn run(input: &str, dist: i64, strand_specific: bool) -> String {
        let mut out = Vec::new();
        cluster(Cursor::new(input), &mut out, dist, strand_specific).unwrap();
        String::from_utf8(out).unwrap()
    }

    #[test]
    fn basic_overlap_cluster() {
        // chr1: 1-10 and 5-15 overlap → same cluster; 20-30 no overlap → new cluster.
        // chr2: 1-5 → new cluster.
        let inp = "chr1\t1\t10\nchr1\t5\t15\nchr1\t20\t30\nchr2\t1\t5\n";
        let out = run(inp, 0, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chr1\t1\t10\t1");
        assert_eq!(lines[1], "chr1\t5\t15\t1");
        assert_eq!(lines[2], "chr1\t20\t30\t2");
        assert_eq!(lines[3], "chr2\t1\t5\t3");
    }

    #[test]
    fn dist_extends_cluster() {
        // Gap of 5 between 10 and 15; with dist=5 all three join.
        let inp = "chr1\t1\t10\nchr1\t15\t25\nchr1\t30\t40\n";
        let out = run(inp, 5, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chr1\t1\t10\t1");
        assert_eq!(lines[1], "chr1\t15\t25\t1");
        assert_eq!(lines[2], "chr1\t30\t40\t1");
    }

    #[test]
    fn adjacent_intervals_join_same_cluster() {
        // bedtools cluster: adjacent (touching) intervals join the same cluster even
        // with no gap. start=10 == prev_end=10 satisfies start <= cur_end + 0.
        let inp = "chr1\t1\t10\nchr1\t10\t20\n";
        let out = run(inp, 0, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chr1\t1\t10\t1");
        assert_eq!(lines[1], "chr1\t10\t20\t1");
    }

    #[test]
    fn gap_of_one_is_different_cluster() {
        // A true gap of 1 bp (start=11, prev_end=10) → different cluster with dist=0.
        let inp = "chr1\t1\t10\nchr1\t11\t20\n";
        let out = run(inp, 0, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chr1\t1\t10\t1");
        assert_eq!(lines[1], "chr1\t11\t20\t2");
    }

    #[test]
    fn strand_specific_separates_strands() {
        // Plus and minus strands on same coords → different clusters.
        let inp = "chr1\t1\t10\t.\t.\t+\nchr1\t5\t15\t.\t.\t-\nchr1\t20\t30\t.\t.\t+\n";
        let out = run(inp, 0, true);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chr1\t1\t10\t.\t.\t+\t1");
        assert_eq!(lines[1], "chr1\t5\t15\t.\t.\t-\t2");
        assert_eq!(lines[2], "chr1\t20\t30\t.\t.\t+\t3");
    }

    #[test]
    fn extra_columns_preserved() {
        let inp = "chr1\t0\t10\tfeat1\nchr1\t5\t15\tfeat2\n";
        let out = run(inp, 0, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chr1\t0\t10\tfeat1\t1");
        assert_eq!(lines[1], "chr1\t5\t15\tfeat2\t1");
    }

    #[test]
    fn header_and_blank_skipped() {
        let inp = "# comment\nchr1\t0\t10\n\nchr1\t20\t30\n";
        let out = run(inp, 0, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "chr1\t0\t10\t1");
        assert_eq!(lines[1], "chr1\t20\t30\t2");
    }

    #[test]
    fn cluster_extends_via_bridge() {
        // A overlaps B, B overlaps C → all same cluster (even if A and C don't overlap).
        let inp = "chr1\t1\t10\nchr1\t5\t20\nchr1\t15\t25\n";
        let out = run(inp, 0, false);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines[0], "chr1\t1\t10\t1");
        assert_eq!(lines[1], "chr1\t5\t20\t1");
        assert_eq!(lines[2], "chr1\t15\t25\t1");
    }
}
