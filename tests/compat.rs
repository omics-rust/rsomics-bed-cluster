//! Golden tests (always) + bedtools compat tests (when bedtools is available).

use rsomics_bed_cluster::cluster;
use std::io::Cursor;
use std::path::Path;
use std::process::Command;

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn read_golden(name: &str) -> Vec<u8> {
    std::fs::read(golden(name)).unwrap_or_else(|e| panic!("cannot read golden {name}: {e}"))
}

fn bedtools_available() -> bool {
    Command::new("bedtools").arg("--version").output().is_ok()
}

// ── Golden tests ────────────────────────────────────────────────────────────

#[test]
fn golden_basic_cluster() {
    let input = read_golden("input.bed");
    let expected = read_golden("basic_cluster.expected");
    let mut got = Vec::new();
    cluster(Cursor::new(&input), &mut got, 0, false).unwrap();
    assert_eq!(
        String::from_utf8(got).unwrap(),
        String::from_utf8(expected).unwrap(),
        "basic_cluster golden mismatch"
    );
}

#[test]
fn golden_dist5() {
    let input = read_golden("input.bed");
    let expected = read_golden("dist5.expected");
    let mut got = Vec::new();
    cluster(Cursor::new(&input), &mut got, 5, false).unwrap();
    assert_eq!(
        String::from_utf8(got).unwrap(),
        String::from_utf8(expected).unwrap(),
        "dist5 golden mismatch"
    );
}

// ── bedtools compat tests ────────────────────────────────────────────────────

#[test]
fn compat_basic_cluster() {
    if !bedtools_available() {
        eprintln!("bedtools not found, skipping compat test");
        return;
    }
    let input_path = golden("input.bed");
    let bt_out = Command::new("bedtools")
        .args(["cluster", "-i"])
        .arg(&input_path)
        .output()
        .expect("bedtools cluster failed");
    assert!(bt_out.status.success(), "bedtools cluster exited non-zero");

    let input = read_golden("input.bed");
    let mut ours = Vec::new();
    cluster(Cursor::new(&input), &mut ours, 0, false).unwrap();

    assert_eq!(
        String::from_utf8(ours).unwrap(),
        String::from_utf8(bt_out.stdout).unwrap(),
        "compat_basic_cluster: output differs from bedtools"
    );
}

#[test]
fn compat_dist5() {
    if !bedtools_available() {
        eprintln!("bedtools not found, skipping compat test");
        return;
    }
    let input_path = golden("input.bed");
    let bt_out = Command::new("bedtools")
        .args(["cluster", "-i"])
        .arg(&input_path)
        .args(["-d", "5"])
        .output()
        .expect("bedtools cluster -d 5 failed");
    assert!(
        bt_out.status.success(),
        "bedtools cluster -d 5 exited non-zero"
    );

    let input = read_golden("input.bed");
    let mut ours = Vec::new();
    cluster(Cursor::new(&input), &mut ours, 5, false).unwrap();

    assert_eq!(
        String::from_utf8(ours).unwrap(),
        String::from_utf8(bt_out.stdout).unwrap(),
        "compat_dist5: output differs from bedtools"
    );
}
