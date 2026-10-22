//! End-to-end smoke test of the `piscem` binary: build a tiny index, map
//! bulk and single-cell reads, and check the output contract that downstream
//! tools (alevin-fry, piscem-infer, simpleaf) depend on.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

fn piscem() -> Command {
    Command::new(env!("CARGO_BIN_EXE_piscem"))
}

/// Deterministic xorshift so the test needs no extra dependencies.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn seq(&mut self, n: usize) -> String {
        (0..n)
            .map(|_| b"ACGT"[(self.next() % 4) as usize] as char)
            .collect()
    }
}

/// A fresh directory under cargo's per-target scratch space, so the test
/// writes nothing to the source tree or the system temp directory.
fn scratch_dir(name: &str) -> PathBuf {
    let pid = std::process::id();
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("piscem-smoke-{name}-{pid}"));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn run_ok(cmd: &mut Command) {
    let out = cmd.output().expect("failed to spawn piscem");
    assert!(
        out.status.success(),
        "piscem failed: {:?}\nstderr:\n{}",
        cmd,
        String::from_utf8_lossy(&out.stderr)
    );
}

fn assert_rad_v2(path: &Path) {
    let bytes = fs::read(path).unwrap();
    assert!(bytes.len() > 10, "{} is too short", path.display());
    assert_eq!(
        &bytes[..8],
        b"RAD_FILE",
        "{} lacks the RAD v2 magic",
        path.display()
    );
    assert_eq!(bytes[8], 2, "{} is not RAD spec major 2", path.display());
}

fn assert_all_mapped(map_info: &Path, nreads: usize) {
    let v: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(map_info).unwrap()).unwrap();
    for key in ["num_reads", "num_mapped"] {
        assert_eq!(
            v.get(key),
            Some(&serde_json::json!(nreads)),
            "{key} in {}",
            map_info.display()
        );
    }
}

#[test]
fn build_help_lists_mphf_options() {
    let out = piscem().args(["build", "--help"]).output().unwrap();
    assert!(out.status.success(), "`piscem build --help` failed");
    let help = String::from_utf8_lossy(&out.stdout);
    // Match the option's own `--flag <VALUE>` line: the doc text of other
    // options also mentions these flags by name.
    for flag in ["--mphf", "--mphf-offline", "--tmp-dir", "--work-dir"] {
        assert!(
            help.contains(&format!("{flag} <")),
            "`piscem build --help` lacks {flag}"
        );
    }
}

#[test]
fn build_and_map_write_rad_v2() {
    let dir = scratch_dir("e2e");
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let refs: Vec<String> = (0..40).map(|_| rng.seq(1500)).collect();

    let fa = dir.join("ref.fa");
    let mut f = fs::File::create(&fa).unwrap();
    for (i, s) in refs.iter().enumerate() {
        writeln!(f, ">t{i}\n{s}").unwrap();
    }
    drop(f);

    let idx = dir.join("idx/idx");
    let work = dir.join("work");
    run_ok(
        piscem()
            .args(["build", "-t", "2", "-s"])
            .arg(&fa)
            .arg("-o")
            .arg(&idx)
            .arg("-w")
            .arg(&work),
    );
    for ext in [".ctab", ".refinfo", "_ver.json"] {
        let p = PathBuf::from(format!("{}{}", idx.display(), ext));
        assert!(p.exists(), "missing index file {}", p.display());
    }

    // Bulk: exact 100 bp substrings of the references, so every read maps.
    let nreads = 200;
    let reads = dir.join("bulk.fq");
    let mut r = fs::File::create(&reads).unwrap();
    for i in 0..nreads {
        let s = &refs[i % refs.len()];
        let off = (rng.next() % (s.len() as u64 - 100)) as usize;
        let q = "I".repeat(100);
        writeln!(r, "@r{i}\n{}\n+\n{q}", &s[off..off + 100]).unwrap();
    }
    drop(r);
    let stem = dir.join("bulk/sample");
    run_ok(
        piscem()
            .args(["map-bulk", "-t", "2", "-i"])
            .arg(&idx)
            .arg("-r")
            .arg(&reads)
            .arg("-o")
            .arg(&stem),
    );
    assert_rad_v2(&dir.join("bulk/sample.rad"));
    assert_all_mapped(&dir.join("bulk/sample.map_info.json"), nreads);

    // Single-cell (chromium_v3: 16 bp barcode + 12 bp UMI in read 1).
    let (r1p, r2p) = (dir.join("r1.fq"), dir.join("r2.fq"));
    let (mut r1, mut r2) = (
        fs::File::create(&r1p).unwrap(),
        fs::File::create(&r2p).unwrap(),
    );
    for i in 0..nreads {
        let s = &refs[i % refs.len()];
        let off = (rng.next() % (s.len() as u64 - 91)) as usize;
        let bcumi = rng.seq(28);
        writeln!(r1, "@r{i}\n{bcumi}\n+\n{}", "I".repeat(28)).unwrap();
        writeln!(r2, "@r{i}\n{}\n+\n{}", &s[off..off + 91], "I".repeat(91)).unwrap();
    }
    drop((r1, r2));
    let sc = dir.join("sc");
    run_ok(
        piscem()
            .args(["map-sc", "-t", "2", "-g", "chromium_v3", "-i"])
            .arg(&idx)
            .arg("-1")
            .arg(&r1p)
            .arg("-2")
            .arg(&r2p)
            .arg("-o")
            .arg(&sc),
    );
    assert_rad_v2(&sc.join("map.rad"));
    assert_all_mapped(&sc.join("map_info.json"), nreads);
    assert!(sc.join("unmapped_bc_count.bin").exists());

    let _ = fs::remove_dir_all(&dir);
}
