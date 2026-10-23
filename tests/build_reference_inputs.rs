//! `piscem build` reference input handling: `--polya-clip-length` reaches the
//! index, and several `--ref-dirs` index only their FASTA/FASTQ files.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const CLIP: usize = 10;

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
    let d = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
        "piscem-refinputs-{}-{}",
        name,
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn fasta(records: &[(String, String)]) -> String {
    let mut s = String::new();
    for (name, seq) in records {
        s.push_str(&format!(">{} some description\n", name));
        for line in seq.as_bytes().chunks(60) {
            s.push_str(std::str::from_utf8(line).unwrap());
            s.push('\n');
        }
    }
    s
}

fn write_fasta(path: &Path, records: &[(String, String)]) {
    fs::write(path, fasta(records)).unwrap();
}

/// `data` as a gzip file holding one stored (uncompressed) deflate block, so
/// the test needs no compression crate.
fn gzip_stored(data: &[u8]) -> Vec<u8> {
    let len = u16::try_from(data.len()).unwrap();
    let crc = !data.iter().fold(!0u32, |mut c, &b| {
        c ^= u32::from(b);
        for _ in 0..8 {
            c = if c & 1 == 1 {
                (c >> 1) ^ 0xedb8_8320
            } else {
                c >> 1
            };
        }
        c
    });
    let mut gz = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff];
    gz.push(1); // final block, stored
    gz.extend(len.to_le_bytes());
    gz.extend((!len).to_le_bytes());
    gz.extend(data);
    gz.extend(crc.to_le_bytes());
    gz.extend((data.len() as u32).to_le_bytes());
    gz
}

/// `n` random references named `<prefix><i>`; every third one ends in a
/// 40-base polyA tail, the others end in a base other than `A`.
fn references(rng: &mut Rng, prefix: &str, n: usize) -> Vec<(String, String)> {
    (0..n)
        .map(|i| {
            let len = 200 + (rng.next() % 600) as usize;
            let mut seq = rng.seq(len).trim_end_matches('A').to_string();
            if i % 3 == 0 {
                seq.push_str(&"A".repeat(40));
            }
            (format!("{}{}", prefix, i), seq)
        })
        .collect()
}

/// Run `piscem build -t 1` (deterministic output) with `args`, writing the
/// index to `<dir>/<name>`; returns the index prefix and the process output.
fn build(dir: &Path, name: &str, args: &[&str]) -> (PathBuf, Output) {
    let prefix = dir.join(name);
    let out = piscem()
        .arg("build")
        .args(["-t", "1", "-k", "31", "-m", "19"])
        .arg("-o")
        .arg(&prefix)
        .arg("-w")
        .arg(dir.join(format!("{}_work", name)))
        .args(args)
        .output()
        .expect("failed to spawn piscem");
    (prefix, out)
}

fn build_ok(dir: &Path, name: &str, args: &[&str]) -> PathBuf {
    let (prefix, out) = build(dir, name, args);
    assert!(
        out.status.success(),
        "piscem build {:?} failed\nstderr:\n{}",
        args,
        String::from_utf8_lossy(&out.stderr)
    );
    prefix
}

fn index_file(prefix: &Path, ext: &str) -> PathBuf {
    PathBuf::from(format!("{}.{}", prefix.display(), ext))
}

fn read_index_file(prefix: &Path, ext: &str) -> Vec<u8> {
    let p = index_file(prefix, ext);
    fs::read(&p).unwrap_or_else(|e| panic!("reading {}: {}", p.display(), e))
}

/// Reference names and lengths from `<prefix>.refinfo` (piscem-rs `PRFINF01`).
fn ref_lengths(prefix: &Path) -> HashMap<String, u64> {
    let b = read_index_file(prefix, "refinfo");
    assert_eq!(&b[..8], b"PRFINF01");
    let u64_at = |i: usize| u64::from_le_bytes(b[i..i + 8].try_into().unwrap());
    let n = u64_at(8) as usize;
    let mut pos = 16;
    let mut names = Vec::with_capacity(n);
    for _ in 0..n {
        let len = u32::from_le_bytes(b[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;
        names.push(String::from_utf8(b[pos..pos + len].to_vec()).unwrap());
        pos += len;
    }
    names
        .into_iter()
        .enumerate()
        .map(|(i, name)| (name, u64_at(pos + 8 * i)))
        .collect()
}

fn signatures(prefix: &Path) -> serde_json::Value {
    serde_json::from_slice(&read_index_file(prefix, "sigs.json")).unwrap()
}

/// The index files that depend on the reference sequences.
const INDEX_EXTS: [&str; 5] = ["refinfo", "ctab", "ectab", "ssi", "ssi.mphf"];

#[test]
fn polya_clip_length_reaches_the_index() {
    let dir = scratch_dir("polya");
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let refs_a = references(&mut rng, "a", 60);
    let refs_b = references(&mut rng, "b", 30);
    let clip = |refs: &[(String, String)]| -> Vec<(String, String)> {
        refs.iter()
            .map(|(n, s)| (n.clone(), s.trim_end_matches('A').to_string()))
            .collect()
    };
    write_fasta(&dir.join("a.fa"), &refs_a);
    write_fasta(&dir.join("b.fa"), &refs_b);
    write_fasta(&dir.join("a_clipped.fa"), &clip(&refs_a));
    write_fasta(&dir.join("b_clipped.fa"), &clip(&refs_b));
    let (a, b) = (dir.join("a.fa"), dir.join("b.fa"));
    let inputs = format!("{},{}", a.display(), b.display());
    let clipped_inputs = format!(
        "{},{}",
        dir.join("a_clipped.fa").display(),
        dir.join("b_clipped.fa").display()
    );
    let clip_arg = CLIP.to_string();

    let plain = build_ok(&dir, "plain", &["-s", &inputs]);
    let clipped = build_ok(
        &dir,
        "clipped",
        &["-s", &inputs, "--polya-clip-length", &clip_arg],
    );
    let expected = build_ok(&dir, "expected", &["-s", &clipped_inputs]);

    // The flag yields the index of the clipped references...
    for ext in INDEX_EXTS {
        assert_eq!(
            read_index_file(&clipped, ext),
            read_index_file(&expected, ext),
            "{} differs from the index of the pre-clipped references",
            ext
        );
    }
    let plain_lens = ref_lengths(&plain);
    let clipped_lens = ref_lengths(&clipped);
    assert_eq!(plain_lens.len(), refs_a.len() + refs_b.len());
    let mut n_tailed = 0;
    for (name, seq) in refs_a.iter().chain(&refs_b) {
        let tail = seq.len() - seq.trim_end_matches('A').len();
        assert_eq!(plain_lens[name], seq.len() as u64, "{}", name);
        let want = if tail >= CLIP {
            seq.len() - tail
        } else {
            seq.len()
        };
        assert_eq!(clipped_lens[name], want as u64, "{}", name);
        n_tailed += usize::from(tail >= CLIP);
    }
    assert!(n_tailed > 0);

    // ...while the signatures stay those of the unclipped input.
    let (sp, sc) = (signatures(&plain), signatures(&clipped));
    for key in ["sha256_names", "sha256_seqs", "sha512_names", "sha512_seqs"] {
        assert_eq!(sp[key], sc[key], "{}", key);
    }
    assert_eq!(sc["num_trimmed_polya"], n_tailed);

    // The clipped FASTA is an intermediate, removed by default.
    assert!(!index_file(&clipped, "fa.gz").exists());
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn polya_clip_length_keeps_clipped_reference_with_keep_intermediate_dbg() {
    let dir = scratch_dir("polya-keep");
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    write_fasta(&dir.join("a.fa"), &references(&mut rng, "a", 30));
    let a = dir.join("a.fa");
    let clip_arg = CLIP.to_string();
    let idx = build_ok(
        &dir,
        "kept",
        &[
            "-s",
            a.to_str().unwrap(),
            "--polya-clip-length",
            &clip_arg,
            "--keep-intermediate-dbg",
        ],
    );
    let clipped_fa = index_file(&idx, "fa.gz");
    assert!(clipped_fa.exists(), "{} was removed", clipped_fa.display());
    // gzip magic: this is prepare_fasta's clipped output, not a leftover.
    assert_eq!(&fs::read(&clipped_fa).unwrap()[..2], &[0x1f, 0x8b]);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn polya_clip_length_refuses_to_overwrite_its_input() {
    let dir = scratch_dir("polya-input");
    let mut rng = Rng(0x1234_5678_9abc_def1);
    // `-s ref.fa.gz -o ref` would make prepare_fasta truncate the input
    // before reading it. The guard runs before anything is read, so the
    // content need not be gzipped.
    let input = dir.join("ref.fa.gz");
    write_fasta(&input, &references(&mut rng, "a", 10));
    let before = fs::read(&input).unwrap();
    let clip_arg = CLIP.to_string();
    let (_, out) = build(
        &dir,
        "ref",
        &[
            "-s",
            input.to_str().unwrap(),
            "--polya-clip-length",
            &clip_arg,
        ],
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success(), "build succeeded:\n{}", stderr);
    assert!(
        stderr.contains("choose another --output prefix"),
        "{}",
        stderr
    );
    assert_eq!(fs::read(&input).unwrap(), before, "the input was modified");

    // Likewise for a decoy file, which the build would also have removed.
    let refs = dir.join("tx.fa");
    write_fasta(&refs, &references(&mut rng, "t", 10));
    let decoy = dir.join("idx.fa.gz");
    write_fasta(&decoy, &references(&mut rng, "d", 5));
    let mut cases = vec![(decoy, "idx")];
    // On Unix, also a decoy file that is `<output>.fa.gz` under another name
    // (a hard link, so its canonical path differs).
    #[cfg(unix)]
    {
        let decoy = dir.join("dec.fa");
        write_fasta(&decoy, &references(&mut rng, "h", 5));
        fs::hard_link(&decoy, dir.join("hl.fa.gz")).unwrap();
        cases.push((decoy, "hl"));
    }
    for (decoy, name) in cases {
        let before = fs::read(&decoy).unwrap();
        let (prefix, out) = build(
            &dir,
            name,
            &[
                "-s",
                refs.to_str().unwrap(),
                "--decoy-paths",
                decoy.to_str().unwrap(),
                "--polya-clip-length",
                &clip_arg,
            ],
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            !out.status.success(),
            "-o {} build succeeded:\n{}",
            name,
            stderr
        );
        assert!(
            stderr.contains("which is the decoy file")
                && stderr.contains("choose another --output prefix"),
            "{}",
            stderr
        );
        assert_eq!(
            fs::read(&decoy).unwrap(),
            before,
            "the decoy file {} was modified",
            decoy.display()
        );
        assert!(!index_file(&prefix, "refinfo").exists());
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn polya_clip_length_is_rejected_with_ref_lists_and_ref_dirs() {
    let dir = scratch_dir("polya-reject");
    let mut rng = Rng(0x0bad_cafe_f00d_d00d);
    let refs = dir.join("refs");
    fs::create_dir_all(&refs).unwrap();
    write_fasta(&refs.join("a.fa"), &references(&mut rng, "a", 10));
    let list = dir.join("list.txt");
    fs::write(&list, format!("{}\n", refs.join("a.fa").display())).unwrap();
    let clip_arg = CLIP.to_string();
    for (flag, path) in [("-l", &list), ("-d", &refs)] {
        let (prefix, out) = build(
            &dir,
            "idx",
            &[
                flag,
                path.to_str().unwrap(),
                "--polya-clip-length",
                &clip_arg,
            ],
        );
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "{} build succeeded", flag);
        assert!(
            stderr.contains("cannot be used with")
                && stderr.contains("--polya-clip-length")
                && !stderr.contains("panicked"),
            "{}: {}",
            flag,
            stderr
        );
        assert!(!index_file(&prefix, "refinfo").exists());
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn multiple_ref_dirs_index_only_sequence_files() {
    let dir = scratch_dir("refdirs");
    let mut rng = Rng(0x5851_f42d_4c95_7f2d);
    let dirs: Vec<PathBuf> = ["d1", "d2", "d3"].iter().map(|d| dir.join(d)).collect();
    for d in &dirs {
        fs::create_dir_all(d).unwrap();
    }
    // One FASTA per directory, so the order in which each is listed does not
    // matter. Their names do not either: a FASTA with no extension, or a
    // gzipped one, is indexed as it was before.
    let a = dirs[0].join("a.fa");
    let c = dirs[1].join("chr1");
    let b = dirs[2].join("b.fa.gz");
    write_fasta(&a, &references(&mut rng, "a", 40));
    write_fasta(&c, &references(&mut rng, "c", 20));
    let b_fasta = fasta(&references(&mut rng, "b", 20));
    fs::write(&b, gzip_stored(b_fasta.as_bytes())).unwrap();
    // None is FASTA/FASTQ; handed to cf1-rs, each panics the build.
    let skipped = [
        dirs[0].join("a.fa.fai"),
        dirs[0].join("empty.fa"),
        dirs[1].join("README.txt"),
        dirs[2].join("genes.gtf.gz"),
    ];
    fs::write(&skipped[0], "a0\t300\t40\t60\t61\n").unwrap();
    fs::write(&skipped[1], "").unwrap();
    fs::write(&skipped[2], "these are references\n").unwrap();
    let gtf = "c0\tsrc\texon\t1\t100\t.\t+\t.\tgene_id \"g0\";\n";
    fs::write(&skipped[3], gzip_stored(gtf.as_bytes())).unwrap();

    let dir_list = dirs
        .iter()
        .map(|d| d.display().to_string())
        .collect::<Vec<_>>()
        .join(",");
    let (from_dirs, out) = build(&dir, "dirs", &["-d", &dir_list]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "build failed:\n{}", stderr);
    for f in &skipped {
        assert!(
            stderr.contains(&format!(
                "skipping {}: empty, or not a FASTA/FASTQ file",
                f.display()
            )),
            "no warning for {}:\n{}",
            f.display(),
            stderr
        );
    }
    // ...and only for them.
    assert_eq!(
        stderr.matches(": empty, or not a FASTA/FASTQ file").count(),
        skipped.len(),
        "{}",
        stderr
    );

    // Same index as naming the three FASTA files directly.
    let seqs = format!("{},{},{}", a.display(), c.display(), b.display());
    let from_seqs = build_ok(&dir, "seqs", &["-s", &seqs]);
    for ext in INDEX_EXTS {
        assert_eq!(
            read_index_file(&from_dirs, ext),
            read_index_file(&from_seqs, ext),
            "{} differs",
            ext
        );
    }
    assert_eq!(ref_lengths(&from_dirs).len(), 80);
    let _ = fs::remove_dir_all(&dir);
}
