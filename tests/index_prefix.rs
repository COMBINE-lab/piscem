//! Where `piscem build` puts an index.
//!
//! Every file of an index is named after the full output prefix: piscem-rs's
//! files and the signatures are `<prefix>.<ext>`, the wrapper's own are
//! `<prefix>_cfish.*` and `<prefix>_ver.json`. A dot in the prefix's file name
//! (`-o idx/gen.v2`) is part of the name, so two such indexes in one directory
//! do not collide. A prefix that names a directory (`-o idx/`), or that is not
//! UTF-8, is refused before anything is written.

use std::collections::{BTreeMap, BTreeSet};
use std::ffi::OsStr;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// A per-test scratch directory under cargo's per-target scratch space
/// (`CARGO_TARGET_TMPDIR`), so the test writes nothing to the source tree or
/// the system temp directory. It is removed afterwards unless the test
/// failed, so a failure can be inspected.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!(
            "piscem-index-prefix-{}-{}",
            name,
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let s = Scratch(dir);
        s.write_inputs();
        s
    }

    fn path(&self, rel: &str) -> PathBuf {
        self.0.join(rel)
    }

    /// Two unrelated random references (refA has 20 records, refB 30), decoys
    /// that embed refA segments in random flanks (so they yield poison
    /// k-mers), and 100 bp reads cut from each reference.
    fn write_inputs(&self) {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let ref_a: Vec<String> = (0..20).map(|_| rng.seq(1000)).collect();
        let ref_b: Vec<String> = (0..30).map(|_| rng.seq(1000)).collect();
        write_fasta(&self.path("refA.fa"), "A", &ref_a);
        write_fasta(&self.path("refB.fa"), "B", &ref_b);
        let decoys: Vec<String> = ref_a
            .iter()
            .flat_map(|s| (0..900).step_by(300).map(|j| s[j..j + 60].to_owned()))
            .map(|core| format!("{}{}{}", rng.seq(40), core, rng.seq(40)))
            .collect();
        write_fasta(&self.path("decoys.fa"), "D", &decoys);
        write_reads(&self.path("readsA.fq"), &ref_a, NREADS, &mut rng);
        write_reads(&self.path("readsB.fq"), &ref_b, NREADS, &mut rng);
    }

    /// Run `piscem` with `args` in this directory.
    fn piscem<S: AsRef<OsStr> + std::fmt::Debug>(&self, args: &[S]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_piscem"))
            .current_dir(&self.0)
            .env("NO_COLOR", "1")
            .args(args)
            .output()
            .unwrap()
    }

    fn piscem_ok(&self, args: &[&str]) {
        let out = self.piscem(args);
        assert!(
            out.status.success(),
            "piscem {:?} failed:\n{}",
            args,
            String::from_utf8_lossy(&out.stderr)
        );
    }

    /// Path relative to this directory -> contents (`None` for a directory)
    /// of everything below it. Compare two with [`changes`], which names the
    /// paths that differ instead of printing the files.
    fn tree(&self) -> Tree {
        fn walk(root: &Path, dir: &Path, out: &mut Tree) {
            for e in fs::read_dir(dir).unwrap() {
                let p = e.unwrap().path();
                let rel = p.strip_prefix(root).unwrap().to_path_buf();
                if p.is_dir() {
                    out.insert(rel, None);
                    walk(root, &p, out);
                } else {
                    out.insert(rel, Some(fs::read(&p).unwrap()));
                }
            }
        }
        let mut out = BTreeMap::new();
        walk(&self.0, &self.0, &mut out);
        out
    }

    /// Name -> contents of every file in the directory `rel`.
    fn files(&self, rel: &str) -> BTreeMap<String, Vec<u8>> {
        fs::read_dir(self.path(rel))
            .unwrap()
            .map(|e| {
                let e = e.unwrap();
                (
                    e.file_name().to_string_lossy().into_owned(),
                    fs::read(e.path()).unwrap(),
                )
            })
            .collect()
    }

    /// `num_mapped` from the `map_info.json` of the `map-bulk` stem `stem`.
    fn num_mapped(&self, stem: &str) -> u64 {
        let info = fs::read_to_string(self.path(&format!("{stem}.map_info.json"))).unwrap();
        let v: serde_json::Value = serde_json::from_str(&info).unwrap();
        assert_eq!(v["num_reads"], NREADS as u64, "num_reads in {stem}");
        v["num_mapped"].as_u64().unwrap()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
}

const NREADS: usize = 200;

type Tree = BTreeMap<PathBuf, Option<Vec<u8>>>;

/// The paths added, removed or changed from `before` to `after`.
fn changes(before: &Tree, after: &Tree) -> Vec<String> {
    let paths: BTreeSet<&PathBuf> = before.keys().chain(after.keys()).collect();
    paths
        .into_iter()
        .filter_map(|p| match (before.get(p), after.get(p)) {
            (None, _) => Some(format!("added {}", p.display())),
            (_, None) => Some(format!("removed {}", p.display())),
            (Some(b), Some(a)) if a != b => Some(format!("changed {}", p.display())),
            _ => None,
        })
        .collect()
}

/// xorshift64: deterministic sequences without a `rand` dependency.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn seq(&mut self, len: usize) -> String {
        (0..len)
            .map(|_| b"ACGT"[(self.next() >> 32) as usize % 4] as char)
            .collect()
    }
}

fn write_fasta(path: &Path, prefix: &str, seqs: &[String]) {
    let mut s = String::new();
    for (i, seq) in seqs.iter().enumerate() {
        s.push_str(&format!(">{prefix}{i}\n{seq}\n"));
    }
    fs::write(path, s).unwrap();
}

/// `n` exact 100 bp substrings of `refs`, so every read maps to its reference.
fn write_reads(path: &Path, refs: &[String], n: usize, rng: &mut Rng) {
    let mut s = String::new();
    for i in 0..n {
        let r = &refs[i % refs.len()];
        let off = (rng.next() % (r.len() as u64 - 100)) as usize;
        s.push_str(&format!(
            "@r{i}\n{}\n+\n{}\n",
            &r[off..off + 100],
            "I".repeat(100)
        ));
    }
    fs::write(path, s).unwrap();
}

/// Every file a build with all the optional files enabled writes for a prefix
/// whose file name is `name`.
fn expected_files(name: &str) -> BTreeSet<String> {
    [
        ".ssi",
        ".ssi.mphf",
        ".ctab",
        ".refinfo",
        ".ectab",
        ".tdct",
        ".tct",
        ".poison",
        ".poison.json",
        ".sigs.json",
        "_cfish.json",
        "_cfish.cf_seg",
        "_cfish.cf_seq",
        "_ver.json",
    ]
    .iter()
    .map(|suffix| format!("{name}{suffix}"))
    .collect()
}

/// `build -o <prefix>` fails, names the prefix (as the error prints it,
/// `shown`), and writes nothing at all: no signatures, no cf1-rs files, no
/// output directory, and no work directory (there is no `-w`, so the default
/// `./workdir.noindex` would be created in the scratch directory).
fn assert_build_refused(s: &Scratch, prefix: &OsStr, shown: &str) {
    let before = s.tree();
    let mut args: Vec<&OsStr> = [
        "build",
        "-s",
        "refA.fa",
        "--decoy-paths",
        "decoys.fa",
        "-t",
        "1",
        "-o",
    ]
    .into_iter()
    .map(OsStr::new)
    .collect();
    args.push(prefix);
    let out = s.piscem(&args);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !out.status.success(),
        "`build -o {shown}` should fail:\n{stderr}"
    );
    assert!(
        stderr.contains("invalid index prefix") && stderr.contains(&format!("\"{shown}\"")),
        "`build -o {shown}`: the error should name the prefix:\n{stderr}"
    );
    let changed = changes(&before, &s.tree());
    assert!(
        changed.is_empty(),
        "`build -o {shown}` wrote something: {changed:#?}"
    );
}

/// `-o <dir>/` names a directory, not an index, so it is refused before
/// anything is written, whether or not the directory exists.
#[test]
fn build_refuses_a_directory_prefix_before_writing_anything() {
    let s = Scratch::new("dir");
    fs::create_dir(s.path("existing")).unwrap();
    for prefix in ["existing/", "missing/", "missing/deeper/", "."] {
        assert_build_refused(&s, OsStr::new(prefix), prefix);
    }
}

/// A prefix that is not UTF-8 is refused the same way: cf1-rs cannot write
/// the graph to it (it panics), and `map-* -i` could not load the index.
#[cfg(unix)]
#[test]
fn build_refuses_a_non_utf8_prefix_before_writing_anything() {
    use std::os::unix::ffi::OsStrExt;
    let s = Scratch::new("utf8");
    assert_build_refused(&s, OsStr::from_bytes(b"idx/g\xffn"), r"idx/g\xFFn");
}

/// An undotted prefix (`idx/gen`, like simpleaf's `<dir>/piscem_idx`) and a
/// dotted one sharing its stem (`idx/gen.v2`) in one directory: each build
/// writes exactly the `<prefix>.*` / `<prefix>_*` files, neither touches the
/// other's, and each index maps the reads of its own reference.
#[test]
fn undotted_and_dotted_prefixes_name_every_file_after_the_full_prefix() {
    let s = Scratch::new("names");
    let build = |reference: &str, prefix: &str| {
        s.piscem_ok(&[
            "build",
            "-s",
            reference,
            "--decoy-paths",
            "decoys.fa",
            "--dict",
            "tiny",
            "--keep-intermediate-dbg",
            "-t",
            "1",
            "-w",
            "work",
            "-o",
            prefix,
        ]);
    };

    build("refA.fa", "idx/gen");
    let gen_files = s.files("idx");
    let names: BTreeSet<String> = gen_files.keys().cloned().collect();
    assert_eq!(names, expected_files("gen"));

    build("refB.fa", "idx/gen.v2");
    let mut both = expected_files("gen");
    both.extend(expected_files("gen.v2"));
    let after = s.files("idx");
    let names: BTreeSet<String> = after.keys().cloned().collect();
    assert_eq!(names, both);
    for (name, bytes) in &gen_files {
        assert!(
            after.get(name) == Some(bytes),
            "building idx/gen.v2 changed {name}"
        );
    }

    for (index, reads, stem, mapped) in [
        ("idx/gen", "readsA.fq", "map/A_on_gen", NREADS as u64),
        ("idx/gen", "readsB.fq", "map/B_on_gen", 0),
        ("idx/gen.v2", "readsB.fq", "map/B_on_gen_v2", NREADS as u64),
        ("idx/gen.v2", "readsA.fq", "map/A_on_gen_v2", 0),
    ] {
        s.piscem_ok(&["map-bulk", "-t", "1", "-i", index, "-r", reads, "-o", stem]);
        assert_eq!(s.num_mapped(stem), mapped, "{reads} against {index}");
    }
}
