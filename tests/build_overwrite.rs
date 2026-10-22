//! `piscem build` into an output prefix that already holds a previous build.
//!
//! piscem-rs loads the optional index files (EC table, Tiny dictionary, poison
//! table) whenever they exist, so a rebuild must never leave one from an
//! earlier build behind: `--overwrite` removes all of them once cf1-rs has
//! read the reference inputs (so a bad one costs nothing), and without it the
//! build refuses to start before writing anything. Every file of an index is named
//! after the full prefix (`<prefix>.<ext>`, `<prefix>_*`), so neither touches
//! another index in the same directory.

use std::collections::BTreeMap;
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
            "piscem-build-overwrite-{}-{}",
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

    /// Two unrelated random references (refA has 20 records, refB 40), decoys
    /// that embed refA segments in random flanks (so they yield poison
    /// k-mers), and a list file naming refB.
    fn write_inputs(&self) {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let ref_a: Vec<String> = (0..20).map(|_| rng.seq(1000)).collect();
        let ref_b: Vec<String> = (0..40).map(|_| rng.seq(1000)).collect();
        write_fasta(&self.path("refA.fa"), "A", &ref_a);
        write_fasta(&self.path("refB.fa"), "B", &ref_b);
        let decoys: Vec<String> = ref_a
            .iter()
            .flat_map(|s| (0..900).step_by(300).map(|j| s[j..j + 60].to_owned()))
            .map(|core| format!("{}{}{}", rng.seq(40), core, rng.seq(40)))
            .collect();
        write_fasta(&self.path("decoys.fa"), "D", &decoys);
        fs::write(
            self.path("listB.txt"),
            format!("{}\n", self.path("refB.fa").display()),
        )
        .unwrap();
    }

    /// Run `piscem build` with `args` (plus `-t 1 -w work`) in this directory.
    fn build(&self, args: &[&str]) -> Output {
        self.build_in(args, "work")
    }

    /// Run `piscem build` with `args` (plus `-t 1 -w <work>`) in this
    /// directory.
    fn build_in(&self, args: &[&str], work: &str) -> Output {
        self.build_command(args, work).output().unwrap()
    }

    /// [`Scratch::build`], with `input` written to the build's standard input.
    #[cfg(unix)]
    fn build_with_stdin(&self, args: &[&str], input: &str) -> Output {
        use std::io::Write;
        use std::process::Stdio;
        let mut child = self
            .build_command(args, "work")
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        // The pipe's buffer holds `input` whether or not the build has
        // started reading; dropping the handle closes the pipe. A write error
        // means the build exited early, which the caller's status check shows.
        let _ = child.stdin.take().unwrap().write_all(input.as_bytes());
        child.wait_with_output().unwrap()
    }

    /// `piscem build` with `args` (plus `-t 1 -w <work>`), to run in this
    /// directory.
    fn build_command(&self, args: &[&str], work: &str) -> Command {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_piscem"));
        cmd.current_dir(&self.0)
            .env("NO_COLOR", "1")
            .arg("build")
            .args(args)
            .args(["-t", "1", "-w", work]);
        cmd
    }

    fn build_ok(&self, args: &[&str]) {
        let out = self.build(args);
        assert!(
            out.status.success(),
            "piscem build {:?} failed:\n{}",
            args,
            String::from_utf8_lossy(&out.stderr)
        );
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

    /// Names of the entries in the directory `rel`, sorted.
    fn names(&self, rel: &str) -> Vec<String> {
        let mut names: Vec<String> = fs::read_dir(self.path(rel))
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        names.sort();
        names
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            let _ = fs::remove_dir_all(&self.0);
        }
    }
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

fn num_records(sigs: &Path) -> serde_json::Value {
    let sigs: serde_json::Value = serde_json::from_slice(&fs::read(sigs).unwrap()).unwrap();
    sigs["num_records"].clone()
}

fn assert_refused(out: &Output) {
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !out.status.success(),
        "build over an existing index without --overwrite should fail:\n{stderr}"
    );
    assert!(
        stderr.contains("Cannot write over existing index without the --overwrite flag"),
        "unexpected error:\n{stderr}"
    );
}

#[test]
fn overwrite_without_decoys_removes_stale_poison_table() {
    let s = Scratch::new("poison");
    s.build_ok(&["-s", "refA.fa", "--decoy-paths", "decoys.fa", "-o", "idx/i"]);
    assert!(s.path("idx/i.poison").exists());
    assert!(s.path("idx/i.poison.json").exists());

    s.build_ok(&["-s", "refA.fa", "--overwrite", "-o", "idx/i"]);
    assert!(!s.path("idx/i.poison").exists());
    assert!(!s.path("idx/i.poison.json").exists());
}

#[test]
fn overwrite_with_sshash_dict_removes_stale_tiny_artifacts() {
    let s = Scratch::new("tiny");
    s.build_ok(&["-s", "refA.fa", "--dict", "tiny", "-o", "idx/i"]);
    assert!(s.path("idx/i.tdct").exists());
    assert!(s.path("idx/i.tct").exists());

    s.build_ok(&[
        "-s",
        "refB.fa",
        "--dict",
        "sshash",
        "--overwrite",
        "-o",
        "idx/i",
    ]);
    assert!(!s.path("idx/i.tdct").exists());
    assert!(!s.path("idx/i.tct").exists());
    // The signature pass has just rewritten `.sigs.json` for refB; --overwrite
    // must keep it, not remove it with the rest of the previous build.
    assert_eq!(num_records(&s.path("idx/i.sigs.json")), 40);
}

/// After `--overwrite`, the prefix holds exactly what a build of the same
/// command into an empty directory produces: nothing from the previous build
/// survives, whichever optional files it had.
#[test]
fn overwrite_leaves_only_the_new_build() {
    let s = Scratch::new("all");
    let everything = |prefix: &str| {
        s.build_ok(&[
            "-s",
            "refA.fa",
            "--decoy-paths",
            "decoys.fa",
            "--dict",
            "tiny",
            "--keep-intermediate-dbg",
            "-o",
            prefix,
        ]);
    };
    let minimal = |prefix: &str, overwrite: bool| {
        let mut args = vec!["-l", "listB.txt", "--dict", "sshash", "--no-ec-table"];
        if overwrite {
            args.push("--overwrite");
        }
        args.extend(["-o", prefix]);
        s.build_ok(&args);
    };

    // An undotted prefix, and a dotted one (`dot/gen.v2.ssi`, ...).
    for (prefix, fresh) in [("idx/i", "fresh/i"), ("dot/gen.v2", "dotfresh/gen.v2")] {
        let dir = prefix.split('/').next().unwrap();
        let fresh_dir = fresh.split('/').next().unwrap();
        everything(prefix);
        minimal(prefix, true);
        minimal(fresh, false);
        assert_eq!(s.names(dir), s.names(fresh_dir), "files under {dir}/");
        let refinfo = |p: &str| fs::read(format!("{}.refinfo", s.path(p).display())).unwrap();
        assert_eq!(refinfo(prefix), refinfo(fresh));
    }
    assert!(!s.path("dot/gen.v2.tdct").exists());
    assert!(!s.path("dot/gen.v2.ectab").exists());
    assert!(!s.path("dot/gen.v2.poison").exists());
    assert!(!s.path("idx/i.sigs.json").exists());
}

/// A refused rebuild must fail before it writes anything: previously the
/// signature pass had already replaced `.sigs.json` (and written the clipped
/// `.fa.gz`) by the time the overwrite check ran. It is refused even once
/// `_cfish.json` has been deleted, as long as other files of the index are
/// there.
#[test]
fn refused_rebuild_leaves_previous_build_untouched() {
    let s = Scratch::new("refused");
    s.build_ok(&["-s", "refA.fa", "-o", "idx/i"]);
    let rebuild = || {
        Command::new(env!("CARGO_BIN_EXE_piscem"))
            .current_dir(&s.0)
            .env("NO_COLOR", "1")
            .args(["build", "-s", "refB.fa", "--polya-clip-length", "10"])
            .args(["-o", "idx/i", "-t", "1", "-w", "newwork"])
            .output()
            .unwrap()
    };

    for delete_cfish_json in [false, true] {
        if delete_cfish_json {
            fs::remove_file(s.path("idx/i_cfish.json")).unwrap();
        }
        let before = s.files("idx");
        assert_refused(&rebuild());
        assert_eq!(
            s.files("idx"),
            before,
            "_cfish.json deleted: {delete_cfish_json}"
        );
        assert!(!s.path("newwork").exists());
        assert_eq!(num_records(&s.path("idx/i.sigs.json")), 20);
    }
}

/// The guard used to fire only when the cDBG json existed without its
/// `.cf_seg`/`.cf_seq`, so a `--keep-intermediate-dbg` index was silently
/// replaced by a rebuild without `--overwrite`.
#[test]
fn refused_rebuild_after_keep_intermediate_dbg() {
    let s = Scratch::new("keep");
    s.build_ok(&["-s", "refA.fa", "--keep-intermediate-dbg", "-o", "idx/i"]);
    let before = s.files("idx");

    assert_refused(&s.build(&["-s", "refB.fa", "-o", "idx/i"]));
    assert_eq!(s.files("idx"), before);
}

/// `--overwrite` removes nothing until cf1-rs has read every reference input,
/// so a bad one does not cost the previous index. A mistyped path, or a list
/// naming a file that does not exist (alone; after one that exists, in the
/// second of two lists; before one, in the first), is caught up front; a list
/// naming a file that is not FASTA or a directory without FASTA by cf1-rs; an
/// `-s` file that is not FASTA by the signature pass.
#[test]
fn overwrite_with_bad_input_keeps_previous_build() {
    let s = Scratch::new("badinput");
    s.build_ok(&["-s", "refA.fa", "--decoy-paths", "decoys.fa", "-o", "idx/i"]);
    let before = s.files("idx");
    assert!(before.contains_key("i.poison"));

    fs::write(
        s.path("listMissing.txt"),
        format!("{}\n", s.path("missing.fa").display()),
    )
    .unwrap();
    fs::write(
        s.path("listBMissing.txt"),
        format!(
            "{}\n{}\n",
            s.path("refB.fa").display(),
            s.path("missing.fa").display()
        ),
    )
    .unwrap();
    fs::write(
        s.path("listMissingB.txt"),
        format!(
            "{}\n{}\n",
            s.path("missing.fa").display(),
            s.path("refB.fa").display()
        ),
    )
    .unwrap();
    fs::create_dir(s.path("nofasta")).unwrap();
    fs::write(s.path("nofasta/README"), "no FASTA here\n").unwrap();
    fs::write(s.path("notfa.txt"), "not a FASTA file\n").unwrap();
    fs::write(
        s.path("listNotFa.txt"),
        format!("{}\n", s.path("notfa.txt").display()),
    )
    .unwrap();

    for (input, name) in [
        (["-l", "nolist.txt"], "nolist.txt"),
        (
            ["-l", "listB.txt,listBMissing.txt"],
            "missing.fa (listed in listBMissing.txt)",
        ),
        (
            ["-l", "listMissingB.txt,listB.txt"],
            "missing.fa (listed in listMissingB.txt)",
        ),
        (
            ["-l", "listMissing.txt"],
            "missing.fa (listed in listMissing.txt)",
        ),
        (["-l", "listNotFa.txt"], "failed to open FASTA file"),
        (["-d", "nofasta"], "No input files found"),
        (["-s", "notfa.txt"], "Expected '@' or '>'"),
    ] {
        let mut args = input.to_vec();
        args.extend(["--overwrite", "-o", "idx/i"]);
        let out = s.build(&args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "{input:?} should fail:\n{stderr}");
        assert!(
            stderr.contains(name),
            "{input:?}: unexpected error:\n{stderr}"
        );
        assert_eq!(
            s.files("idx"),
            before,
            "{input:?} changed the previous build"
        );
    }
}

/// The build writes a regular file at each of its paths, so anything else
/// there is refused before anything is written or removed: a directory would
/// fail part-way through `--overwrite`, and a symbolic link would be written
/// through.
#[test]
fn overwrite_refuses_a_non_file_at_an_index_path() {
    let s = Scratch::new("nonfile");
    s.build_ok(&["-s", "refA.fa", "-o", "idx/i"]);
    fs::remove_file(s.path("idx/i.ctab")).unwrap();
    fs::create_dir(s.path("idx/i.ctab")).unwrap();
    let before = s.names("idx");

    let out = s.build(&["-s", "refB.fa", "--overwrite", "-o", "idx/i"]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(stderr.contains("idx/i.ctab is a directory"), "{stderr}");
    assert_eq!(s.names("idx"), before);
    assert_eq!(num_records(&s.path("idx/i.sigs.json")), 20);

    #[cfg(unix)]
    {
        fs::create_dir(s.path("link")).unwrap();
        std::os::unix::fs::symlink("../elsewhere.ssi", s.path("link/i.ssi")).unwrap();
        let out = s.build(&["-s", "refA.fa", "-o", "link/i"]);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success());
        assert!(stderr.contains("link/i.ssi is a symbolic link"), "{stderr}");
        assert_eq!(s.names("link"), ["i.ssi"]);
        assert!(!s.path("elsewhere.ssi").exists());
    }
}

/// A build at `leg/gen` writes over the core files of an index that
/// piscem-rs <= 0.11.0 built at the dotted prefix `leg/gen.v1` (`leg/gen.ssi`,
/// ...), so it names that index in a warning, whether it is refused or, with
/// `--overwrite`, goes ahead. So it does for one known only by the poison
/// table `build-poison` wrote at its prefix. Not named: an index in the
/// current layout, what a failed build left (`.sigs.json`, `_cfish.json`),
/// and a `_ver.json` that records a later piscem-rs, even though each of them
/// would also load the core files at `leg/gen`.
#[test]
fn build_over_the_core_of_a_legacy_dotted_index_warns() {
    let s = Scratch::new("legacy");
    let ver = |piscem_rs: &str| format!("{{\"piscem-rs\": \"{piscem_rs}\"}}");
    s.build_ok(&["-s", "refA.fa", "-o", "leg/gen.v1"]);
    // Move the core files to where piscem-rs <= 0.11.0 wrote them, and record
    // that version.
    for ext in ["ssi", "ssi.mphf", "ctab", "refinfo", "ectab", "tdct", "tct"] {
        let from = s.path(&format!("leg/gen.v1.{ext}"));
        if from.exists() {
            fs::rename(from, s.path(&format!("leg/gen.{ext}"))).unwrap();
        }
    }
    fs::write(s.path("leg/gen.v1_ver.json"), ver("0.11.0")).unwrap();
    s.build_ok(&["-s", "refA.fa", "-o", "leg/gen.v3"]);
    for f in [".sigs.json", "_cfish.json"] {
        fs::copy(
            s.path(&format!("leg/gen.v1{f}")),
            s.path(&format!("leg/gen.v4{f}")),
        )
        .unwrap();
    }
    fs::write(s.path("leg/gen.v5_ver.json"), ver("0.11.1")).unwrap();
    fs::write(s.path("leg/gen.v6.poison"), "").unwrap();

    for overwrite in [false, true] {
        let mut args = vec!["-s", "refB.fa", "-o", "leg/gen"];
        if overwrite {
            args.push("--overwrite");
        }
        let out = s.build(&args);
        let stderr = String::from_utf8_lossy(&out.stderr);
        if overwrite {
            assert!(out.status.success(), "{stderr}");
        } else {
            assert_refused(&out);
        }
        for sib in ["leg/gen.v1", "leg/gen.v6"] {
            let warning =
                format!("are also the core files of the index {sib}, which piscem-rs <= 0.11.0");
            assert!(stderr.contains(&warning), "{args:?}: {sib}:\n{stderr}");
        }
        for sib in ["leg/gen.v3", "leg/gen.v4", "leg/gen.v5"] {
            assert!(!stderr.contains(sib), "{args:?}: {sib}:\n{stderr}");
        }
    }
}

/// An output or work directory that cannot be created, because a file or a
/// dangling symbolic link is at its path or at one of its parents', fails
/// before either directory is created.
#[test]
fn a_file_in_the_way_of_the_output_or_work_directory_writes_nothing() {
    let s = Scratch::new("dirisfile");
    fs::write(s.path("afile"), "").unwrap();
    let not_dir = "exists but is not a directory";
    // (-o, -w, error, the other directory, which must not be created)
    #[cfg_attr(not(unix), allow(unused_mut))]
    let mut cases = vec![
        (
            "out/i",
            "afile",
            format!("the work directory afile {not_dir}"),
            "out",
        ),
        (
            "out/i",
            "afile/sub",
            format!("cannot create the work directory afile/sub: afile {not_dir}"),
            "out",
        ),
        (
            "afile/i",
            "newwork",
            format!("the output directory afile {not_dir}"),
            "newwork",
        ),
        (
            "afile/sub/i",
            "newwork",
            format!("cannot create the output directory afile/sub: afile {not_dir}"),
            "newwork",
        ),
    ];
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink("nowhere", s.path("dangling")).unwrap();
        let not_dir = "is a symbolic link, but not to a directory";
        cases.push((
            "out/i",
            "dangling",
            format!("the work directory dangling {not_dir}"),
            "out",
        ));
        cases.push((
            "dangling/i",
            "newwork",
            format!("the output directory dangling {not_dir}"),
            "newwork",
        ));
    }
    for (prefix, work, error, other) in cases {
        let out = s.build_in(&["-s", "refA.fa", "-o", prefix], work);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "-o {prefix} -w {work} should fail");
        assert!(
            stderr.contains(&error),
            "-o {prefix} -w {work}: unexpected error:\n{stderr}"
        );
        assert!(
            !s.path(other).exists(),
            "-o {prefix} -w {work} created {other}"
        );
    }
}

/// The work directory is created before the signature pass. A build that fails
/// there (on an `-s` file that is not FASTA) removes the directories it created
/// for it again, parents included, so it leaves no new, empty work directory
/// behind. A directory that existed is kept, empty or not, also when `-w`
/// reaches it through a parent the build created (`missing/../empty`) or
/// (unix) is a symbolic link to it; so is the output directory the build
/// created before it, and so is the work directory once cf1-rs has it,
/// whether the build then succeeds or fails in cf1-rs. A file at the work
/// directory's path that only a created parent reveals (`missing/../afile`)
/// fails before the signature pass and leaves no parent behind.
#[test]
fn a_build_that_fails_before_cf1_rs_removes_the_work_directory_it_created() {
    let s = Scratch::new("newwork");
    fs::write(s.path("notfa.txt"), "not a FASTA file\n").unwrap();
    fs::create_dir(s.path("nofasta")).unwrap();
    fs::create_dir(s.path("empty")).unwrap();
    fs::create_dir(s.path("full")).unwrap();
    fs::write(s.path("full/keep.txt"), "kept\n").unwrap();
    let fails = |input: &[&str], work: &str, error: &str| {
        let mut args = input.to_vec();
        args.extend(["-o", "idx/i"]);
        let out = s.build_in(&args, work);
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(!out.status.success(), "{input:?} -w {work} should fail");
        assert!(
            stderr.contains(error),
            "{input:?} -w {work}: unexpected error:\n{stderr}"
        );
    };
    let (not_fasta, parse_error) = (["-s", "notfa.txt"], "Expected '@' or '>'");

    // (-w, the outermost directory it creates)
    for (work, created) in [
        ("newwork", "newwork"),
        ("new/deeper/work", "new"),
        ("empty/new/work", "empty/new"),
    ] {
        fails(&not_fasta, work, parse_error);
        assert!(!s.path(created).exists(), "-w {work} left {created} behind");
    }
    for work in ["empty", "full"] {
        fails(&not_fasta, work, parse_error);
        assert!(s.path(work).is_dir(), "-w {work} was removed");
    }
    // `missing/../empty` names `empty` only once `missing` has been created:
    // `missing` is removed again, `empty` (or `full`) is kept.
    for existing in ["empty", "full"] {
        let work = format!("missing/../{existing}");
        fails(&not_fasta, &work, parse_error);
        assert!(s.path(existing).is_dir(), "-w {work} removed {existing}");
        assert!(!s.path("missing").exists(), "-w {work} left missing behind");
    }
    assert_eq!(s.names("empty"), Vec::<String>::new());
    assert_eq!(s.names("full"), ["keep.txt"]);
    // A file reached through a parent the build creates (`missing/../afile`)
    // is found only when the work directory is created: that fails, before
    // the signature pass, and removes `missing` again.
    fs::write(s.path("afile"), "").unwrap();
    fails(
        &["-s", "refA.fa"],
        "missing/../afile",
        "failed to create working directory missing/../afile",
    );
    assert!(
        !s.path("missing").exists(),
        "-w missing/../afile left missing behind"
    );
    assert!(!s.path("idx/i.sigs.json").exists());
    // The output directory is created first and kept, as on main; only the
    // work directory created in it is removed.
    let out = s.build_in(&["-s", "notfa.txt", "-o", "newout/i"], "newout/work");
    assert!(!out.status.success());
    assert_eq!(s.names("newout"), Vec::<String>::new());
    // A symbolic link to a directory is used as that directory: kept, as is
    // what it holds; a directory created through the link is removed.
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink("full", s.path("lnkfull")).unwrap();
        for work in ["lnkfull", "lnkfull/new"] {
            fails(&not_fasta, work, parse_error);
        }
        assert!(s.path("lnkfull").is_symlink());
        assert_eq!(s.names("full"), ["keep.txt"]);
    }

    let out = s.build_in(&["-s", "refA.fa", "-o", "ok/i"], "okwork");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    fails(&["-d", "nofasta"], "cfwork", "No input files found");
    for work in ["okwork", "cfwork"] {
        assert!(s.path(work).is_dir(), "-w {work} was removed");
    }
}

/// Each `-l` list is read once, so one that can be read only once, a pipe
/// such as `-l /dev/stdin` or `-l <(ls *.fa)`, works too, alone or beside
/// another list, and builds exactly the index the same list in a file does.
/// It is read as cf1-rs reads a list: one path per line, trimmed, blank lines
/// skipped. And every list is used, in order, each path once: the lists build
/// the index that the FASTA files they name build when given to `-s`.
#[cfg(unix)]
#[test]
fn a_list_read_from_a_pipe_builds_the_same_index() {
    let s = Scratch::new("pipe");
    fs::write(
        s.path("listA.txt"),
        format!("{}\n", s.path("refA.fa").display()),
    )
    .unwrap();
    // listB.txt, with blank lines and blanks around the path.
    let piped = format!("\n  {}\t\r\n\n", s.path("refB.fa").display());
    for other in [vec![], vec!["-l", "listA.txt"]] {
        let mut args = vec!["-l", "listB.txt"];
        args.extend(&other);
        args.extend(["-o", "idx/i"]);
        s.build_ok(&args);
        fs::rename(s.path("idx"), s.path("plain")).unwrap();
        // The same index as the FASTA files the lists name, given to `-s`.
        let seqs = if other.is_empty() {
            "refB.fa"
        } else {
            "refB.fa,refA.fa"
        };
        s.build_ok(&["-s", seqs, "-o", "direct/i"]);
        for ext in ["refinfo", "ctab"] {
            assert!(
                fs::read(s.path(&format!("plain/i.{ext}"))).unwrap()
                    == fs::read(s.path(&format!("direct/i.{ext}"))).unwrap(),
                "{args:?}: .{ext} differs from -s {seqs}"
            );
        }
        fs::remove_dir_all(s.path("direct")).unwrap();

        args[1] = "/dev/stdin";
        let out = s.build_with_stdin(&args, &piped);
        assert!(
            out.status.success(),
            "piscem build {:?} failed:\n{}",
            args,
            String::from_utf8_lossy(&out.stderr)
        );
        let (plain, from_pipe) = (s.files("plain"), s.files("idx"));
        assert_eq!(
            plain.keys().collect::<Vec<_>>(),
            from_pipe.keys().collect::<Vec<_>>(),
            "{args:?}"
        );
        for (name, bytes) in &plain {
            assert!(from_pipe[name] == *bytes, "{args:?}: {name} differs");
        }
        fs::remove_dir_all(s.path("plain")).unwrap();
        fs::remove_dir_all(s.path("idx")).unwrap();
    }
}

/// `--overwrite` at a dotted prefix removes only that index's files: an index
/// at the stem (`sib/gen`) shares the names piscem-rs <= 0.11.0 gave the core
/// files of `sib/gen.v2`, and must survive, as must a sibling `sib/gen.v1`.
#[test]
fn overwrite_at_a_dotted_prefix_keeps_sibling_indexes() {
    let s = Scratch::new("sibling");
    s.build_ok(&["-s", "refA.fa", "-o", "sib/gen"]);
    s.build_ok(&["-s", "refA.fa", "-o", "sib/gen.v1"]);
    let before = s.files("sib");
    s.build_ok(&["-s", "refB.fa", "-o", "sib/gen.v2"]);
    s.build_ok(&[
        "-l",
        "listB.txt",
        "--dict",
        "sshash",
        "-o",
        "sib/gen.v2",
        "--overwrite",
    ]);
    let after = s.files("sib");
    for (name, bytes) in &before {
        assert_eq!(after.get(name), Some(bytes), "{name} changed");
    }
}
