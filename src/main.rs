use mimalloc::MiMalloc;
use std::io;
use std::path::{Path, PathBuf};

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use anyhow::{Context, Result, bail};
use cf1_rs::{CfInput, cf_build as cf1_build};
use clap::{Parser, Subcommand};
use serde_json::json;
use tracing::{Level, info, warn};

mod piscem_commands;
use piscem_commands::*;

use piscem_rs::cli::build as rs_build;
use piscem_rs::cli::build::BuildArgs;
use piscem_rs::cli::map_bulk;
use piscem_rs::cli::map_bulk::MapBulkArgs;
use piscem_rs::cli::map_scatac;
use piscem_rs::cli::map_scatac::MapScatacArgs;
use piscem_rs::cli::map_scrna;
use piscem_rs::cli::map_scrna::MapScrnaArgs;
use piscem_rs::cli::poison as rs_poison;
use piscem_rs::cli::poison::BuildPoisonArgs;
use piscem_rs::index::paths::{
    BUILD_ARTIFACT_EXTENSIONS, CTAB_EXT, DICT_EXT, IndexLayout, IndexPaths,
    POISON_ARTIFACT_EXTENSIONS, POISON_EXT, SIGS_EXT, index_artifact_path, validate_index_prefix,
};

/// Indexing and mapping to compacted colored de Bruijn graphs
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
struct Cli {
    /// be quiet: print only warnings and errors.
    #[arg(short, long)]
    quiet: bool,
    #[command(subcommand)]
    command: Commands,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// Index a reference sequence
    #[command(arg_required_else_help = true)]
    Build(BuildOpts),

    /// map reads for single-cell processing
    #[command(arg_required_else_help = true)]
    MapSC(MapSCOpts),

    /// map reads for bulk processing
    #[command(arg_required_else_help = true)]
    MapBulk(MapBulkOpts),

    /// map reads for scAtac processing
    #[command(arg_required_else_help = true)]
    MapSCAtac(MapSCAtacOpts),
}

/// Whether a file found in one of several `--ref-dirs` should be indexed.
///
/// cf1-rs picks the files of a single directory itself, by extension; several
/// directories are listed here, and cf1-rs panics on any file it is handed
/// that needletail cannot open as FASTA/FASTQ (a `.fai`, a `README.txt`, a
/// `genes.gtf.gz`, an empty file). So keep exactly the files needletail can
/// open, whatever their name: a build that worked before indexes the same
/// files. A file that cannot be read is an error, not a file to skip.
fn is_reference_file(path: &Path) -> Result<bool> {
    use needletail::errors::ParseErrorKind::{EmptyFile, UnknownFormat};
    use std::io::Read;
    let context = || format!("reading reference file {}", path.display());
    match needletail::parse_fastx_file(path) {
        Ok(_) => Ok(true),
        Err(e) if matches!(e.kind, UnknownFormat) => Ok(false),
        // needletail reports any error on its first 2-byte read as EmptyFile;
        // read them again so that an I/O error is not taken for a short file.
        Err(e) if matches!(e.kind, EmptyFile) => {
            let mut head = [0u8; 2];
            match std::fs::File::open(path).and_then(|mut f| f.read_exact(&mut head)) {
                Err(err) if err.kind() != io::ErrorKind::UnexpectedEof => {
                    Err(err).with_context(context)
                }
                _ => Ok(false),
            }
        }
        Err(e) => Err(e).with_context(context),
    }
}

/// What identifies the file at `path` (following symlinks), or `None` if there
/// is none: two paths get the same value exactly when they name the same file.
/// That is the (device, inode) pair on Unix, so a hard link is caught as well
/// as a symlink; elsewhere it is the canonical path.
#[cfg(unix)]
fn file_id(path: &Path) -> Option<(u64, u64)> {
    use std::os::unix::fs::MetadataExt;
    std::fs::metadata(path).ok().map(|m| (m.dev(), m.ino()))
}

#[cfg(not(unix))]
fn file_id(path: &Path) -> Option<PathBuf> {
    path.canonicalize().ok()
}

/// Validate a requested thread count against the parallelism actually available.
///
/// Zero is an error, but over-requesting is only a warning: the effective
/// budget is capped at `ncpus` (for the mapping commands by piscem-rs, which
/// records both the requested and effective values in `map_info.json`; for
/// `build` by the wrapper), so refusing outright would reject runs that work
/// fine.
fn check_threads(threads: usize, ncpus: usize) -> Result<()> {
    if threads == 0 {
        bail!(
            "the number of provided threads ({}) must be greater than 0.",
            threads
        );
    }
    if threads > ncpus {
        warn!(
            "the number of provided threads ({}) exceeds the available parallelism ({}); \
             the effective budget will be {}.",
            threads, ncpus, ncpus
        );
    }
    Ok(())
}

/// `prefix` with `suffix` appended to its file name, byte for byte, as in the
/// wrapper's own `<prefix>_cfish` and `<prefix>_ver.json`. Unlike
/// [`index_artifact_path`], this inserts no `.`.
fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name = prefix.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

/// `<prefix>_cfish`, the prefix cf1-rs is given. It appends [`CF1_SUFFIXES`].
fn cf_prefix(output: &Path) -> PathBuf {
    with_suffix(output, "_cfish")
}

/// What cf1-rs appends to [`cf_prefix`]: the cDBG's json, segment and tiling
/// files.
const CF1_SUFFIXES: [&str; 3] = [".json", ".cf_seg", ".cf_seq"];

/// `<prefix>_ver.json`, the versions of the tools that built the index.
fn ver_file(output: &Path) -> PathBuf {
    with_suffix(output, "_ver.json")
}

/// Every file `piscem build` can write for the output prefix `output`, named
/// exactly as its writer names it: the files of piscem-rs's index and poison
/// table and the reference signatures, `<prefix>.<ext>`
/// ([`index_artifact_path`]), and the wrapper's `<prefix>_cfish.*`
/// ([`cf_prefix`]) and `<prefix>_ver.json` ([`ver_file`]).
///
/// `<prefix>.fa.gz` (the `--polya-clip-length` output) is left out on purpose:
/// it is a natural name for the user's own input FASTA. So are the names
/// piscem-rs <= 0.11.0 gave the core files of a dotted prefix (`idx/gen.ssi`
/// for `-o idx/gen.v2`, its last dotted part replaced): they may belong to the
/// index at the stem (`idx/gen`) or to another old index sharing that stem
/// (`idx/gen.v1`), and piscem-rs never loads them for an index that has the
/// current `<prefix>.<ext>` files.
fn build_artifacts(output: &Path) -> Vec<PathBuf> {
    let cf_out = cf_prefix(output);
    let mut files: Vec<PathBuf> = CF1_SUFFIXES
        .into_iter()
        .map(|ext| with_suffix(&cf_out, ext))
        .collect();
    files.push(ver_file(output));
    files.extend(
        BUILD_ARTIFACT_EXTENSIONS
            .iter()
            .chain(POISON_ARTIFACT_EXTENSIONS)
            .chain(&[SIGS_EXT])
            .map(|ext| index_artifact_path(output, ext)),
    );
    files
}

/// The first of `dir` and its ancestors that exists (a dangling symbolic link
/// counts) but is not a directory (following links), so that creating `dir`
/// must fail; `None` if nothing in its path is known to be in the way.
fn non_directory_in_path(dir: &Path) -> Option<&Path> {
    dir.ancestors()
        .filter(|a| !a.as_os_str().is_empty())
        .find(|a| std::fs::symlink_metadata(a).is_ok())
        .filter(|a| !a.is_dir())
}

/// Directories a build has created, removed again, leaf first, if this is
/// dropped before it is cleared. `remove_dir` removes only an empty
/// directory, so nothing put in one is lost.
struct NewDirs(Vec<PathBuf>);

impl NewDirs {
    /// Create `dir` and each missing parent, as `create_dir_all(dir)` does,
    /// recording only the directories this call itself created. A directory
    /// that existed is never recorded, even under a name that does not
    /// resolve before a parent is created (`missing/../existing`).
    fn create(dir: &Path) -> io::Result<Self> {
        let mut created = NewDirs(Vec::new());
        let missing: Vec<&Path> = dir
            .ancestors()
            .take_while(|a| !a.as_os_str().is_empty() && !a.is_dir())
            .collect();
        for d in missing.into_iter().rev() {
            match std::fs::create_dir(d) {
                Ok(()) => created.0.insert(0, d.to_path_buf()),
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists && d.is_dir() => {}
                // Dropping `created` removes what this call has created.
                Err(e) => return Err(e),
            }
        }
        Ok(created)
    }
}

impl Drop for NewDirs {
    fn drop(&mut self) {
        for dir in &self.0 {
            let _ = std::fs::remove_dir(dir);
        }
    }
}

/// The piscem-rs version (major, minor, patch) that `<prefix>_ver.json`
/// records, if it exists and records one.
fn recorded_piscem_rs_version(prefix: &Path) -> Option<(u64, u64, u64)> {
    let ver: serde_json::Value =
        serde_json::from_slice(&std::fs::read(ver_file(prefix)).ok()?).ok()?;
    let mut parts = ver
        .get("piscem-rs")?
        .as_str()?
        .split(['.', '-', '+'])
        .map(|p| p.parse().ok());
    Some((parts.next()??, parts.next()??, parts.next()??))
}

/// Indexes beside `output` that piscem-rs <= 0.11.0 built at a dotted prefix
/// `<output>.<x>`, and whose core files are therefore `<output>.ssi`,
/// `<output>.ctab`, ...: the files a build at `output` writes.
///
/// Candidates are the prefixes `<output>.<x>` with a file that only a
/// finished build writes at its own prefix: the wrapper's `_ver.json` (last)
/// or `build-poison`'s `.poison` (after the core files). A failed build leaves
/// at most `.sigs.json` and `_cfish.*`, so it is not one. A candidate is kept
/// if [`IndexPaths::resolve`] resolves it to the [`IndexLayout::Legacy`]
/// layout, i.e. to the core files at `output`, and either its `_ver.json`
/// records piscem-rs <= 0.11.0 (a later one wrote the current layout, so its
/// core files cannot be at `output`) or, without a `_ver.json`, it has a
/// `.poison`. An index that piscem-rs's own `build` made without a poison
/// table has neither, so it is never named.
fn legacy_siblings(output: &Path) -> Vec<PathBuf> {
    let Some(name) = output.file_name() else {
        return Vec::new();
    };
    let dir = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let poison_suffix = format!(".{POISON_EXT}");
    let mut found: Vec<PathBuf> = entries
        .filter_map(|e| {
            let file = e.ok()?.file_name();
            let file = file.to_str()?;
            let sib = ["_ver.json", &poison_suffix]
                .into_iter()
                .find_map(|suffix| file.strip_suffix(suffix))?;
            // `<output>.<x>` with no further dot: its pre-0.11.1 core files
            // are at its prefix with the last extension dropped, `output`.
            (Path::new(sib).with_extension("").as_os_str() == name && sib != name)
                .then(|| output.with_file_name(sib))
        })
        .filter(|sib| {
            IndexPaths::resolve(sib).layout() == IndexLayout::Legacy
                && if ver_file(sib).exists() {
                    recorded_piscem_rs_version(sib).is_some_and(|v| v <= (0, 11, 0))
                } else {
                    index_artifact_path(sib, POISON_EXT).is_file()
                }
        })
        .collect();
    found.sort();
    found.dedup();
    found
}

fn main() -> Result<(), anyhow::Error> {
    let cli_args = Cli::parse();

    let quiet = cli_args.quiet;
    if quiet {
        tracing_subscriber::fmt()
            .with_max_level(Level::WARN)
            .with_writer(io::stderr)
            .init();
    } else {
        tracing_subscriber::fmt()
            .with_max_level(Level::INFO)
            .with_writer(io::stderr)
            .init();
    }

    // `available_parallelism` honours cpuset/cgroup limits, which
    // `num_cpus::get()` does not; piscem-rs caps the effective budget the same
    // way, so agreeing here keeps the wrapper from accepting a number that gets
    // silently clamped one layer down.
    let detected_ncpus = std::thread::available_parallelism().map(|n| n.get()).ok();
    let ncpus = detected_ncpus.unwrap_or(1);

    match cli_args.command {
        Commands::Build(BuildOpts {
            ref_seqs,
            ref_lists,
            ref_dirs,
            klen,
            mlen,
            threads,
            output,
            keep_intermediate_dbg,
            work_dir,
            overwrite,
            no_ec_table,
            polya_clip_length,
            decoy_paths,
            seed,
            dict,
            mphf,
            mphf_offline,
            tmp_dir,
            ram_limit_gib,
        }) => {
            info!("starting piscem build");
            // First, before anything (signatures, cf1-rs files, the work
            // directory) is written: every file of the index is named
            // `<output>.<ext>` or `<output>_<suffix>`, so `output` must end in
            // a file name, not name a directory as `-o idx/` does.
            validate_index_prefix(&output)?;
            // cf1-rs panics on a prefix that is not UTF-8 (after the signature
            // pass has written its file), and `map-* -i` takes only UTF-8, so
            // such a prefix can never become a usable index.
            if output.to_str().is_none() {
                bail!(
                    "invalid index prefix {:?}: it must be valid UTF-8, because cf1-rs and \
                     `piscem map-* -i` only accept UTF-8 paths",
                    output
                );
            }
            check_threads(threads, ncpus)?;
            // Unlike the mapping paths, nothing below the wrapper caps the
            // build's thread count: cf1-rs, SSHash and the poison-table builder
            // each spawn as many threads as they are given. Cap it here, as
            // `check_threads` has just reported. If the parallelism can't be
            // detected, `ncpus` is only a fallback, so keep the request as
            // given, as piscem-rs's mapping budget does.
            let threads = detected_ncpus.map_or(threads, |n| threads.min(n));
            if mlen >= klen {
                bail!(
                    "minimizer length ({}) must be < k-mer length ({})",
                    mlen,
                    klen
                );
            }

            // if the decoy sequences are provided, ensure they are valid paths
            if let Some(ref decoys) = decoy_paths {
                for d in decoys {
                    match d.try_exists() {
                        Ok(true) => {}
                        Ok(false) => {
                            bail!(
                                "Path for decoy file {} seems not to point to a valid file",
                                d.display()
                            );
                        }
                        Err(e) => {
                            bail!(
                                "Error {} when checking the existence of decoy file {}",
                                e,
                                d.display()
                            );
                        }
                    }
                }
            }

            // Likewise for the reference inputs, so that a mistyped path fails
            // here, before anything is written (or removed by --overwrite).
            for (inputs, what) in [
                (&ref_seqs, "reference file"),
                (&ref_lists, "reference list file"),
                (&ref_dirs, "reference directory"),
            ] {
                for p in inputs.iter().flatten() {
                    if !Path::new(p)
                        .try_exists()
                        .with_context(|| format!("checking the existence of {} {}", what, p))?
                    {
                        bail!("Path for {} {} does not exist", what, p);
                    }
                }
            }
            // And the files each list names: one per line, trimmed, blank lines
            // skipped, as cf1-rs reads a list; it would panic on a missing file
            // without naming it. Each list is read only here, and cf1-rs is
            // given the paths: a list may be a pipe (`-l <(ls *.fa)`,
            // `-l /dev/stdin`), which a second read would find empty.
            let mut listed: Vec<PathBuf> = Vec::new();
            for list_path in ref_lists.iter().flatten() {
                let contents = std::fs::read_to_string(list_path)
                    .with_context(|| format!("reading list file {}", list_path))?;
                for p in contents.lines().map(str::trim).filter(|l| !l.is_empty()) {
                    listed.push(PathBuf::from(p));
                    if !Path::new(p).try_exists().with_context(|| {
                        format!(
                            "checking the existence of reference file {} (listed in {})",
                            p, list_path
                        )
                    })? {
                        bail!(
                            "Path for reference file {} (listed in {}) does not exist",
                            p,
                            list_path
                        );
                    }
                }
            }

            // Every file a previous build may have left at this prefix. Without
            // --overwrite, refuse to touch any of them; with it, they are all
            // removed once cf1-rs has succeeded (below), because piscem-rs
            // loads the optional ones (EC table, Tiny dictionary, poison
            // table) whenever they exist, and a stale one would silently
            // change, or crash, mapping against the new index. The build
            // writes a regular file at each of these paths, so anything else
            // there (a directory, a symbolic link, which would be written
            // through) is refused before anything is written or removed.
            let artifacts = build_artifacts(&output);
            let sigs_file = index_artifact_path(&output, SIGS_EXT);
            let mut existing: Vec<&PathBuf> = Vec::new();
            for f in &artifacts {
                match std::fs::symlink_metadata(f) {
                    Ok(md) if md.is_file() => existing.push(f),
                    Ok(md) => {
                        let what = if md.is_dir() {
                            "is a directory"
                        } else if md.is_symlink() {
                            "is a symbolic link"
                        } else {
                            "is not a regular file"
                        };
                        bail!(
                            "{} {}, but piscem build writes an index file there; \
                             remove it or choose another output prefix.",
                            f.display(),
                            what
                        );
                    }
                    // Absent. A parent that is not a directory is reported
                    // below.
                    Err(e)
                        if matches!(
                            e.kind(),
                            io::ErrorKind::NotFound | io::ErrorKind::NotADirectory
                        ) => {}
                    Err(e) => {
                        return Err(e).with_context(|| format!("failed to check {}", f.display()));
                    }
                }
            }
            for sib in legacy_siblings(&output) {
                warn!(
                    "{}, {}, ... are also the core files of the index {sib}, which \
                     piscem-rs <= 0.11.0 built in its pre-0.11.1 layout. Building {} \
                     replaces them, and {sib} would then load the new index's files \
                     instead of its own: rebuild {sib}, which gives it its own {sib}.* \
                     files, or remove it.",
                    index_artifact_path(&output, DICT_EXT).display(),
                    index_artifact_path(&output, CTAB_EXT).display(),
                    output.display(),
                    sib = sib.display(),
                );
            }
            if !overwrite && !existing.is_empty() {
                warn!(
                    "The prefix you have chosen for output already holds files from a previous index build:"
                );
                for f in &existing {
                    warn!("  {}", f.display());
                }
                warn!(
                    "Please either delete these files, choose another output prefix, or use the --overwrite flag."
                );
                bail!("Cannot write over existing index without the --overwrite flag.");
            }

            // An output or work directory that cannot be created, because
            // something other than a directory is at its path or at one of its
            // parents' (`-w afile`, `-w afile/sub`, a dangling symbolic link),
            // fails here, before either directory is created.
            let out_dir = output.parent().filter(|p| !p.as_os_str().is_empty());
            for (dir, what) in [(out_dir, "output"), (Some(work_dir.as_path()), "work")] {
                if let Some(dir) = dir
                    && let Some(blocker) = non_directory_in_path(dir)
                {
                    let problem = if blocker.is_symlink() {
                        "is a symbolic link, but not to a directory"
                    } else {
                        "exists but is not a directory"
                    };
                    if blocker == dir {
                        bail!("the {} directory {} {}", what, dir.display(), problem);
                    }
                    bail!(
                        "cannot create the {} directory {}: {} {}",
                        what,
                        dir.display(),
                        blocker.display(),
                        problem
                    );
                }
            }

            // Ensure the parent directory of the output prefix exists.
            if let Some(parent) = output.parent()
                && !parent.as_os_str().is_empty()
                && !parent.exists()
            {
                std::fs::create_dir_all(parent).with_context(|| {
                    format!("failed to create output directory {}", parent.display())
                })?;
                info!(
                    "created output directory {} (did not previously exist)",
                    parent.display()
                );
            }

            // Ensure work directory exists. A build that fails before cf1-rs
            // is given it (on an `-s` file that is not FASTA, say) removes the
            // directories created here again, so it leaves no new, empty work
            // directory behind.
            let mut new_work_dirs = if work_dir.is_dir() {
                info!(
                    "will use {} as the work directory for temporary files.",
                    work_dir.display()
                );
                NewDirs(Vec::new())
            } else {
                NewDirs::create(&work_dir).with_context(|| {
                    format!(
                        "failed to create working directory {} for index construction",
                        work_dir.display()
                    )
                })?
            };

            // Build CfInput from CLI args, using native variants where possible.
            // With --polya-clip-length, cf1-rs reads the clipped reference
            // instead; it is removed along with the cf1-rs intermediates.
            let mut clipped_ref = None;
            let cf_input = if let Some(ref seqs) = ref_seqs {
                if !seqs.is_empty() {
                    // prepare_fasta appends `.json` to this stem, so it writes
                    // `<output>.sigs.json`, the file piscem-rs reads.
                    let out_stem = index_artifact_path(&output, SIGS_EXT).with_extension("");
                    if polya_clip_length.is_some() {
                        // prepare_fasta truncates `<out>.fa.gz` before it reads
                        // the inputs (and it is removed after indexing), so it
                        // must be neither a reference nor a decoy file.
                        let clipped = out_stem.with_extension("fa.gz");
                        if let Some(clipped_id) = file_id(&clipped) {
                            let is_clipped = |p: &Path| file_id(p).as_ref() == Some(&clipped_id);
                            if let Some(s) = seqs.iter().find(|s| is_clipped(Path::new(s))) {
                                bail!(
                                    "--polya-clip-length writes the clipped reference to {}, \
                                     which is the input {}; choose another --output prefix.",
                                    clipped.display(),
                                    s
                                );
                            }
                            if let Some(d) = decoy_paths.iter().flatten().find(|d| is_clipped(d)) {
                                bail!(
                                    "--polya-clip-length writes the clipped reference to {}, \
                                     which is the decoy file {}; choose another --output prefix.",
                                    clipped.display(),
                                    d.display()
                                );
                            }
                            warn!(
                                "{} already exists; it will be overwritten with the clipped reference.",
                                clipped.display()
                            );
                        }
                        clipped_ref = Some(clipped);
                    }
                    let configs = prepare_fasta::RecordParseConfig {
                        input: seqs.clone(),
                        output_stem: out_stem,
                        polya_clip_length,
                    };
                    info!("Computing and recording reference signatures...");
                    // Without clipping these are the input files; with it, the
                    // one clipped file (the signatures are of the unclipped input).
                    let files = prepare_fasta::parse_records(configs).with_context(|| {
                        format!(
                            "failed to record reference signatures in {}",
                            sigs_file.display()
                        )
                    })?;
                    info!("done.");
                    CfInput::Files(files)
                } else {
                    bail!("--ref-seqs provided but empty");
                }
            } else if ref_lists.is_some() {
                // The paths every list names, in order, read above.
                CfInput::Files(listed)
            } else if let Some(ref dirs) = ref_dirs {
                if dirs.len() == 1 {
                    CfInput::Directory(PathBuf::from(&dirs[0]))
                } else {
                    // Multiple directories — resolve all to individual file paths.
                    let mut files = Vec::new();
                    for dir_path in dirs {
                        for entry in std::fs::read_dir(dir_path)
                            .with_context(|| format!("reading directory {}", dir_path))?
                        {
                            let path = entry
                                .with_context(|| format!("reading directory {}", dir_path))?
                                .path();
                            if path.is_file() {
                                if is_reference_file(&path)? {
                                    files.push(path);
                                } else {
                                    warn!(
                                        "skipping {}: empty, or not a FASTA/FASTQ file",
                                        path.display()
                                    );
                                }
                            }
                        }
                    }
                    CfInput::Files(files)
                }
            } else {
                bail!("Input (via --ref-seqs, --ref-lists, or --ref-dirs) must be provided.");
            };

            // cf1-rs appends its extensions to the prefix it is given.
            let cf_out = cf_prefix(&output);

            // Ensure the output parent directory exists.
            if let Some(parent_path) = cf_out.parent()
                && !parent_path.as_os_str().is_empty()
                && !parent_path.exists()
            {
                std::fs::create_dir_all(parent_path).with_context(|| {
                    format!(
                        "failed to create output directory {}",
                        parent_path.display()
                    )
                })?;
                info!(
                    "directory {} did not already exist; creating it.",
                    parent_path.display()
                );
            }

            info!("starting cDBG construction with cf1-rs");
            // cf1-rs has the work directory from here on; it is kept.
            new_work_dirs.0.clear();
            let cf_result = cf1_build()
                .input(cf_input)
                .output_prefix(cf_out.clone())
                .k(klen)
                .threads(threads)
                .work_dir(work_dir.clone())
                .call()?;

            // Only now, with every reference input read and the cDBG built,
            // clear the previous build, so a build that fails on its reference
            // inputs (a list naming a file that is not FASTA, a directory
            // holding none) keeps the previous index. Kept: the files cf1-rs
            // has just rewritten and, for -s, the signatures the pass above
            // has just rewritten. The decoys are read only by build-poison,
            // after piscem-rs has built the new index, so a bad decoy file
            // leaves that index without a poison table.
            if overwrite {
                let rewritten = [
                    &cf_result.json_file,
                    &cf_result.seg_file,
                    &cf_result.seq_file,
                ];
                for f in existing {
                    if rewritten.contains(&f) || (ref_seqs.is_some() && *f == sigs_file) {
                        continue;
                    }
                    info!(
                        "--overwrite: removing {} from the previous build",
                        f.display()
                    );
                    std::fs::remove_file(f).with_context(|| {
                        format!("failed to remove {} from the previous build", f.display())
                    })?;
                }
            }

            info!(
                "cDBG construction complete: {} unitigs, {} vertices",
                cf_result.unitig_count, cf_result.vertex_count
            );

            // Build piscem-rs index from cf1-rs output
            rs_build::run(BuildArgs {
                input: cf_out.clone(),
                output: output.clone(),
                klen,
                mlen,
                threads,
                build_ec_table: !no_ec_table,
                // Deprecated in piscem-rs 0.10 (the index is always
                // canonical); false avoids the accept-and-warn notice.
                canonical: false,
                seed,
                mphf,
                mphf_offline,
                single_mphf: false,
                dict,
                // SSHash's temporary files (external sort, offline MPHF
                // records) go inside the work directory unless placed
                // explicitly, not in `sshash_tmp` under the current directory.
                tmp_dir: Some(tmp_dir.unwrap_or_else(|| work_dir.join("sshash_tmp"))),
                ram_limit_gib,
            })?;

            // Build poison table if decoys were provided
            if let Some(decoy_pathbufs) = decoy_paths {
                rs_poison::run(BuildPoisonArgs {
                    index: output.clone(),
                    decoys: decoy_pathbufs,
                    threads,
                })?;
            }

            if !keep_intermediate_dbg {
                info!("removing intermediate cdBG files produced by cf1-rs.");
                for (path, label) in [
                    (&cf_result.seg_file, "segment"),
                    (&cf_result.seq_file, "tiling"),
                ] {
                    match std::fs::remove_file(path) {
                        Ok(_) => info!("removed {} file {}", label, path.display()),
                        Err(e) => warn!("cannot remove {}, error: {:?}", path.display(), e),
                    }
                }
                // Keep the json file — it's small and may contain useful info.
                if let Some(path) = &clipped_ref {
                    match std::fs::remove_file(path) {
                        Ok(_) => info!("removed polyA-clipped reference {}", path.display()),
                        Err(e) => warn!("cannot remove {}, error: {:?}", path.display(), e),
                    }
                }
            }

            let piscem_rs_ver = piscem_rs::VERSION;
            let cf1_rs_ver = env!("CF1_RS_VERSION");
            let piscem_ver = clap::crate_version!();

            let version_json = json!({
                "piscem-rs": piscem_rs_ver,
                "cf1-rs": cf1_rs_ver,
                "piscem": piscem_ver
            });

            let fname = ver_file(&output);
            let ver_out = std::fs::File::create(&fname)
                .with_context(|| format!("failed to create {}", fname.display()))?;
            serde_json::to_writer_pretty(ver_out, &version_json)
                .with_context(|| format!("failed to write {}", fname.display()))?;

            info!("piscem build finished.");
        }

        Commands::MapSC(sc_opts) => {
            check_threads(sc_opts.threads, ncpus)?;

            let args = MapScrnaArgs {
                index: PathBuf::from(&sc_opts.index),
                read1: sc_opts.read1.iter().map(PathBuf::from).collect(),
                read2: sc_opts.read2.iter().map(PathBuf::from).collect(),
                geometry: sc_opts.geometry,
                output: sc_opts.output.clone(),
                threads: sc_opts.threads,
                decoder: sc_opts.decode.decoder,
                thread_policy: sc_opts.decode.thread_policy,
                skipping_strategy: sc_opts.skipping_strategy,
                no_poison: sc_opts.no_poison,
                struct_constraints: sc_opts.struct_constraints,
                ignore_ambig_hits: sc_opts.ignore_ambig_hits,
                max_ec_card: sc_opts.max_ec_card,
                max_hit_occ: sc_opts.max_hit_occ as usize,
                max_hit_occ_recover: sc_opts.max_hit_occ_recover as usize,
                max_read_occ: sc_opts.max_read_occ as usize,
                with_position: sc_opts.with_position,
                quiet,
                dict: sc_opts.dict,
            };
            map_scrna::run(args)?;
        }

        Commands::MapSCAtac(scatac_opts) => {
            check_threads(scatac_opts.threads, ncpus)?;

            for (flag, name) in [
                (scatac_opts.sam_format, "--sam-format"),
                (scatac_opts.bed_format, "--bed-format"),
                (scatac_opts.use_chr, "--use-chr"),
                (scatac_opts.check_kmer_orphan, "--check-kmer-orphan"),
                (scatac_opts.struct_constraints, "--struct-constraints"),
            ] {
                if flag {
                    warn!("{} is not supported by piscem-rs and will be ignored", name);
                }
            }
            if scatac_opts.skipping_strategy.is_some() {
                warn!(
                    "--skipping-strategy has no effect for scATAC (every k-mer is used) and will be ignored"
                );
            }

            let barcode = scatac_opts.barcode.unwrap_or_default();
            let args = MapScatacArgs {
                index: PathBuf::from(&scatac_opts.index),
                reads: scatac_opts
                    .reads
                    .unwrap_or_default()
                    .iter()
                    .map(PathBuf::from)
                    .collect(),
                read1: scatac_opts
                    .read1
                    .unwrap_or_default()
                    .iter()
                    .map(PathBuf::from)
                    .collect(),
                read2: scatac_opts
                    .read2
                    .unwrap_or_default()
                    .iter()
                    .map(PathBuf::from)
                    .collect(),
                barcode: barcode.iter().map(PathBuf::from).collect(),
                output: scatac_opts.output.clone(),
                threads: scatac_opts.threads,
                decoder: scatac_opts.decode.decoder,
                thread_policy: scatac_opts.decode.thread_policy,
                bc_len: scatac_opts.bclen as usize,
                no_tn5_shift: scatac_opts.no_tn5_shift,
                no_poison: scatac_opts.no_poison,
                check_ambig_hits: !scatac_opts.ignore_ambig_hits,
                max_ec_card: scatac_opts.max_ec_card,
                max_hit_occ: scatac_opts.max_hit_occ as usize,
                max_hit_occ_recover: scatac_opts.max_hit_occ_recover as usize,
                max_read_occ: scatac_opts.max_read_occ as usize,
                end_cache_capacity: scatac_opts.end_cache_capacity,
                bin_size: scatac_opts.bin_size as u64,
                bin_overlap: scatac_opts.bin_overlap as u64,
                thr: scatac_opts.thr,
                min_overlap: 30,
                skipping_strategy: None,
                quiet,
                dict: scatac_opts.dict,
            };
            map_scatac::run(args)?;
        }

        Commands::MapBulk(bulk_opts) => {
            check_threads(bulk_opts.threads, ncpus)?;

            let args = MapBulkArgs {
                index: PathBuf::from(&bulk_opts.index),
                reads: bulk_opts
                    .reads
                    .unwrap_or_default()
                    .iter()
                    .map(PathBuf::from)
                    .collect(),
                read1: bulk_opts
                    .read1
                    .unwrap_or_default()
                    .iter()
                    .map(PathBuf::from)
                    .collect(),
                read2: bulk_opts
                    .read2
                    .unwrap_or_default()
                    .iter()
                    .map(PathBuf::from)
                    .collect(),
                output: bulk_opts.output.clone(),
                threads: bulk_opts.threads,
                decoder: bulk_opts.decode.decoder,
                thread_policy: bulk_opts.decode.thread_policy,
                skipping_strategy: bulk_opts.skipping_strategy,
                no_poison: bulk_opts.no_poison,
                struct_constraints: bulk_opts.struct_constraints,
                ignore_ambig_hits: bulk_opts.ignore_ambig_hits,
                max_ec_card: bulk_opts.max_ec_card,
                max_hit_occ: bulk_opts.max_hit_occ as usize,
                max_hit_occ_recover: bulk_opts.max_hit_occ_recover as usize,
                max_read_occ: bulk_opts.max_read_occ as usize,
                quiet,
                dict: bulk_opts.dict,
            };
            map_bulk::run(args)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::Cli;
    use clap::CommandFactory;

    /// Catches clap misconfigurations (conflicting short flags, bad defaults,
    /// ...) that would otherwise only surface when the affected subcommand is
    /// parsed. clap only runs these checks with debug assertions enabled, so
    /// this guards `cargo test`, not `cargo test --release`.
    #[test]
    fn cli_is_consistent() {
        Cli::command().debug_assert();
    }
}
