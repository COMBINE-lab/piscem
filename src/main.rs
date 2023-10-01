use std::ffi::CString;
use std::ffi::{OsStr, OsString};
use std::io;
use std::os::raw::{c_char, c_int};
use std::path::PathBuf;

use prepare_fasta;
use anyhow::{bail, Result};
use clap::{Parser, Subcommand};
use tracing::{error, info, warn, Level};

mod piscem_commands;
use piscem_commands::*;

#[link(name = "pesc_static", kind = "static")]
extern "C" {
    pub fn run_pesc_sc(args: c_int, argsv: *const *const c_char) -> c_int;
    pub fn run_pesc_bulk(args: c_int, argsv: *const *const c_char) -> c_int;
}

#[link(name = "build_static", kind = "static")]
extern "C" {
    pub fn run_build(args: c_int, argsv: *const *const c_char) -> c_int;
    pub fn run_build_poison_table(args: c_int, argsv: *const *const c_char) -> c_int;
}

#[link(name = "cfcore_static", kind = "static", modifiers = "+whole-archive")]
extern "C" {
    pub fn cf_build(args: c_int, argsv: *const *const c_char) -> c_int;
}

/// Indexing and mapping to compacted colored de Bruijn graphs
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
struct Cli {
    /// be quiet (no effect yet for cDBG building phase of indexing).
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
}

// from: https://stackoverflow.com/questions/74322541/how-to-append-to-pathbuf
fn append_to_path(p: impl Into<OsString>, s: impl AsRef<OsStr>) -> PathBuf {
    let mut p = p.into();
    p.push(s);
    p.into()
}

fn main() -> Result<(), anyhow::Error> {
    let cli_args = Cli::parse();
    //env_logger::Builder::from_env(Env::default().default_filter_or("warn")).init();

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

    let ncpus = num_cpus::get();

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
            decoy_paths,
        }) => {
            info!("starting piscem build");
            if threads == 0 {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    threads
                );
            }
            if threads > ncpus {
                bail!("the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    threads, ncpus);
            }
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

            let mut args: Vec<CString> = vec![];

            let cf_out = PathBuf::from(output.as_path().to_string_lossy().into_owned() + "_cfish");
            let cf_base_path = cf_out.as_path();
            let seg_file = append_to_path(cf_base_path, ".cf_seg");
            let seq_file = append_to_path(cf_base_path, ".cf_seq");
            let struct_file = append_to_path(cf_base_path, ".json");
            let mut build_ret;

            if overwrite {
                if struct_file.exists() {
                    std::fs::remove_file(struct_file.clone())?;
                }
                if seg_file.exists() {
                    std::fs::remove_file(seg_file.clone())?;
                }
                if seq_file.exists() {
                    std::fs::remove_file(seq_file.clone())?;
                }
            }

            if struct_file.exists() && (!seq_file.exists() || !seg_file.exists()) {
                warn!("The prefix you have chosen for output already corresponds to an existing cDBG structure file {:?}.", struct_file.display());
                warn!("However, the corresponding seq and seg files do not exist. Please either delete this structure file, choose another output prefix, or use the --overwrite flag.");
                bail!("Cannot write over existing index without the --overwrite flag.");
            }

            args.push(CString::new("cdbg_builder").unwrap());

            // We can treat the different input options independently
            // here because the argument parser should have enforced
            // their exclusivity.
            let mut has_input = false;

            if let Some(seqs) = ref_seqs {
                if !seqs.is_empty() {
                    let out_stem = PathBuf::from(output.as_path().to_string_lossy().into_owned() + ".sigs");
                    let configs = prepare_fasta::RecordParseConfig{
                            input: seqs.clone(),
                            output_stem: out_stem,
                            polya_clip_length: None
                        };
                    prepare_fasta::parse_records(configs)?;
                    args.push(CString::new("--seq").unwrap());
                    let reflist = seqs.join(",");
                    args.push(CString::new(reflist.as_str()).unwrap());
                    has_input = true;
                }
            }

            if let Some(lists) = ref_lists {
                if !lists.is_empty() {
                    args.push(CString::new("--list").unwrap());
                    let reflist = lists.join(",");
                    args.push(CString::new(reflist.as_str()).unwrap());
                    has_input = true;
                }
            }

            if let Some(dirs) = ref_dirs {
                if !dirs.is_empty() {
                    args.push(CString::new("--dir").unwrap());
                    let reflist = dirs.join(",");
                    args.push(CString::new(reflist.as_str()).unwrap());
                    has_input = true;
                }
            }

            assert!(
                has_input,
                "Input (via --ref-seqs, --ref-lists, or --ref-dirs) must be provided."
            );

            args.push(CString::new("-k").unwrap());
            args.push(CString::new(klen.to_string()).unwrap());
            args.push(CString::new("--track-short-seqs").unwrap());
            args.push(CString::new("--poly-N-stretch").unwrap());

            // check if the provided work directory exists.
            // If not, then try and create it.
            match work_dir.try_exists() {
                Ok(true) => {
                    info!(
                        "will use {} as the work directory for temporary files.",
                        work_dir.display()
                    );
                }
                Ok(false) => {
                    // try to create it
                    match std::fs::create_dir_all(&work_dir) {
                        Ok(_) => {}
                        Err(e) => {
                            error!("when attempting to create working directory {}, encountered error {:#?}", &work_dir.display(), e);
                            bail!(
                                "Failed to create working directory {} for index construction : {:#?}",
                                &work_dir.display(),
                                e
                            );
                        }
                    }
                }
                Err(e) => {
                    error!("when checking existence of working directory {:#?}", e);
                    bail!(
                        "Failed to create working directory for index construction : {:#?}",
                        e
                    );
                }
            }

            // check if the provided output path is more than just a prefix
            // if so, check if the specified directory exists and create it
            // if it doesn't.
            if let Some(parent_path) = cf_out.parent() {
                if !parent_path.exists() {
                    std::fs::create_dir_all(parent_path)?;
                    info!(
                        "directory {} did not already exist; creating it.",
                        parent_path.display()
                    );
                }
            }

            args.push(CString::new("-o").unwrap());
            args.push(CString::new(cf_out.as_path().to_string_lossy().into_owned()).unwrap());

            args.push(CString::new("-t").unwrap());
            args.push(CString::new(threads.to_string()).unwrap());
            // output format
            args.push(CString::new("-f").unwrap());
            args.push(CString::new("3").unwrap());
            // work directory
            args.push(CString::new("-w").unwrap());
            args.push(CString::new(work_dir.as_path().to_string_lossy().into_owned()).unwrap());

            info!("args = {:?}", args);
            {
                let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
                let args_len: c_int = args.len() as c_int;
                build_ret = unsafe { cf_build(args_len, arg_ptrs.as_ptr()) };
            }

            if build_ret != 0 {
                bail!(
                    "cDBG constructor returned exit code {}; failure.",
                    build_ret
                );
            }

            args.clear();
            args.push(CString::new("ref_index_builder").unwrap());

            args.push(CString::new("-i").unwrap());
            args.push(CString::new(cf_out.as_path().to_string_lossy().into_owned()).unwrap());
            args.push(CString::new("-k").unwrap());
            args.push(CString::new(klen.to_string()).unwrap());
            args.push(CString::new("-m").unwrap());
            args.push(CString::new(mlen.to_string()).unwrap()); // minimizer length

            args.push(CString::new("--canonical-parsing").unwrap());
            if !no_ec_table {
                args.push(CString::new("--build-ec-table").unwrap());
            }
            args.push(CString::new("-o").unwrap());
            args.push(CString::new(output.as_path().to_string_lossy().into_owned()).unwrap());

            args.push(CString::new("-d").unwrap());
            args.push(CString::new(work_dir.as_path().to_string_lossy().into_owned()).unwrap());

            args.push(CString::new("-t").unwrap());
            args.push(CString::new(threads.to_string()).unwrap());

            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            {
                println!("{:?}", args);
                let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
                let args_len: c_int = args.len() as c_int;
                build_ret = unsafe { run_build(args_len, arg_ptrs.as_ptr()) };
            }

            if build_ret != 0 {
                bail!("indexer returned exit code {}; failure.", build_ret);
            }

            // now, build the poison table if there are decoys
            if let Some(decoy_pathbufs) = decoy_paths {
                args.clear();
                args.push(CString::new("poison_table_builder").unwrap());

                // index is the one we just built
                args.push(CString::new("-i").unwrap());
                args.push(CString::new(output.as_path().to_string_lossy().into_owned()).unwrap());

                args.push(CString::new("-t").unwrap());
                args.push(CString::new(threads.to_string()).unwrap());

                if overwrite {
                    args.push(CString::new("--overwrite").unwrap());
                }

                let path_args = decoy_pathbufs
                    .into_iter()
                    .map(|x| x.to_string_lossy().into_owned())
                    .collect::<Vec<String>>()
                    .join(",");
                args.push(CString::new("-d").unwrap());
                args.push(CString::new(path_args).unwrap());

                if quiet {
                    args.push(CString::new("--quiet").unwrap());
                }

                {
                    println!("{:?}", args);
                    let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
                    let args_len: c_int = args.len() as c_int;
                    build_ret = unsafe { run_build_poison_table(args_len, arg_ptrs.as_ptr()) };
                }
                if build_ret != 0 {
                    bail!(
                        "building poison table returned exit code {}; failure.",
                        build_ret
                    );
                }
            }

            if !keep_intermediate_dbg {
                info!("removing intermediate cdBG files produced by cuttlefish.");

                match std::fs::remove_file(seg_file.clone()) {
                    Ok(_) => {
                        info!("removed segment file {}", seg_file.display());
                    }
                    Err(e) => {
                        warn!(
                            "cannot remove {}, encountered error {:?}!",
                            seg_file.display(),
                            e
                        );
                    }
                };

                match std::fs::remove_file(seq_file.clone()) {
                    Ok(_) => {
                        info!("removed tiling file {}", seq_file.display());
                    }
                    Err(e) => {
                        warn!(
                            "cannot remove {}, encountered error {:?}!",
                            seq_file.display(),
                            e
                        );
                    }
                };
                // for now, let the json file stick around. It's
                // generally very small and may contain useful information
                // about the references being indexed.
            }

            info!("piscem build finished");
        }

        Commands::MapSC(sc_opts) => {
            if sc_opts.threads == 0 {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    sc_opts.threads
                );
            }
            if sc_opts.threads > ncpus {
                bail!("the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    sc_opts.threads, ncpus);
            }

            let mut args = sc_opts.as_argv()?;
            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            info!("cmd: {:?}", args);
            let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
            let args_len: c_int = args.len() as c_int;

            let map_ret = unsafe { run_pesc_sc(args_len, arg_ptrs.as_ptr()) };
            if map_ret != 0 {
                bail!("mapper returned exit code {}; failure", map_ret);
            }
        }

        Commands::MapBulk(bulk_opts) => {
            if bulk_opts.threads == 0 {
                bail!(
                    "the number of provided threads ({}) must be greater than 0.",
                    bulk_opts.threads
                );
            }
            if bulk_opts.threads > ncpus {
                bail!("the number of provided threads ({}) should be <= the number of logical CPUs ({}).",
                    bulk_opts.threads, ncpus);
            }

            let mut args = bulk_opts.as_argv()?;

            if quiet {
                args.push(CString::new("--quiet").unwrap());
            }

            let arg_ptrs: Vec<*const c_char> = args.iter().map(|s| s.as_ptr()).collect();
            let args_len: c_int = args.len() as c_int;

            let map_ret = unsafe { run_pesc_bulk(args_len, arg_ptrs.as_ptr()) };
            if map_ret != 0 {
                bail!("mapper returned exit code {}; failure", map_ret);
            }
        }
    }
    Ok(())
}
